"""GraphML files saved by the network functions: valid, typed and readable by Gephi and Cytoscape."""
import io
import json
import math
import re
import xml.etree.ElementTree as ET

import matplotlib
import networkx as nx
import numpy as np
import pytest

import omicscope as omics
from omicscope.General.GraphExport import graph_to_graphml

matplotlib.use('Agg')

NS = '{http://graphml.graphdrawing.org/xmlns}'
NAME = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')
HEX = re.compile(r'^#[0-9a-f]{6}$')
GEPHI_RESERVED_IDS = {'d3', 'd7', 'label', 'x', 'y', 'z', 'size', 'r', 'g', 'b', 'color', 'weight', 'id', 'edgeid'}


def check(path_or_text, layout=True):
    text = open(path_or_text, encoding='utf-8').read() if not path_or_text.startswith('<?xml') else path_or_text
    assert text.startswith('<?xml version="1.0" encoding="UTF-8"?>')
    root = ET.fromstring(text.encode('utf-8'))
    keys = {}
    for k in root.findall(NS + 'key'):
        assert k.get('id') not in keys and k.get('id') not in GEPHI_RESERVED_IDS
        assert k.get('for') in ('node', 'edge') and NAME.match(k.get('attr.name'))
        assert k.get('attr.type') in ('string', 'int', 'long', 'double', 'boolean')
        keys[k.get('id')] = (k.get('for'), k.get('attr.name'), k.get('attr.type'))
    graph = root.find(NS + 'graph')
    assert graph.get('edgedefault') == 'undirected'
    nodes = {}
    for element in graph.iter():
        if element.tag not in (NS + 'node', NS + 'edge'):
            continue
        attrs = {}
        for d in element.findall(NS + 'data'):
            domain, name, kind = keys[d.get('key')]
            assert element.tag == NS + domain
            value = d.text or ''
            assert value not in ('', 'nan', 'NaN', 'None', 'inf', '-inf') and not value.startswith('[')
            if kind == 'double':
                assert math.isfinite(float(value))
            elif kind in ('int', 'long'):
                int(value)
            elif kind == 'boolean':
                assert value in ('true', 'false')
            attrs[name] = value
        if element.tag == NS + 'node':
            assert attrs['name'] == attrs['label'] == element.get('id')
            assert (('x' in attrs) and ('y' in attrs)) == layout
            if 'color' in attrs:
                assert HEX.match(attrs['color']) and {'r', 'g', 'b'} <= set(attrs)
            nodes[element.get('id')] = attrs
        else:
            assert element.get('source') in nodes and element.get('target') in nodes
            assert 'weight' not in attrs or keys['e_weight'][2] == 'double'
    G = nx.read_graphml(io.BytesIO(text.encode('utf-8')))
    assert not G.is_directed() and G.number_of_nodes() == len(nodes)
    return G, {name: kind for _, name, kind in keys.values()}


def test_graph_to_graphml_cleans_attributes():
    G = nx.Graph()
    G.add_node('A', Color='#E4001B', Size=np.float64(12.5), Module=np.int64(2), **{'-log10(pAdj)': 3.2})
    G.add_node('B', Color='red', Size=np.nan, genes=['X', 'Y'], flag=np.bool_(True), label='old')
    G.add_edge('A', 'B', weight=np.int64(1))
    G.add_edge('A', 'A', weight=1.0)
    G2, kinds = check(graph_to_graphml(G, pos={'A': (0, 0), 'B': (1, 2)}))
    assert kinds['color'] == 'string' and kinds['size'] == 'double' and kinds['module'] == 'int'
    assert kinds['minus_log10_pAdj'] == 'double' and kinds['flag'] == 'boolean' and kinds['weight'] == 'double'
    assert G2.nodes['B']['color'] == '#ff0000' and G2.nodes['B']['genes'] == 'X;Y' and 'size' not in G2.nodes['B']
    assert G2.nodes['B']['source_label'] == 'old' and G2.number_of_edges() == 1       # self-loop left out
    assert G2.nodes['B']['y'] > G2.nodes['A']['y']                                   # y keeps pointing up


@pytest.fixture(scope='module')
def progenesis():
    return omics.OmicScope(Table='tests/data/proteins/progenesis.xls', Method='Progenesis')


@pytest.fixture(scope='module')
def enrichment(progenesis):
    q = progenesis.quant_data.dropna(subset=['gene_name'])
    regulated = list(q[q[progenesis.pvalue] <= 0.05].sort_values(progenesis.pvalue)['gene_name'].drop_duplicates())
    other = list(q[q[progenesis.pvalue] > 0.5]['gene_name'].drop_duplicates())
    gene_sets = {f'Pathway {i}': regulated[i * 8:i * 8 + 20] + other[i * 3:i * 3 + 4] for i in range(10)}
    return omics.EnrichmentScope(progenesis, Analysis='ORA', dbs=[gene_sets], background=True)


def test_enrichment_graphml(enrichment, tmp_path):
    prefix = str(tmp_path) + '/'
    enrichment.enrichment_map(save=prefix)
    G, kinds = check(prefix + 'PathMap.graphml')
    assert kinds['module'] == 'int' and kinds['size'] == 'double'
    enrichment.enrichment_network(save=prefix)
    G, kinds = check(prefix + 'PPNetwork_custom.graphml')
    assert {d['node_type'] for _, d in G.nodes(data=True)} == {'pathway', 'protein'}


def test_ppi_graphml(progenesis, tmp_path, monkeypatch):
    from omicscope.General import GeneralVisualization as gv
    genes = list(progenesis.quant_data.sort_values(progenesis.pvalue)['gene_name'].dropna().drop_duplicates()[:40])
    rng = np.random.default_rng(0)
    rows = [{'preferredName_A': a, 'preferredName_B': b, 'score': float(rng.uniform(0.6, 1))}
            for i, a in enumerate(genes) for b in genes[i + 1:] if rng.random() < 0.1]

    class Response:
        text = json.dumps(rows)
    monkeypatch.setattr(gv.requests, 'post', lambda *a, **k: Response())
    progenesis.PPInteractions(save=str(tmp_path) + '/', modules=True)
    G, kinds = check(str(tmp_path) + '/PPinteraction.graphml')
    assert kinds['weight'] == 'double' and kinds['log2FC'] == 'double'


def test_nebula_graphml(tmp_path):
    neb = omics.Nebula('tests/data/MultipleGroups/omics_file')
    prefix = str(tmp_path) + '/'
    neb.whole_network(save=prefix)
    assert all('Experiment' not in df.columns for df in neb.original)
    G, kinds = check(prefix + 'PPNetwork.graphml')
    assert {d['node_type'] for _, d in G.nodes(data=True)} == {'experiment', 'protein'}
    neb.similarity_network(save=prefix)
    check(prefix + 'Similarity_network.graphml')
    neb.stat_network(save=prefix)
    G, kinds = check(prefix + 'groupNetwork.graphml')
    assert all(d['weight'] >= 0 for *_, d in G.edges(data=True))
