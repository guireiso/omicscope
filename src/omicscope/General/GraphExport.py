"""GraphML files for Gephi and Cytoscape.

write_graphml(G, path, pos=None) saves a networkx graph so that it imports cleanly in both tools:
- an XML declaration (UTF-8), the GraphML namespace and `edgedefault` set from the graph;
- one `<key>` per attribute with its real type (string, int, long, double, boolean) and an id that never
  collides with the key ids Gephi maps to node properties (networkx's d3/d7 become Gephi labels);
- attribute names reduced to letters, digits and underscores (so "-log10(pAdj)" becomes
  "minus_log10_pAdj"), the colour and size attributes OmicScope uses (Color/color, Size) written as
  `color` (hex, plus integer `r`, `g`, `b`) and `size`, which Gephi applies on import;
- `name` and `label` with the node name (Cytoscape shows `name`, Gephi shows `label`);
- missing values (None, NaN, inf) left out instead of written as "nan"; lists joined with ';';
- node positions, when given, as `x` and `y`, so Gephi opens the network with the same layout.
"""
import math
import re
import xml.etree.ElementTree as ET

import numpy as np

GRAPHML_NS = 'http://graphml.graphdrawing.org/xmlns'
XSI_NS = 'http://www.w3.org/2001/XMLSchema-instance'
SCHEMA = GRAPHML_NS + ' http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd'
VALID_NAME = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*$')
_INVALID_XML = re.compile('[^\t\n\r -퟿-�\U00010000-\U0010ffff]')
_HEX = re.compile(r'^#[0-9a-f]{6}([0-9a-f]{2})?$')
_INT32 = 2 ** 31
_RENAMED = {'color': 'color', 'size': 'size', 'module': 'module', 'weight': 'weight', 'label': 'source_label',
            'name': 'source_name', 'x': 'source_x', 'y': 'source_y', 'r': 'source_r', 'g': 'source_g', 'b': 'source_b'}


def attribute_name(name):
    """GraphML/Cytoscape/Gephi-safe column name: letters, digits and underscores."""
    text = str(name).strip()
    text = re.sub(r'^-', 'minus_', text)
    text = re.sub(r'[^A-Za-z0-9_]+', '_', text).strip('_')
    if not text or text[0].isdigit():
        text = 'a_' + text
    return _RENAMED.get(text.lower(), text)


def _scalar(value):
    """Python scalar ready for GraphML, or None when the value is missing."""
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        value = float(value)
        return value if math.isfinite(value) else None
    if isinstance(value, str):
        return _INVALID_XML.sub('', value)
    if isinstance(value, (list, tuple, set, np.ndarray)):
        return ';'.join(str(v) for v in value if _scalar(v) is not None)
    try:
        import pandas as pd
        if value is pd.NA or value is pd.NaT:
            return None
    except ImportError:  # pragma: no cover
        pass
    return _INVALID_XML.sub('', str(value))


def _type_of(values):
    kinds = {type(v) for v in values}
    if kinds == {bool}:
        return 'boolean'
    if kinds <= {int}:
        return 'int' if all(-_INT32 <= v < _INT32 for v in values) else 'long'
    if kinds <= {int, float}:
        return 'double'
    return 'string'


def _text(value, kind):
    if kind == 'boolean':
        return 'true' if value else 'false'
    if kind == 'double':
        return repr(float(value))
    if kind in ('int', 'long'):
        return str(int(value))
    return str(value).lower() if isinstance(value, bool) else str(value)


def _clean_attributes(attrs):
    row = {}
    for name, value in attrs.items():
        key = attribute_name(name)
        if key == 'color':
            color = str(value).strip().lower()
            if not _HEX.match(color):
                try:
                    import matplotlib.colors as mcolors
                    color = mcolors.to_hex(value)
                except (ValueError, TypeError):
                    continue
            row['color'] = color[:7]
            row['r'], row['g'], row['b'] = (int(row['color'][i:i + 2], 16) for i in (1, 3, 5))
            continue
        value = _scalar(value)
        if value is None or value == '':
            continue
        row[key] = value
    return row


def graphml_text(nodes, edges, directed=False):
    """GraphML text. nodes: [(id, {name: scalar})]; edges: [(source, target, {name: scalar})];
    attribute names must already be valid (see attribute_name)."""
    clean = lambda v: _INVALID_XML.sub('', str(v))
    nodes = [(clean(v), a) for v, a in nodes]
    edges = [(clean(s), clean(t), a) for s, t, a in edges]
    if len({v for v, _ in nodes}) != len(nodes):
        raise ValueError('duplicate node ids')
    keys = {}
    for domain, rows in (('node', [a for _, a in nodes]), ('edge', [a for _, _, a in edges])):
        seen = {}
        for row in rows:
            for name, value in row.items():
                if not VALID_NAME.match(name):
                    raise ValueError(f'invalid {domain} attribute name for GraphML: {name!r}')
                seen.setdefault(name, []).append(value)
        keys[domain] = {name: _type_of(values) for name, values in seen.items()}
    ET.register_namespace('', GRAPHML_NS)
    ET.register_namespace('xsi', XSI_NS)
    q = lambda tag: f'{{{GRAPHML_NS}}}{tag}'
    root = ET.Element(q('graphml'), {f'{{{XSI_NS}}}schemaLocation': SCHEMA})
    for domain, prefix in (('node', 'n_'), ('edge', 'e_')):
        for name, kind in keys[domain].items():
            ET.SubElement(root, q('key'), {'id': prefix + name, 'for': domain, 'attr.name': name, 'attr.type': kind})
    graph = ET.SubElement(root, q('graph'), {'id': 'G', 'edgedefault': 'directed' if directed else 'undirected'})
    for v, attrs in nodes:
        node = ET.SubElement(graph, q('node'), {'id': v})
        for name, value in attrs.items():
            ET.SubElement(node, q('data'), {'key': 'n_' + name}).text = _text(value, keys['node'][name])
    for i, (s, t, attrs) in enumerate(edges):
        edge = ET.SubElement(graph, q('edge'), {'id': f'e{i}', 'source': s, 'target': t})
        for name, value in attrs.items():
            ET.SubElement(edge, q('data'), {'key': 'e_' + name}).text = _text(value, keys['edge'][name])
    ET.indent(root, space='  ')
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + ET.tostring(root, encoding='unicode') + '\n'


def graph_to_graphml(G, pos=None):
    """GraphML text of a networkx graph (see the module docstring); self-loops are left out."""
    scale, center = 1.0, np.zeros(2)
    if pos:
        xy = np.array([pos[v] for v in G.nodes() if v in pos], dtype=float).reshape(-1, 2)
        if len(xy):
            center = (xy.max(axis=0) + xy.min(axis=0)) / 2
            span = float(np.ptp(xy, axis=0).max())
            scale = 1000.0 / span if span > 0 else 1.0
    nodes = []
    for v, d in G.nodes(data=True):
        attrs = {'name': str(v), 'label': str(v), **_clean_attributes(d)}
        if pos and v in pos:
            attrs['x'] = round(float((pos[v][0] - center[0]) * scale), 3)
            attrs['y'] = round(float((pos[v][1] - center[1]) * scale), 3)
        nodes.append((v, attrs))
    edges = []
    for a, b, d in G.edges(data=True):
        if a == b:
            continue
        attrs = _clean_attributes(d)
        if 'weight' in attrs and not isinstance(attrs['weight'], str):
            attrs['weight'] = float(attrs['weight'])
        edges.append((a, b, attrs))
    return graphml_text(nodes, edges, directed=G.is_directed())


def write_graphml(G, path, pos=None):
    """Save G as GraphML for Gephi and Cytoscape (UTF-8)."""
    with open(path, 'w', encoding='utf-8', newline='\n') as handle:
        handle.write(graph_to_graphml(G, pos))
