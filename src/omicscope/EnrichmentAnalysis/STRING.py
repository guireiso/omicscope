"""Over-representation analysis with the STRING database.

STRING (https://string-db.org) annotates the proteomes of about 12,500 organisms with GO, KEGG, Reactome,
WikiPathways, protein domains, tissues, diseases and more. One request returns every category, so the
choice of categories is a filter on that answer, not a new analysis. Genes are first mapped to STRING
identifiers (one request), which also maps the background and brings the results back to the user's
gene names. Organisms are NCBI taxonomy IDs (9606 = human, 10090 = mouse, ...); `string_species()`
lists them.

Only over-representation is offered: STRING's ranked (GSEA-like) analysis needs an API key and
asynchronous jobs.
"""

import io

import numpy as np
import pandas as pd

from .ORA import ENRICHR_COLUMNS

STRING_API = 'https://version12.string-db.org/api'
SPECIES_URL = 'https://stringdb-downloads.org/download/species.v12.0.txt'
CALLER = 'omicscope'

# STRING category code -> name shown in results and figures
CATEGORIES = {
    'Process': 'GO Biological Process',
    'Function': 'GO Molecular Function',
    'Component': 'GO Cellular Component',
    'KEGG': 'KEGG Pathways',
    'RCTM': 'Reactome Pathways',
    'WikiPathways': 'WikiPathways',
    'Hallmark': 'MSigDB Hallmark',
    'Keyword': 'UniProt Keywords',
    'Pfam': 'Pfam Domains',
    'InterPro': 'InterPro Domains',
    'SMART': 'SMART Domains',
    'COMPARTMENTS': 'Subcellular Localization (COMPARTMENTS)',
    'TISSUES': 'Tissue Expression (TISSUES)',
    'DISEASES': 'Disease-gene Associations (DISEASES)',
    'HPO': 'Human Phenotype (HPO)',
    'GWAS': 'GWAS Catalog',
    'NetworkNeighborAL': 'Local Network Cluster (STRING)',
    'PMID': 'Reference Publications (PubMed)',
}
DEFAULT_CATEGORIES = ['GO Biological Process', 'KEGG Pathways', 'Reactome Pathways']
COMMON_SPECIES = {9606: 'Homo sapiens', 10090: 'Mus musculus', 10116: 'Rattus norvegicus',
                  7955: 'Danio rerio', 7227: 'Drosophila melanogaster', 6239: 'Caenorhabditis elegans',
                  559292: 'Saccharomyces cerevisiae', 3702: 'Arabidopsis thaliana', 9913: 'Bos taurus',
                  9823: 'Sus scrofa', 9031: 'Gallus gallus', 8364: 'Xenopus tropicalis',
                  9544: 'Macaca mulatta', 511145: 'Escherichia coli K-12 MG1655'}


def _post(method, data, timeout=120):
    import requests
    response = requests.post(f'{STRING_API}/json/{method}', data=dict(data, caller_identity=CALLER),
                             timeout=timeout)
    response.raise_for_status()
    return response.json()


def category_names(dbs):
    """Result names for the requested categories (codes or names); None means every category."""
    if dbs is None:
        return list(CATEGORIES.values())
    dbs = [dbs] if isinstance(dbs, str) else list(dbs)
    names = []
    for db in dbs:
        name = CATEGORIES.get(db, db)
        if name not in CATEGORIES.values():
            raise ValueError(f'Unknown STRING category {db!r}. Use one of: {", ".join(CATEGORIES.values())}.')
        names.append(name)
    return names


def string_species():
    """STRING organisms: DataFrame with taxon_id, name and type ('core' organisms have the most complete
    annotation). Falls back to a short list of model organisms when the STRING site is unreachable."""
    import requests
    try:
        text = requests.get(SPECIES_URL, timeout=60).text
        table = pd.read_csv(io.StringIO(text), sep='\t')
        table = table.rename(columns={'#taxon_id': 'taxon_id', 'official_name_NCBI': 'name', 'STRING_type': 'type'})
        return table[['taxon_id', 'name', 'type']].sort_values('name').reset_index(drop=True)
    except Exception:
        return pd.DataFrame({'taxon_id': list(COMMON_SPECIES), 'name': list(COMMON_SPECIES.values()),
                             'type': 'core'})


def map_string_ids(genes, species):
    """Map gene names to STRING identifiers (best match per name): DataFrame query, stringId, preferredName."""
    genes = list(dict.fromkeys(str(g) for g in genes if pd.notna(g) and str(g) != ''))
    if not genes:
        return pd.DataFrame(columns=['query', 'stringId', 'preferredName'])
    mapped = _post('get_string_ids', {'identifiers': '\r'.join(genes), 'species': int(species),
                                      'limit': 1, 'echo_query': 1})
    table = pd.DataFrame(mapped, columns=['queryItem', 'stringId', 'preferredName'])
    return table.rename(columns={'queryItem': 'query'}).drop_duplicates('query')


def string_terms(string_ids, species, background_ids=None):
    """Raw STRING enrichment of STRING identifiers (every category); one request."""
    data = {'identifiers': '\r'.join(string_ids), 'species': int(species)}
    if background_ids is not None:
        data['background_string_identifiers'] = '\r'.join(background_ids)
    return _post('enrichment', data)


def format_string_terms(terms, id_to_gene, dbs=None, padjust_cutoff=0.05, foldchange=None):
    """STRING enrichment rows -> the table EnrichmentScope builds from Enrichr (same columns), keeping
    the requested categories and the user's gene names."""
    foldchange = {} if foldchange is None else foldchange
    keep = set(category_names(dbs))
    rows = []
    for t in terms:
        name = CATEGORIES.get(t['category'], t['category'])
        if name not in keep:
            continue
        genes = [id_to_gene.get(i, p) for i, p in zip(t['inputGenes'], t['preferredNames'])]
        description = t.get('description') or t['term']
        rows.append({'Gene_set': name,
                     'Term': description if description == t['term'] else f"{description} ({t['term']})",
                     'Overlap': f"{t['number_of_genes']}/{t['number_of_genes_in_background']}",
                     'P-value': t['p_value'], 'Adjusted P-value': t['fdr'],
                     'Odds Ratio': np.nan, 'Combined Score': -np.log10(max(t['p_value'], 1e-300)),
                     'Genes': genes})
    df = pd.DataFrame(rows, columns=ENRICHR_COLUMNS)
    df = df[df['Adjusted P-value'] <= padjust_cutoff].copy()
    df['Adjusted P-value'] = df['Adjusted P-value'].clip(lower=1e-300)
    df['-log10(pAdj)'] = -np.log10(df['Adjusted P-value'].astype(float))
    df['N_Proteins'] = df['Genes'].apply(len)
    df['regulation'] = df['Genes'].apply(lambda x: [foldchange.get(str(i).upper(), np.nan) for i in x])
    df['down-regulated'] = df['regulation'].apply(lambda x: len([i for i in x if i < 0]))
    df['up-regulated'] = df['regulation'].apply(lambda x: len([i for i in x if i > 0]))
    df = df.sort_values(['Adjusted P-value', 'Combined Score'], ascending=[True, False])
    return df.reset_index(drop=True).reset_index()


def run_string(genes, species=9606, dbs=None, background=None, padjust_cutoff=0.05, foldchange=None):
    """Over-representation of a gene list with STRING.

    Args:
        genes (list): gene names (or UniProt accessions) to test.
        species (int): NCBI taxonomy ID of the organism (9606 = human). See string_species().
        dbs (list, optional): STRING categories, as codes ('Process', 'KEGG', 'RCTM', ...) or names
            ('GO Biological Process', ...). None keeps every category.
        background (list, optional): background genes (e.g. every quantified protein). None uses the whole
            proteome of the organism.
        padjust_cutoff (float, optional): FDR cutoff. Defaults to 0.05.
        foldchange (dict, optional): upper-case gene name -> log2(fc), for the regulation columns.

    Returns:
        pd.DataFrame: same columns as the Enrichr results.
    """
    background = _background_list(background)
    return StringSession(list(genes) + (background or []), species).enrich(
        genes, dbs=dbs, background=background, padjust_cutoff=padjust_cutoff, foldchange=foldchange)


def _background_list(background):
    """STRING needs the background proteins themselves: a list (or Series) of genes, or None."""
    if background is None:
        return None
    if isinstance(background, (int, float)):
        raise ValueError('STRING needs the background genes, not their number: pass a list of genes or '
                         'background=True (every quantified protein).')
    return [str(g) for g in background if pd.notna(g)]


class StringSession:
    """Maps every gene once (one request) and then runs any number of enrichments on subsets of them,
    e.g. one per cluster, each costing a single request."""

    def __init__(self, genes, species):
        self.species = int(species)
        mapping = map_string_ids(genes, species)
        self.gene_to_id = dict(zip(mapping['query'], mapping['stringId']))
        self.id_to_gene = dict(zip(mapping['stringId'], mapping['query']))

    def _ids(self, genes):
        return list(dict.fromkeys(self.gene_to_id[str(g)] for g in genes if str(g) in self.gene_to_id))

    def enrich(self, genes, dbs=None, background=None, padjust_cutoff=0.05, foldchange=None):
        ids = self._ids(genes)
        background = _background_list(background)
        background_ids = None if background is None else self._ids(background + list(genes))
        terms = string_terms(ids, self.species, background_ids) if ids else []
        return format_string_terms(terms, self.id_to_gene, dbs=dbs, padjust_cutoff=padjust_cutoff,
                                   foldchange=foldchange)
