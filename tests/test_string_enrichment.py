"""STRING over-representation (EnrichmentScope(source='STRING'))."""
import os

import pandas as pd
import pytest

import omicscope as omics
from omicscope.EnrichmentAnalysis.STRING import (CATEGORIES, _background_list, category_names,
                                                 format_string_terms)

EX = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'docs', 'example_inputs')

# two rows as STRING's enrichment endpoint returns them
TERMS = [
    {'category': 'KEGG', 'term': 'hsa04610', 'description': 'Complement and coagulation cascades',
     'number_of_genes': 3, 'number_of_genes_in_background': 85, 'inputGenes': ['9606.A', '9606.B', '9606.C'],
     'preferredNames': ['HP', 'FGA', 'C3'], 'p_value': 1e-8, 'fdr': 1e-6},
    {'category': 'PMID', 'term': 'PMID:1', 'description': 'A paper', 'number_of_genes': 2,
     'number_of_genes_in_background': 10, 'inputGenes': ['9606.A', '9606.B'], 'preferredNames': ['HP', 'FGA'],
     'p_value': 1e-3, 'fdr': 0.2},
]


def test_category_names():
    assert category_names(['KEGG', 'Reactome Pathways']) == ['KEGG Pathways', 'Reactome Pathways']
    assert category_names(None) == list(CATEGORIES.values())
    with pytest.raises(ValueError):
        category_names(['Not a category'])


def test_format_like_enrichr():
    df = format_string_terms(TERMS, {'9606.A': 'Hp', '9606.B': 'FGA'}, dbs=None, padjust_cutoff=0.05,
                             foldchange={'HP': 2.0, 'FGA': 1.0, 'C3': -1.0})
    assert len(df) == 1                                   # the PMID row is above the FDR cutoff
    row = df.iloc[0]
    assert row['Gene_set'] == 'KEGG Pathways'
    assert row['Term'] == 'Complement and coagulation cascades (hsa04610)'
    assert row['Overlap'] == '3/85'
    assert row['Genes'] == ['Hp', 'FGA', 'C3']            # the user's names, STRING's name when unmapped
    assert (row['N_Proteins'], row['up-regulated'], row['down-regulated']) == (3, 2, 1)
    assert {'-log10(pAdj)', 'regulation', 'Adjusted P-value', 'Combined Score'} <= set(df.columns)
    assert format_string_terms(TERMS, {}, dbs=['Pfam']).empty


def test_background_needs_genes():
    assert _background_list(None) is None
    assert _background_list(pd.Series(['A', None, 'B'])) == ['A', 'B']
    with pytest.raises(ValueError):
        _background_list(20000)


def test_string_gsea_is_rejected():
    data = omics.OmicScope(os.path.join(EX, 'MaxQuant_proteinGroups.txt'), Method='MaxQuant',
                           pdata=os.path.join(EX, 'pdata.xlsx'), ControlGroup='Control')
    with pytest.raises(ValueError):
        omics.EnrichmentScope(data, Analysis='GSEA', dbs=['KEGG'], organism=9606, source='STRING')


def test_string_finds_planted_biology():
    """Online: the example's up-regulated proteins are complement/coagulation, the down-regulated ones
    oxidative phosphorylation."""
    data = omics.OmicScope(os.path.join(EX, 'MaxQuant_proteinGroups.txt'), Method='MaxQuant',
                           pdata=os.path.join(EX, 'pdata.xlsx'), ControlGroup='Control', FoldChange_cutoff=0.5)
    try:
        enrichment = omics.EnrichmentScope(data, Analysis='ORA', dbs=['KEGG'], organism=9606, source='STRING',
                                           background=True)
    except Exception as error:      # STRING unreachable
        pytest.skip(f'STRING not reachable: {error}')
    res = enrichment.results
    assert enrichment.dbs == ['KEGG Pathways'] and set(res['Gene_set']) == {'KEGG Pathways'}
    complement = res[res['Term'].str.contains('Complement and coagulation')]
    oxphos = res[res['Term'].str.contains('Oxidative phosphorylation')]
    assert len(complement) == 1 and complement['up-regulated'].iloc[0] > complement['down-regulated'].iloc[0]
    assert len(oxphos) == 1 and oxphos['down-regulated'].iloc[0] > oxphos['up-regulated'].iloc[0]
    assert "source='STRING'" in data.reproducible_script(table='x', enrichment=enrichment)
