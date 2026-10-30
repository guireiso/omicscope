"""The mock example inputs (docs/example_inputs) give the planted answer with every import method."""
import filecmp
import os
import sys

import pandas as pd
import pytest

import omicscope as omics
from omicscope.General.Snapshot import Omicscope_Snapshot

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, '..', 'docs', 'example_inputs')
PDATA = os.path.join(EX, 'pdata.xlsx')
EXPECTED = pd.read_csv(os.path.join(EX, 'expected_results.csv'))
UP = set(EXPECTED.loc[EXPECTED['planted_log2FC'] > 0, 'gene_name'])
DOWN = set(EXPECTED.loc[EXPECTED['planted_log2FC'] < 0, 'gene_name'])
ON_OFF = {'IFIT1': 'exclusive to Treated', 'COX7A1': 'exclusive to Control'}
REMOVED = {'KRT1', 'TRYP', 'ALB', 'ACTB'}       # contaminants and decoys written to some files

# method, file, needs pdata, statistics computed by OmicScope (False = imported from the file)
CASES = [
    ('General', 'General_perform_statistics.xlsx', False, True),
    ('General', 'General_import_statistics.xlsx', False, False),
    ('Progenesis', 'Progenesis_proteins.csv', False, False),
    ('PatternLab', 'PatternLab.plp', False, True),
    ('MaxQuant', 'MaxQuant_proteinGroups.txt', True, True),
    ('DIA-NN', 'DIANN_report.pg_matrix.tsv', True, True),
    ('DIA-NN', 'DIANN_report.tsv', True, True),
    ('Spectronaut', 'Spectronaut_Report_Pivot.tsv', True, True),
    ('Spectronaut', 'Spectronaut_Report.tsv', False, True),
    ('FragPipe', 'FragPipe_combined_protein.tsv', True, True),
    ('ProteomeDiscoverer', 'ProteomeDiscoverer_proteins.xlsx', False, True),
    ('RPKM', 'RPKM_matrix.tsv', True, True),
]


@pytest.mark.parametrize('method,name,needs_pdata,own_stats', CASES, ids=[c[1] for c in CASES])
def test_example_input_gives_planted_answer(method, name, needs_pdata, own_stats):
    kwargs = {'pdata': PDATA} if needs_pdata else {}
    data = omics.OmicScope(os.path.join(EX, name), Method=method, ControlGroup='Control',
                           FoldChange_cutoff=0.5, **kwargs)
    quant = data.quant_data
    assert set(data.Conditions) == {'Control', 'Treated'}
    assert len(quant) == 200
    assert not REMOVED & set(quant['gene_name'])
    deps = data.deps
    found_up = set(deps.loc[deps['log2(fc)'] > 0, 'gene_name']) - set(ON_OFF)
    found_down = set(deps.loc[deps['log2(fc)'] < 0, 'gene_name']) - set(ON_OFF)
    # with a fold-change cutoff of 0.5 the result is exactly the planted list
    assert found_up == UP
    assert found_down == DOWN
    if method == 'RPKM':
        # a transcript at 0 RPKM is log2(0 + 1) = 0, a measured value: on/off genes are ordinary changes
        assert {'IFIT1'} <= set(deps.loc[deps['log2(fc)'] > 0, 'gene_name'])
        assert {'COX7A1'} <= set(deps.loc[deps['log2(fc)'] < 0, 'gene_name'])
    elif own_stats:
        labels = quant.set_index('gene_name').loc[list(ON_OFF), 'Regulation'].to_dict()
        assert labels == ON_OFF


def test_snapshot_example():
    data = Omicscope_Snapshot(os.path.join(EX, 'Snapshot.xlsx'))
    assert data.ControlGroup == 'Control'
    assert data.experimental == ['Treated']
    deps = data.deps
    assert set(deps.loc[deps['log2(fc)'] > 0, 'gene_name']) == UP
    assert set(deps.loc[deps['log2(fc)'] < 0, 'gene_name']) == DOWN


def test_example_inputs_are_reproducible(tmp_path):
    sys.path.insert(0, EX)
    try:
        import make_example_inputs
    finally:
        sys.path.remove(EX)
    make_example_inputs.main(str(tmp_path))
    text_files = [f for f in os.listdir(tmp_path) if f.endswith(('.csv', '.tsv', '.txt', '.plp'))]
    assert len(text_files) >= 9
    for f in text_files:        # Excel files carry a timestamp, so only text files are compared
        assert filecmp.cmp(tmp_path / f, os.path.join(EX, f), shallow=False), f
