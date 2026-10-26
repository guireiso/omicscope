"""End-to-end checks of every import method on its example input."""
import io
import os

import numpy as np
import pandas as pd
import pytest

import omicscope as omics
from omicscope.General.Input import Input
from omicscope.General.Input.ProteomeDiscoverer import pd_pdata
from omicscope.General.Input.utils import clean_run_name, match_samples

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, 'data', 'importers')
PROT = os.path.join(HERE, 'data', 'proteins')
DECOY = r'(?:^|;)(?:REV__|CON__|rev_|Reverse_|contam_|contaminant_|Cont_|DECOY_)'

# Method, quantitative file, pdata file (None when the software stores the conditions), conditions
CASES = [
    ('MaxQuant', 'MaxQuant_proteinGroups.txt', 'MaxQuant_pdata.xlsx', {'CoV', 'NC', 'VM'}),
    ('MaxQuant', 'MaxQuant_TMT_proteinGroups.txt', 'MaxQuant_TMT_pdata.xlsx', {'Control', 'Treated'}),
    ('DIA-NN', 'DIANN_report.pg_matrix.tsv', 'DIANN_pdata.xlsx', {'Yeast_15ng', 'Yeast_45ng'}),
    ('DIA-NN', 'DIANN_report.tsv', 'DIANN_pdata.xlsx', {'Yeast_15ng', 'Yeast_45ng'}),
    ('Spectronaut', 'Spectronaut_Report_ProteinQuant_Pivot.tsv', 'Spectronaut_pdata.xlsx',
     {'Yeast_15ng', 'Yeast_45ng'}),
    ('Spectronaut', 'Spectronaut_Report_Normal.tsv', None, {'Yeast_15ng', 'Yeast_45ng'}),
    ('FragPipe', 'FragPipe_combined_protein.tsv', 'FragPipe_pdata.xlsx', {'IDHmut', 'IDHwt'}),
    ('FragPipe', 'FragPipe_TMT_abundance_protein_MD.tsv', 'FragPipe_pdata.xlsx', {'IDHmut', 'IDHwt'}),
    ('Progenesis', 'Progenesis_proteins.csv', None, {'COVID', 'CTRL'}),
    ('Progenesis', 'Progenesis_proteins.csv', 'Progenesis_pdata.xlsx', {'COVID', 'CTRL'}),
    ('ProteomeDiscoverer', 'ProteomeDiscoverer_proteins.xlsx', None, {'CTRL', 'CBD (SCZ)'}),
    ('RPKM', 'RPKM_matrix.tsv', 'RPKM_pdata.xlsx', {'Normal', 'Tumor'}),
    ('PatternLab', os.path.join(PROT, 'patternlab.plp'), None, {'CTL', 'CUP+CLZ'}),
    ('General', os.path.join(PROT, 'general.xlsx'), None, {'COVID', 'CTRL'}),
    ('MaxQuant', os.path.join(PROT, 'MQ.txt'), os.path.join(PROT, 'MQ_pdata.xlsx'), {'CoV', 'NC', 'VM'}),
]


def _path(name):
    return name if os.path.isabs(name) else os.path.join(EX, name)


def _run(method, table, pdata=None, **kwargs):
    if pdata is not None:
        kwargs['pdata'] = pdata
    return omics.OmicScope(table, Method=method, ControlGroup=None, **kwargs)


@pytest.mark.parametrize('method,table,pdata,conditions', CASES,
                         ids=[f'{c[0]}-{os.path.basename(c[1])}-{"pdata" if c[2] else "nopdata"}' for c in CASES])
def test_importer_end_to_end(method, table, pdata, conditions):
    data = _run(method, _path(table), None if pdata is None else _path(pdata))
    quant = data.quant_data
    assert set(data.Conditions) == conditions
    assert len(quant) > 50
    # p-values computed (imported statistics may leave a few blanks)
    assert quant['pvalue'].notna().mean() > 0.99
    assert quant['pAdjusted'].notna().mean() > 0.99
    assert quant['pvalue'].dropna().between(0, 1).all()
    # no decoys or contaminants left
    assert not quant['Accession'].astype(str).str.contains(DECOY).any()
    # gene names parsed (one gene per row)
    if method != 'General':  # General keeps the gene names given by the user
        assert quant['gene_name'].notna().all()
    assert not quant['gene_name'].astype(str).str.contains(';').any()
    # every pdata sample was found and kept
    if pdata is not None and method not in ('Progenesis',):
        expected = set(pd.read_excel(_path(pdata))['Sample'].astype(str))
        used = set(','.join(data.pdata['technical']).split(','))
        assert used == expected


@pytest.mark.parametrize('method,table,pdata', [
    ('DIA-NN', 'DIANN_report.pg_matrix.tsv', 'DIANN_pdata.xlsx'),
    ('Spectronaut', 'Spectronaut_Report_ProteinQuant_Pivot.tsv', 'Spectronaut_pdata.xlsx'),
    ('FragPipe', 'FragPipe_combined_protein.tsv', 'FragPipe_pdata.xlsx'),
    ('MaxQuant', 'MaxQuant_proteinGroups.txt', 'MaxQuant_pdata.xlsx'),
])
def test_input_tables_are_consistent(method, table, pdata):
    data = Input(_path(table), Method=method, pdata=_path(pdata))
    assert list(data.assay.index) == list(data.rdata['Accession'])
    assert data.assay.shape[1] > 0
    assert all(np.issubdtype(t, np.number) for t in data.assay.dtypes)
    assert {'Accession', 'gene_name'}.issubset(data.rdata.columns)


def test_uploaded_files_are_read():
    """The app passes file-like objects; pdata is read twice (importer and OmicScope)."""
    table = io.BytesIO(open(_path('DIANN_report.pg_matrix.tsv'), 'rb').read())
    table.name = 'report.pg_matrix.tsv'
    pdata = io.BytesIO(open(_path('DIANN_pdata.xlsx'), 'rb').read())
    data = _run('DIA-NN', table, pdata)
    assert len(data.quant_data) > 50


@pytest.mark.parametrize('rename', [
    lambda s: 'D:\\Raw\\Two-proteome diaPASEF\\' + s + '.d',
    lambda s: ' ' + s.upper() + '.raw ',
    lambda s: '/data/runs/' + s + '.mzML',
])
def test_pdata_sample_names_tolerant(rename):
    pdata = pd.read_excel(_path('DIANN_pdata.xlsx'))
    pdata['Sample'] = pdata['Sample'].map(rename)
    data = _run('DIA-NN', _path('DIANN_report.pg_matrix.tsv'), pdata)
    assert len(data.quant_data) > 50
    assert data.pdata.shape[0] == 6


def test_pdata_sample_names_wrong_lists_available_names():
    pdata = pd.read_excel(_path('DIANN_pdata.xlsx'))
    pdata.loc[0, 'Sample'] = 'sample_typo'
    with pytest.raises(ValueError, match='sample_typo') as error:
        _run('DIA-NN', _path('DIANN_report.pg_matrix.tsv'), pdata)
    assert '200113_AUR_diaPASEF_HY_200ng_15ng_90min_Slot1-5_1_1630' in str(error.value)


@pytest.mark.parametrize('method,table', [
    ('DIA-NN', 'DIANN_report.pg_matrix.tsv'),
    ('MaxQuant', 'MaxQuant_proteinGroups.txt'),
    ('FragPipe', 'FragPipe_combined_protein.tsv'),
    ('RPKM', 'RPKM_matrix.tsv'),
    ('Spectronaut', 'Spectronaut_Report_ProteinQuant_Pivot.tsv'),
])
def test_missing_pdata_gives_clear_error(method, table):
    with pytest.raises(ValueError, match='pdata file is required'):
        _run(method, _path(table))


def test_spectronaut_long_report_conditions():
    data = Input(_path('Spectronaut_Report_Normal.tsv'), Method='Spectronaut')
    assert set(data.pdata['Condition']) == {'Yeast_15ng', 'Yeast_45ng'}
    assert set(data.pdata['Sample']) == set(data.assay.columns)


def test_spectronaut_filtered_values_are_missing():
    data = Input(_path('Spectronaut_Report_ProteinQuant_Pivot.tsv'), Method='Spectronaut',
                 pdata=_path('Spectronaut_pdata.xlsx'))
    raw = pd.read_table(_path('Spectronaut_Report_ProteinQuant_Pivot.tsv'))
    assert (raw == 'Filtered').any().any()
    assert data.assay.isna().any().any()


def test_maxquant_tmt_uses_reporter_intensities():
    data = Input(_path('MaxQuant_TMT_proteinGroups.txt'), Method='MaxQuant',
                 pdata=_path('MaxQuant_TMT_pdata.xlsx'))
    assert data.Params['Params']['ImportWarning_1'] == 'Quantification using Reporter intensity corrected columns'
    assert sorted(data.assay.columns) == sorted(f'{c} Batch{b}' for b in (1, 2) for c in range(1, 7))


def test_fragpipe_tmt_integrator_back_to_linear_scale():
    data = Input(_path('FragPipe_TMT_abundance_protein_MD.tsv'), Method='FragPipe',
                 pdata=_path('FragPipe_pdata.xlsx'))
    assert data.assay.median().min() > 1000


@pytest.mark.parametrize('method,table', [
    ('ProteomeDiscoverer', 'ProteomeDiscoverer_proteins.xlsx'),
    ('PatternLab', os.path.join(PROT, 'patternlab.plp')),
])
def test_external_pdata_for_methods_with_conditions(method, table):
    raw = Input(_path(table), Method=method)
    data = _run(method, _path(table), raw.pdata)
    assert len(data.quant_data) > 50


def test_proteome_discoverer_pdata_parsing():
    lfq = pd_pdata(['Abundances (Normalized): F17: Sample, SCZ_5, CBD (SCZ)',
                    'Abundances (Normalized): F18: Sample, SCZ_5, CBD (SCZ)'])
    assert list(lfq['Condition']) == ['CBD (SCZ)'] * 2
    assert list(lfq['Biological']) == ['SCZ_5'] * 2
    tmt = pd_pdata(['Abundances (Normalized): F1: 126, Sample, CTRL',
                    'Abundances (Normalized): F1: 127N, Sample, CTRL'])
    assert list(tmt['Condition']) == ['CTRL'] * 2
    assert list(tmt['Biological']) == ['F1: 126', 'F1: 127N']


def test_clean_run_name():
    assert clean_run_name('D:\\Raw\\run_01.d.dia') == 'run_01'
    assert clean_run_name('/home/user/run.02.mzML') == 'run.02'
    assert clean_run_name('run_03') == 'run_03'
    assay = pd.DataFrame(columns=['C:\\x\\A1.raw', 'C:\\x\\B1.raw'])
    renamed, _ = match_samples(assay, pd.DataFrame({'Sample': ['A1', 'B1'], 'Condition': ['a', 'b']}))
    assert list(renamed.columns) == ['A1', 'B1']
