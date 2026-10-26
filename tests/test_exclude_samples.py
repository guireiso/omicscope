"""Excluding samples (e.g. flagged outliers) before the analysis."""
import contextlib
import io

import numpy as np
import pytest

import omicscope as omics
from omicscope.General.Quality import sample_outliers
from test_quality import called, run, simulate, write


def _general(**kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return omics.OmicScope(Table='tests/data/proteins/static_2groups_n3.xlsx', Method='General', **kw)


def test_names_from_the_results_and_from_the_input_are_equivalent():
    by_result = _general(exclude_samples=['BioRep_2.Control'])
    by_input = _general(exclude_samples=['Control_B2_T1', 'Control_B2_T2'])      # its two technical replicates
    assert len(by_result.pdata) == len(by_input.pdata) == 5
    a, b = (r.quant_data.set_index('Accession')['pvalue'] for r in (by_result, by_input))
    assert np.allclose(a, b.loc[a.index])
    assert 'Excluded_Samples' in by_result.Params['Params']
    assert "exclude_samples=['BioRep_2.Control']" in by_result.reproducible_script()


def test_clear_errors():
    with pytest.raises(ValueError, match='Samples not found'):
        _general(exclude_samples=['nope'])
    with pytest.raises(ValueError, match='fewer than two samples'):
        _general(exclude_samples=['BioRep_1.Control', 'BioRep_2.Control'])


def test_imported_statistics_are_recomputed_without_the_excluded_samples():
    with contextlib.redirect_stdout(io.StringIO()):
        res = omics.OmicScope(Table='tests/data/proteins/progenesis.xls', Method='Progenesis',
                              exclude_samples=['VCC_HB_B_1.CTRL'])
    assert 'Excluded_Samples_Stats' in res.Params['Params']
    assert 'Stats_Workflow_3' in res.Params['Params']          # OmicScope ran the test itself
    assert len(res.pdata) == 37


def test_removing_a_flagged_outlier_improves_the_analysis(tmp_path):
    assay, rdata, pdata = simulate(1, n=1500, batch_share=0, outlier='noise')
    path = write(tmp_path / 'outlier.xlsx', assay, rdata, pdata.drop(columns='Batch'))
    full = run(path)
    flagged = list(sample_outliers(full).query('outlier')['Sample'])
    assert len(flagged) == 1
    with contextlib.redirect_stdout(io.StringIO()):
        clean = omics.OmicScope(Table=str(path), Method='General', exclude_samples=flagged)
    power_full, fdr_full = called(full, rdata)
    power_clean, fdr_clean = called(clean, rdata)
    assert power_clean > power_full and fdr_clean <= max(fdr_full, 0.06)
