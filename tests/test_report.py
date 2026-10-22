"""Methods paragraph and reproducible script."""
import contextlib
import io

import numpy as np

import omicscope as omics


def _run(**kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return omics.OmicScope(Table='tests/data/proteins/static_4groups.xlsx', Method='General', **kw)


def test_methods_text_describes_what_was_run():
    text = _run().methods_text()          # simulated conditions differ ~10-fold in variance: Welch chosen
    assert 'Welch ANOVA' in text and 'Games-Howell' in text and 'chosen automatically' in text
    assert 'Benjamini-Krieger-Yekutieli' in text and 'Reis-de-Oliveira' in text
    assert 'detection frequency' not in text        # no on/off proteins in this data set
    text = _run(welch=False).methods_text()
    assert 'moderated' in text and 'Tukey' in text and 'Smyth GK' in text
    text = _run(fdr='fdr_bh', moderated=False, welch=False).methods_text()
    assert 'Benjamini-Hochberg procedure' in text and 'Smyth GK' not in text


def test_reproducible_script_gives_the_same_results():
    data = _run(fdr='qvalue')
    script = data.reproducible_script(table='tests/data/proteins/static_4groups.xlsx')
    namespace = {}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(script.replace("data.quant_data.to_csv('quantitative_results.csv', index=False)", ''), namespace)
    again = namespace['data'].quant_data.set_index('Accession')
    first = data.quant_data.set_index('Accession')
    assert np.allclose(first['pAdjusted'], again.loc[first.index, 'pAdjusted'])
