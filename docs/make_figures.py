"""Regenerate the figures of docs/analyses.rst (run from the repository root: uv run python docs/make_figures.py)."""
import contextlib
import io
import os
import sys
import warnings

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

warnings.filterwarnings('ignore')
sys.path.insert(0, 'tests')
import omicscope as omics  # noqa: E402

OUT = 'docs/figure/analyses/'
os.makedirs(OUT, exist_ok=True)


def load(name, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return omics.OmicScope(Table=f'tests/data/proteins/{name}', Method='General', **kw)


def save(call, name):
    with contextlib.redirect_stdout(io.StringIO()):
        call(save=OUT + name + '_', vector=False, dpi=110)
    plt.close('all')
    produced = [f for f in os.listdir(OUT) if f.startswith(name + '_')]
    for f in produced:
        os.replace(OUT + f, OUT + name + '.png')


two = load('static_2groups_n3.xlsx')
four = load('static_4groups.xlsx')
save(lambda **k: four.missing_values(**k), 'missing_values')
save(lambda **k: four.cv_plot(**k), 'cv_plot')
save(lambda **k: four.outlier_plot(**k), 'outlier_plot')
save(lambda **k: two.effect_plot(top=15, **k), 'effect_plot')
save(lambda **k: four.protein_map(**k), 'protein_map')
save(lambda **k: four.module_trait_plot(**k), 'module_trait_plot')
save(lambda **k: two.plsda_plot(**k), 'plsda_plot')
save(lambda **k: two.roc_plot(**k), 'roc_plot')
save(lambda **k: two.power_plot(**k), 'power_plot')

from test_quality import simulate, write  # noqa: E402
from pathlib import Path  # noqa: E402

path = Path(OUT) / '_batch.xlsx'
write(path, *simulate(5)[:3])
with contextlib.redirect_stdout(io.StringIO()):
    batch = omics.OmicScope(Table=str(path), Method='General')
save(lambda **k: batch.batch_plot(**k), 'batch_plot')
path.unlink()

neb = omics.Nebula('tests/data/MultipleGroups/omics_file')
save(lambda **k: neb.meta_plot(top=20, **k), 'meta_plot')
print(sorted(os.listdir(OUT)))
