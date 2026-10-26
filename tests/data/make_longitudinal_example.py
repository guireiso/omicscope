"""Build the longitudinal example dataset (General input format) shipped with the docs and the app.

    uv run --frozen python tests/data/make_longitudinal_example.py

Design: Control vs Treated, 4 individuals per condition, every individual sampled at 0, 6, 12, 24
and 48 h (repeated measures, 'Individual' column), one injection per sample (40 samples). Protein
identities (Accession, Description with GN=, peptide counts) come from tests/data/proteins/general.xlsx
and each protein keeps the mean abundance it has in that real dataset as its baseline.

On top of individual-to-individual variation and measurement noise, a minority of proteins follow a
planted temporal pattern (log2 effects at 0, 6, 12, 24, 48 h):

* transient_up      Treated rises at 6-12 h and returns to baseline       (condition x time)
* sustained_down    Treated decreases from 6 h on and stays down          (condition x time)
* late_up           Treated increases only at 24-48 h                     (condition x time)
* divergent         Control drifts up while Treated goes down (mirror)    (condition x time)
* shared_time       both conditions change in the same way over time      (time only)
* baseline_offset   Treated is higher at every time point, no time trend  (condition only)
* flat              no change

Values below a detection limit are left blank (missing not at random), plus a few random gaps.

Outputs (tests/data/proteins):

* longitudinal_example.xlsx         sheets assay, rdata, pdata (+ 'truth', ignored by the importer)
* longitudinal_example_assay.xlsx   sheets assay, rdata (use with the separate pdata file)
* longitudinal_example_pdata.xlsx   the pdata on its own

The app ships a copy of longitudinal_example.xlsx in OmicScope_Streamlit/ExampleFiles.
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / 'proteins'
SEED = 20261022
N_PROTEINS = 2000
TIMES = [0, 6, 12, 24, 48]                      # hours
CONDITIONS = ['Control', 'Treated']
N_INDIVIDUALS = 4                               # per condition

# number of proteins per pattern and log2 effect at each time point (Control, Treated)
PATTERNS = {
    'transient_up':    (90,  [0, 0, 0, 0, 0],            [0, 1.4, 1.1, 0.3, 0]),
    'sustained_down':  (100, [0, 0, 0, 0, 0],            [0, -0.7, -1.1, -1.2, -1.2]),
    'late_up':         (80,  [0, 0, 0, 0, 0],            [0, 0, 0.1, 0.7, 1.5]),
    'divergent':       (60,  [0, 0.4, 0.7, 0.9, 1.0],    [0, -0.4, -0.7, -0.9, -1.0]),
    'shared_time':     (90,  [0, 0.4, 0.8, 1.0, 1.1],    [0, 0.4, 0.8, 1.0, 1.1]),
    'baseline_offset': (50,  [0, 0, 0, 0, 0],            [0.8, 0.8, 0.8, 0.8, 0.8]),
}
INTERACTION = ['transient_up', 'sustained_down', 'late_up', 'divergent']   # what pAdjusted should find
# pTime tests the time trend averaged over the conditions: the mirrored 'divergent' profile averages to flat
TIME_EFFECT = ['transient_up', 'sustained_down', 'late_up', 'shared_time']

SD_INDIVIDUAL = 0.45        # log2, person-to-person differences (constant over time)
SD_RESIDUAL = 0.25          # log2, biological + technical noise of each sample
LOD_QUANTILE = 0.05         # detection limit placed at this quantile of all log2 values
LOD_WIDTH = 0.35            # softness of the detection limit (log2)
RANDOM_MISSING = 0.005


def proteins(rng):
    """Real protein annotations and their baseline log2 abundance from general.xlsx."""
    book = pd.read_excel(OUT / 'general.xlsx', sheet_name=None)
    assay, rdata = book['assay'], book['rdata']
    baseline = np.log2(assay.replace(0, np.nan)).mean(axis=1)
    contaminants = set(pd.read_csv(HERE.parents[1] / 'src' / 'omicscope' / 'General' / 'Stats' / 'contaminants.csv')
                       ['Accession'])
    keep = (rdata['Description'].str.contains('GN=', regex=False, na=False)
            & baseline.notna()
            & ~rdata['Accession'].str.split(';').apply(lambda ids: any(a in contaminants for a in ids))
            & ~rdata['Accession'].duplicated())
    rdata = rdata.loc[keep, ['Accession', 'Description', 'Peptide count', 'Unique peptides', 'Confidence score']]
    baseline = baseline[keep]
    pick = np.sort(rng.choice(len(rdata), size=N_PROTEINS, replace=False))
    return rdata.iloc[pick].reset_index(drop=True), baseline.iloc[pick].to_numpy()


def make_pdata():
    rows, bio = [], 0
    for c, cond in enumerate(CONDITIONS):
        for i in range(N_INDIVIDUALS):
            individual = c * N_INDIVIDUALS + i + 1        # unique across conditions
            for t in TIMES:
                bio += 1
                rows.append({'Sample': f'{cond}_P{individual}_{t:02d}h', 'Condition': cond,
                             'Biological': bio, 'TimeCourse': t, 'Individual': individual})
    return pd.DataFrame(rows)


def simulate(rng):
    rdata, baseline = proteins(rng)
    pattern = np.array(['flat'] * N_PROTEINS, dtype=object)
    order = rng.permutation(N_PROTEINS)
    start = 0
    for name, (n, _, _) in PATTERNS.items():
        pattern[order[start:start + n]] = name
        start += n
    pdata = make_pdata()
    individual = {i: rng.normal(0, SD_INDIVIDUAL, N_PROTEINS) for i in pdata['Individual'].unique()}
    log2 = {}
    for row in pdata.itertuples(index=False):
        t = TIMES.index(row.TimeCourse)
        k = CONDITIONS.index(row.Condition) + 1
        effect = np.array([PATTERNS[p][k][t] if p in PATTERNS else 0.0 for p in pattern])
        log2[row.Sample] = (baseline + effect + individual[row.Individual]
                            + rng.normal(0, SD_RESIDUAL, N_PROTEINS))
    log2 = pd.DataFrame(log2)
    lod = np.quantile(log2.to_numpy(), LOD_QUANTILE)
    p_missing = 1 / (1 + np.exp((log2 - lod) / LOD_WIDTH))
    missing = (rng.random(log2.shape) < p_missing) | (rng.random(log2.shape) < RANDOM_MISSING)
    assay = np.round(2 ** log2).mask(missing).astype('Int64')
    truth = pd.DataFrame({'Accession': rdata['Accession'],
                          'gene_name': rdata['Description'].str.split('GN=').str[1].str.split(' ').str[0],
                          'pattern': pattern,
                          'condition_x_time': np.isin(pattern, INTERACTION),
                          'time_effect': np.isin(pattern, TIME_EFFECT)})
    return assay, rdata, pdata, truth


def main():
    rng = np.random.default_rng(SEED)
    assay, rdata, pdata, truth = simulate(rng)
    with pd.ExcelWriter(OUT / 'longitudinal_example.xlsx') as writer:
        assay.to_excel(writer, sheet_name='assay', index=False)
        rdata.to_excel(writer, sheet_name='rdata', index=False)
        pdata.to_excel(writer, sheet_name='pdata', index=False)
        truth.to_excel(writer, sheet_name='truth', index=False)
    with pd.ExcelWriter(OUT / 'longitudinal_example_assay.xlsx') as writer:
        assay.to_excel(writer, sheet_name='assay', index=False)
        rdata.to_excel(writer, sheet_name='rdata', index=False)
    pdata.to_excel(OUT / 'longitudinal_example_pdata.xlsx', sheet_name='pdata', index=False)
    print(assay.shape, f'{assay.isna().mean().mean():.1%} missing', truth['pattern'].value_counts().to_dict())


if __name__ == '__main__':
    main()
