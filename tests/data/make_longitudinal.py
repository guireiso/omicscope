"""Build synthetic longitudinal datasets (General input format) with a known ground truth.

    uv run python tests/data/make_longitudinal.py

Writes two workbooks next to the other test inputs, one per design of the pdata example:

* longitudinal_independent.xlsx - different subjects at every time point
* longitudinal_related.xlsx     - the same subjects followed over time ('Individual' column)

Each workbook has the three sheets read by Method='General': assay (entities x samples, linear
abundances), rdata (Accession, Description with GN=, peptides, ...) and pdata. The rdata sheet also
carries `simulated_pattern`, the profile every entity was generated from, so tests can check what
the statistics recover.
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SEED = 20261021
N_ENTITIES = 1200

# log2 effects per time point (TimeCourse 1, 3, 5, 7) for the treated group; control stays flat
# unless the pattern is a time-only one (same change in both conditions)
PATTERNS = {
    'flat':               dict(share=0.62, control=[0, 0, 0, 0],      treatment=[0, 0, 0, 0]),
    'time_up':            dict(share=0.06, control=[0, .6, 1.2, 1.8], treatment=[0, .6, 1.2, 1.8]),
    'time_down':          dict(share=0.06, control=[0, -.6, -1.2, -1.8], treatment=[0, -.6, -1.2, -1.8]),
    'treatment_up':       dict(share=0.06, control=[0, 0, 0, 0],      treatment=[0, .8, 1.6, 2.4]),
    'treatment_down':     dict(share=0.06, control=[0, 0, 0, 0],      treatment=[0, -.8, -1.6, -2.4]),
    'treatment_peak':     dict(share=0.05, control=[0, 0, 0, 0],      treatment=[0, 2.0, .8, 0]),
    'treatment_late':     dict(share=0.04, control=[0, 0, 0, 0],      treatment=[0, 0, .3, 2.2]),
    'condition_offset':   dict(share=0.05, control=[0, 0, 0, 0],      treatment=[1.2, 1.2, 1.2, 1.2]),
}
SD_BIOLOGICAL, SD_TECHNICAL, SD_INDIVIDUAL = 0.30, 0.10, 0.45
MISSING_RATE = 0.01


def simulate(pdata, related, rng, rdata):
    times = sorted(pdata['TimeCourse'].unique())
    names = list(PATTERNS)
    shares = np.array([PATTERNS[p]['share'] for p in names])
    pattern = rng.choice(names, size=len(rdata), p=shares / shares.sum())
    baseline = rng.normal(22, 1.8, size=len(rdata))                       # log2 abundance

    bio = pdata[['Biological', 'Condition', 'TimeCourse'] + (['Individual'] if related else [])].drop_duplicates('Biological')
    bio_effect = {b: rng.normal(0, SD_BIOLOGICAL, size=len(rdata)) for b in bio['Biological']}
    ind_effect = ({i: rng.normal(0, SD_INDIVIDUAL, size=len(rdata)) for i in pdata['Individual'].unique()}
                  if related else {})

    columns = {}
    for row in pdata.itertuples(index=False):
        t = times.index(row.TimeCourse)
        effect = np.array([PATTERNS[p]['control' if row.Condition == 'Control' else 'treatment'][t] for p in pattern])
        log2 = baseline + effect + bio_effect[row.Biological] + rng.normal(0, SD_TECHNICAL, size=len(rdata))
        if related:
            log2 = log2 + ind_effect[row.Individual]
        value = np.round(2 ** log2, 2)
        value[rng.random(len(rdata)) < MISSING_RATE] = 0                 # a few missing measurements
        columns[row.Sample] = value
    assay = pd.DataFrame(columns)
    out = rdata.copy()
    out['simulated_pattern'] = pattern
    return assay, out


def main():
    rng = np.random.default_rng(SEED)
    source = pd.read_excel(HERE / 'proteins' / 'general.xlsx', sheet_name='rdata')
    rdata = (source.dropna(subset=['Description'])
             .drop_duplicates('Accession')
             .sample(N_ENTITIES, random_state=SEED)
             [['Accession', 'Description', 'Peptide count', 'Unique peptides', 'Confidence score']]
             .reset_index(drop=True))
    pdata_book = pd.read_excel(HERE / 'proteins' / 'longitudinal_pdata_example.xlsx', sheet_name=None)
    for sheet, related, name in (('Independent', False, 'longitudinal_independent.xlsx'),
                                 ('Related', True, 'longitudinal_related.xlsx')):
        pdata = pdata_book[sheet]
        assay, rd = simulate(pdata, related, rng, rdata)
        with pd.ExcelWriter(HERE / 'proteins' / name) as writer:
            assay.to_excel(writer, sheet_name='assay', index=False)
            rd.to_excel(writer, sheet_name='rdata', index=False)
            pdata.to_excel(writer, sheet_name='pdata', index=False)
        print(name, assay.shape, rd['simulated_pattern'].value_counts().to_dict())


if __name__ == '__main__':
    main()
