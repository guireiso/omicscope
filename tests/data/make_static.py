"""Build synthetic static (non-longitudinal) datasets in the General format with a known ground truth.

    uv run python tests/data/make_static.py

* static_4groups.xlsx - 4 conditions x 5 biological x 2 technical replicates, 2000 entities; the
  conditions have different variances (Welch / Games-Howell territory).
* static_2groups_n3.xlsx - 2 conditions x 3 biological x 2 technical replicates, 2000 entities; the
  small-n case where variance moderation matters.
* static_onoff.xlsx / static_onoff_3groups.xlsx - 5 biological x 2 technical replicates per condition,
  missing values driven by a limit of detection (LOD); 4% of the entities are switched off in one
  condition (far below the LOD), truth in `simulated_onoff` (the condition where it is absent).

rdata carries the truth: `simulated_de` (any difference among conditions) and, for 4 groups,
`simulated_shift_<condition>` (log2 shift of each condition against the control 'A').
Missing values are abundance dependent (more frequent in low-abundance entities), as in proteomics.
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SEED = 20261022


def build(conditions, n_bio, n_tech, n_entities, de_share, sd_by_condition, rng, rdata_source):
    rdata = (rdata_source.sample(n_entities, replace=n_entities > len(rdata_source), random_state=int(rng.integers(1e9)))
             .reset_index(drop=True))
    rdata['Accession'] = [f'{a}-{i}' if d else a for i, (a, d) in enumerate(zip(rdata['Accession'], rdata['Accession'].duplicated()))]
    baseline = rng.normal(22, 2, size=n_entities)
    de = rng.random(n_entities) < de_share
    shifts = np.zeros((n_entities, len(conditions)))
    for i in np.where(de)[0]:
        # 1-2 conditions move away from the control, by 0.8-2.5 log2 units, up or down
        moved = rng.choice(np.arange(1, len(conditions)), size=rng.integers(1, min(2, len(conditions) - 1) + 1), replace=False)
        shifts[i, moved] = rng.choice([-1, 1], size=len(moved)) * rng.uniform(0.8, 2.5, size=len(moved))
    samples, pdata = {}, []
    bio_id = 0
    for ci, cond in enumerate(conditions):
        for b in range(n_bio):
            bio_id += 1
            bio = rng.normal(0, sd_by_condition[ci], size=n_entities)
            for t in range(n_tech):
                name = f'{cond}_B{b + 1}_T{t + 1}'
                log2 = baseline + shifts[:, ci] + bio + rng.normal(0, 0.1, size=n_entities)
                # abundance-dependent missingness: up to ~25% for the lowest entities
                p_missing = 0.25 / (1 + np.exp((log2 - 19) * 1.5))
                value = np.round(2 ** log2, 2)
                value[rng.random(n_entities) < p_missing] = 0
                samples[name] = value
                pdata.append({'Sample': name, 'Condition': cond, 'Biological': bio_id})
    rd = rdata[['Accession', 'Description', 'Peptide count', 'Unique peptides', 'Confidence score']].copy()
    rd['simulated_de'] = de
    if len(conditions) > 2:
        for ci, cond in enumerate(conditions[1:], start=1):
            rd[f'simulated_shift_{cond}'] = shifts[:, ci].round(3)
    else:
        rd['simulated_log2fc'] = shifts[:, 1].round(3)
    return pd.DataFrame(samples), rd, pd.DataFrame(pdata)


def build_onoff(conditions, n_bio, n_tech, n_entities, de_share, onoff_share, rng, rdata_source, lod=18.5):
    rdata = rdata_source.sample(n_entities, random_state=int(rng.integers(1e9))).reset_index(drop=True)
    baseline = rng.normal(22, 2, size=n_entities)
    shifts = np.zeros((n_entities, len(conditions)))
    kind = rng.choice(['none', 'de', 'onoff'], size=n_entities, p=[1 - de_share - onoff_share, de_share, onoff_share])
    off_in = np.full(n_entities, '', dtype=object)
    for i in np.where(kind == 'de')[0]:
        c = rng.integers(1, len(conditions))
        shifts[i, c] = rng.choice([-1, 1]) * rng.uniform(0.8, 2.5)
    for i in np.where(kind == 'onoff')[0]:
        c = rng.integers(0, len(conditions))            # the condition where the entity is switched off
        shifts[i, c] = -8                               # far below the limit of detection
        baseline[i] = max(baseline[i], lod + 2.5)       # clearly detected where it is present
        off_in[i] = conditions[c]
    samples, pdata, bio_id = {}, [], 0
    for ci, cond in enumerate(conditions):
        for b in range(n_bio):
            bio_id += 1
            bio = rng.normal(0, 0.3, size=n_entities)
            for t in range(n_tech):
                name = f'{cond}_B{b + 1}_T{t + 1}'
                log2 = baseline + shifts[:, ci] + bio + rng.normal(0, 0.1, size=n_entities)
                # censoring around the LOD (soft edge) plus 1% missing completely at random
                p_missing = 1 / (1 + np.exp((log2 - lod) * 3)) + 0.01
                value = np.round(2 ** log2, 2)
                value[rng.random(n_entities) < p_missing] = 0
                samples[name] = value
                pdata.append({'Sample': name, 'Condition': cond, 'Biological': bio_id})
    rd = rdata[['Accession', 'Description', 'Peptide count', 'Unique peptides', 'Confidence score']].copy()
    rd['simulated_de'] = kind != 'none'
    rd['simulated_onoff'] = off_in
    return pd.DataFrame(samples), rd, pd.DataFrame(pdata)


def main():
    rng = np.random.default_rng(SEED)
    source = (pd.read_excel(HERE / 'proteins' / 'general.xlsx', sheet_name='rdata')
              .dropna(subset=['Description']).drop_duplicates('Accession'))
    sets = {
        'static_4groups.xlsx': dict(conditions=['A', 'B', 'C', 'D'], n_bio=5, n_tech=2, n_entities=2000,
                                    de_share=0.2, sd_by_condition=[0.25, 0.25, 0.5, 0.8]),
        'static_2groups_n3.xlsx': dict(conditions=['Control', 'Treated'], n_bio=3, n_tech=2, n_entities=2000,
                                       de_share=0.2, sd_by_condition=[0.35, 0.35]),
    }
    builds = [(name, build, kw) for name, kw in sets.items()] + [
        ('static_onoff.xlsx', build_onoff, dict(conditions=['Control', 'Treated'], n_bio=5, n_tech=2, n_entities=2000,
                                                 de_share=0.15, onoff_share=0.04)),
        ('static_onoff_3groups.xlsx', build_onoff, dict(conditions=['A', 'B', 'C'], n_bio=5, n_tech=2, n_entities=2000,
                                                         de_share=0.15, onoff_share=0.04)),
    ]
    for name, fn, kw in builds:
        assay, rdata, pdata = fn(rng=rng, rdata_source=source, **kw)
        with pd.ExcelWriter(HERE / 'proteins' / name) as writer:
            assay.to_excel(writer, sheet_name='assay', index=False)
            rdata.to_excel(writer, sheet_name='rdata', index=False)
            pdata.to_excel(writer, sheet_name='pdata', index=False)
        print(name, assay.shape, 'DE:', int(rdata['simulated_de'].sum()), 'missing:', f'{(assay == 0).mean().mean():.1%}')


if __name__ == '__main__':
    main()
