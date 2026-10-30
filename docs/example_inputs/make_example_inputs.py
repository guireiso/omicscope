"""Mock example inputs, one per OmicScope import method.

Every file describes the same small experiment, written in the format of each
software, so that users can open a file, see what OmicScope expects and know
what the result should be:

- 2 conditions, Control and Treated, 4 biological replicates each
  (samples Control_1 ... Control_4 and Treated_1 ... Treated_4);
- 200 real human proteins (protein_catalog.csv):
    13 up-regulated in Treated (acute phase / complement, log2FC 3.0 to 1.0),
    13 down-regulated in Treated (mitochondrial respiratory chain, log2FC -2.5 to -1.0),
    IFIT1 detected only in Treated and COX7A1 only in Control,
    172 unchanged;
- a few contaminants and decoys in the formats whose software reports them
  (OmicScope removes them).

The expected answer is written to expected_results.csv. Abundances, peptide
counts, peptide sequences and UniProt entry names (e.g. HP_HUMAN) are
simulated placeholders; accessions, gene names, descriptions and RefSeq IDs
are real.

Run from this folder:  python make_example_inputs.py
"""

import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 2026
CONDITIONS = ['Control', 'Treated']
REPLICATES = 4
SAMPLES = [f'{c}_{i}' for c in CONDITIONS for i in range(1, REPLICATES + 1)]
RUN_FOLDER = 'C:\\Data\\MyExperiment\\'   # runs keep their folder in some software; OmicScope ignores it

# Proteins that search engines report but OmicScope removes
CONTAMINANTS = [('P04264', 'KRT1', 'Keratin, type II cytoskeletal 1 OS=Homo sapiens OX=9606 GN=KRT1 PE=1 SV=6'),
                ('P00761', 'TRYP', 'Trypsin OS=Sus scrofa OX=9823 PE=1 SV=1')]
DECOYS = [('P02768', 'ALB', 'Serum albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2'),
          ('P60709', 'ACTB', 'Actin, cytoplasmic 1 OS=Homo sapiens OX=9606 GN=ACTB PE=1 SV=1')]


def simulate(catalog, rng):
    """log2 abundances (proteins x samples) with the planted changes; NaN = not detected."""
    n = len(catalog)
    base = rng.normal(24.0, 1.8, size=(n, 1))
    log2 = base + rng.normal(0, 0.15, size=(n, len(SAMPLES)))
    treated = np.array([s.startswith('Treated') for s in SAMPLES])
    log2[:, treated] += catalog['planted_log2FC'].fillna(0).to_numpy()[:, None]
    change = catalog['planted_change'].to_numpy()
    log2[change == 'only in Treated', :] = np.where(treated, log2[change == 'only in Treated', :], np.nan)
    log2[change == 'only in Control', :] = np.where(~treated, log2[change == 'only in Control', :], np.nan)
    return pd.DataFrame(log2, index=catalog['Accession'], columns=SAMPLES)


def extra_rows(rng):
    """Contaminant and decoy rows (accession, gene, description, log2 abundances)."""
    rows = []
    for kind, entries in [('contaminant', CONTAMINANTS), ('decoy', DECOYS)]:
        for acc, gene, desc in entries:
            rows.append((kind, acc, gene, desc, rng.normal(27 if kind == 'contaminant' else 19, 0.4, len(SAMPLES))))
    return rows


def statistics(log2):
    """Welch t-test Treated vs Control on log2 abundances, BH-adjusted (for formats that carry statistics)."""
    control, treated = log2.filter(like='Control'), log2.filter(like='Treated')
    with warnings.catch_warnings():       # on/off proteins have no values in one condition: p is NaN
        warnings.simplefilter('ignore')
        p = stats.ttest_ind(treated, control, axis=1, equal_var=False, nan_policy='omit').pvalue
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    adj = np.full_like(p, np.nan)
    adj[ok] = stats.false_discovery_control(p[ok], method='bh')
    fc = treated.mean(axis=1) - control.mean(axis=1)
    return pd.DataFrame({'log2(fc)': fc.to_numpy(), 'pvalue': p, 'pAdjusted': adj}, index=log2.index)


def protein_name(description):
    return description.split(' OS=')[0]


def pdata():
    return pd.DataFrame({'Sample': SAMPLES,
                         'Condition': [s.split('_')[0] for s in SAMPLES],
                         'Biological': [int(s.split('_')[1]) for s in SAMPLES]})


def peptides(rng, n):
    return rng.integers(2, 40, n)


def random_sequence(rng, length):
    return ''.join(rng.choice(list('ACDEFGHIKLMNPQRSTVWY'), length - 1)) + rng.choice(['K', 'R'])


# --- one writer per import method -------------------------------------------------------------

def write_general(catalog, log2, out):
    intensity = 2 ** log2
    rdata = catalog[['Accession', 'gene_name', 'Description']]
    # 1. statistics computed by OmicScope
    with pd.ExcelWriter(os.path.join(out, 'General_perform_statistics.xlsx')) as xl:
        intensity.reset_index(drop=True).to_excel(xl, sheet_name='assay', index=False)
        rdata.to_excel(xl, sheet_name='rdata', index=False)
        pdata().to_excel(xl, sheet_name='pdata', index=False)
    # 2. statistics imported from the user's own analysis (pvalue and pAdjusted columns in rdata)
    st = statistics(log2).reset_index(drop=True)
    with pd.ExcelWriter(os.path.join(out, 'General_import_statistics.xlsx')) as xl:
        intensity.reset_index(drop=True).to_excel(xl, sheet_name='assay', index=False)
        rdata.assign(pvalue=st['pvalue'], pAdjusted=st['pAdjusted']).to_excel(xl, sheet_name='rdata', index=False)
        pdata().to_excel(xl, sheet_name='pdata', index=False)


def write_snapshot(catalog, log2, out):
    st = statistics(log2).reset_index(drop=True)
    table = pd.concat([catalog[['Accession', 'gene_name']], st], axis=1).dropna(subset=['pvalue'])
    with pd.ExcelWriter(os.path.join(out, 'Snapshot.xlsx')) as xl:
        pd.DataFrame([['ControlGroup: Control'], ['Experimental: Treated']]).to_excel(
            xl, index=False, header=False)
        table.to_excel(xl, index=False, startrow=2)


def write_progenesis(catalog, log2, rng, out):
    """Progenesis QI 'Protein measurements' export: 3 header rows, normalized and raw abundance blocks."""
    intensity = (2 ** log2).fillna(0)            # Progenesis writes 0 when a protein is not measured
    st = statistics(log2)
    means = (2 ** log2).T.groupby([s.split('_')[0] for s in SAMPLES]).mean().T
    info = pd.DataFrame({
        'Accession': catalog['Accession'].to_numpy(),
        'Peptide count': peptides(rng, len(catalog)),
        'Unique peptides': peptides(rng, len(catalog)) // 2 + 1,
        'Confidence score': np.round(rng.uniform(20, 400, len(catalog)), 2),
        'Anova (p)': st['pvalue'].to_numpy(),
        'q Value': st['pAdjusted'].to_numpy(),
        'Max fold change': (means.max(axis=1) / means.min(axis=1)).to_numpy(),
        'Power': 1.0,
        'Highest mean condition': means.idxmax(axis=1).to_numpy(),
        'Lowest mean condition': means.idxmin(axis=1).to_numpy(),
        'Description': catalog['Description'].to_numpy()})
    raw = intensity * rng.uniform(0.8, 1.2, size=(1, len(SAMPLES)))
    n_info = info.shape[1]
    header1 = [''] * n_info + ['Normalized abundance'] + [''] * (len(SAMPLES) - 1) + \
              ['Raw abundance'] + [''] * (len(SAMPLES) - 1)
    groups = []
    for c in CONDITIONS:
        groups += [c] + [''] * (REPLICATES - 1)
    header2 = [''] * n_info + groups + groups
    header3 = list(info.columns) + SAMPLES + SAMPLES
    body = pd.concat([info.reset_index(drop=True), intensity.reset_index(drop=True),
                      raw.reset_index(drop=True)], axis=1)
    body.columns = header3
    with open(os.path.join(out, 'Progenesis_proteins.csv'), 'w', newline='', encoding='utf-8') as f:
        for header in (header1, header2):
            f.write(','.join(header) + '\n')
        body.to_csv(f, index=False)


def write_patternlab(catalog, log2, rng, out):
    """PatternLab V .plp: sparse XIC matrix per run, protein index and normalization factors."""
    factor = 1e10
    loci = list(catalog['Accession'])
    descriptions = list(catalog['Description'])
    values = 2 ** log2.to_numpy()
    for kind, acc, gene, desc, row in extra_rows(rng):
        loci.append(('contaminant_' if kind == 'contaminant' else 'Reverse_') + acc)
        descriptions.append(desc)
        values = np.vstack([values, 2 ** row])
    lines = ['###Description : XIC Quant :: mock example (Control vs Treated, 4 replicates each)',
             '###SparseMatrix']
    lines += [f'#ClassDescription\t{i}\t{c}' for i, c in enumerate(CONDITIONS, 1)]
    for j, sample in enumerate(SAMPLES):
        lines.append(f'#{RUN_FOLDER}{sample}.sepr2')
        cls = CONDITIONS.index(sample.split('_')[0]) + 1
        cells = [f'{i}:{v / factor:.10g}' for i, v in enumerate(values[:, j]) if not np.isnan(v)]
        lines.append(f'{cls} ' + ' '.join(cells))
    lines.append('###Index')
    lines += [f'{i}\t{locus}\t{desc}' for i, (locus, desc) in enumerate(zip(loci, descriptions))]
    lines.append('###SecondaryLabels')
    lines += [f'{CONDITIONS.index(s.split("_")[0]) + 1}\t{RUN_FOLDER}{s}.sepr2\t{factor}' for s in SAMPLES]
    with open(os.path.join(out, 'PatternLab.plp'), 'w', newline='\r\n', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


def write_maxquant(catalog, log2, rng, out):
    """MaxQuant proteinGroups.txt with LFQ and raw intensities (0 = not quantified)."""
    rows = [(acc, gene, desc, log2.loc[acc].to_numpy(), '', '')
            for acc, gene, desc in catalog[['Accession', 'gene_name', 'Description']].itertuples(index=False)]
    for kind, acc, gene, desc, row in extra_rows(rng):
        rows.append(('CON__' + acc if kind == 'contaminant' else 'REV__' + acc, gene, desc, row,
                     '+' if kind == 'decoy' else '', '+' if kind == 'contaminant' else ''))
    n = len(rows)
    pep = peptides(rng, n)
    table = pd.DataFrame({
        'Protein IDs': [r[0] for r in rows],
        'Majority protein IDs': [r[0] for r in rows],
        'Peptide counts (all)': pep, 'Peptide counts (razor+unique)': pep, 'Peptide counts (unique)': pep,
        'Protein names': [protein_name(r[2]) for r in rows],
        'Gene names': [r[1] for r in rows],
        'Fasta headers': [f'{r[0].replace("CON__", "").replace("REV__", "")} {r[2]}' for r in rows],
        'Number of proteins': 1, 'Peptides': pep, 'Razor + unique peptides': pep, 'Unique peptides': pep,
        'Sequence coverage [%]': np.round(rng.uniform(5, 60, n), 1),
        'Mol. weight [kDa]': np.round(rng.uniform(10, 200, n), 3),
        'Sequence length': rng.integers(100, 1800, n),
        'Q-value': 0.0,
        'Score': np.round(rng.uniform(10, 300, n), 3)})
    intensity = np.nan_to_num(2 ** np.vstack([r[3] for r in rows]))
    table['Intensity'] = intensity.sum(axis=1)
    for j, s in enumerate(SAMPLES):
        table[f'Intensity {s}'] = np.round(intensity[:, j] * rng.uniform(0.9, 1.1))
    for j, s in enumerate(SAMPLES):
        table[f'LFQ intensity {s}'] = np.round(intensity[:, j])
    table['MS/MS count'] = pep * 3
    table['Only identified by site'] = ''
    table['Reverse'] = [r[4] for r in rows]
    table['Potential contaminant'] = [r[5] for r in rows]
    table['id'] = range(n)
    table.to_csv(os.path.join(out, 'MaxQuant_proteinGroups.txt'), sep='\t', index=False)


def write_diann(catalog, log2, rng, out):
    """DIA-NN report.pg_matrix.tsv (one column per run) and a small main report (report.tsv)."""
    accs, genes, descs, values = (list(catalog['Accession']), list(catalog['gene_name']),
                                  list(catalog['Description']), 2 ** log2.to_numpy())
    for kind, acc, gene, desc, row in extra_rows(rng):
        if kind == 'contaminant':                 # DIA-NN 2.x tags its contaminant library with 'Cont_'
            accs.append('Cont_' + acc), genes.append(gene), descs.append(desc)
            values = np.vstack([values, 2 ** row])
    runs = [f'{RUN_FOLDER}{s}.raw' for s in SAMPLES]
    matrix = pd.DataFrame({'Protein.Group': accs, 'Protein.Ids': accs,
                           'Protein.Names': [f'{g}_HUMAN' for g in genes], 'Genes': genes,
                           'First.Protein.Description': [protein_name(d) for d in descs]})
    matrix = pd.concat([matrix, pd.DataFrame(values, columns=runs)], axis=1)
    matrix.to_csv(os.path.join(out, 'DIANN_report.pg_matrix.tsv'), sep='\t', index=False)

    rows = []
    for i, (acc, gene) in enumerate(zip(accs, genes)):
        for k in range(2):                       # two precursors per protein group
            sequence = random_sequence(rng, int(rng.integers(8, 20)))
            for j, (sample, run) in enumerate(zip(SAMPLES, runs)):
                if np.isnan(values[i, j]):
                    continue
                rows.append({'File.Name': run, 'Run': sample, 'Protein.Group': acc, 'Protein.Ids': acc,
                             'Protein.Names': f'{gene}_HUMAN', 'Genes': gene,
                             'PG.Quantity': values[i, j], 'PG.Normalised': values[i, j],
                             'PG.MaxLFQ': values[i, j], 'Genes.Quantity': values[i, j],
                             'Genes.Normalised': values[i, j], 'Genes.MaxLFQ': values[i, j],
                             'Modified.Sequence': sequence, 'Stripped.Sequence': sequence,
                             'Precursor.Id': f'{sequence}2', 'Precursor.Charge': 2,
                             'Q.Value': 0.001, 'PG.Q.Value': 0.001,
                             'Precursor.Quantity': values[i, j] / 2 * rng.uniform(0.8, 1.2)})
    pd.DataFrame(rows).to_csv(os.path.join(out, 'DIANN_report.tsv'), sep='\t', index=False)


def write_spectronaut(catalog, log2, out):
    """Spectronaut pivot report ('[n] <run>.PG.Quantity' columns) and long report (one row per protein and run)."""
    values = 2 ** log2
    info = pd.DataFrame({'PG.ProteinGroups': catalog['Accession'].to_numpy(),
                         'PG.Genes': catalog['gene_name'].to_numpy(),
                         'PG.ProteinDescriptions': [protein_name(d) for d in catalog['Description']],
                         'PG.Organisms': 'Homo sapiens'})
    quant = values.reset_index(drop=True).astype(object).where(values.notna().to_numpy(), 'Filtered')
    quant.columns = [f'[{i}] {s}.PG.Quantity' for i, s in enumerate(SAMPLES, 1)]
    pd.concat([info, quant], axis=1).to_csv(
        os.path.join(out, 'Spectronaut_Report_Pivot.tsv'), sep='\t', index=False)

    rows = []
    for i, row in info.iterrows():
        for s in SAMPLES:
            v = values.iloc[i][s]
            rows.append({'R.Condition': s.split('_')[0], 'R.FileName': s, 'R.Replicate': int(s.split('_')[1]),
                         **row.to_dict(), 'PG.Quantity': 'Filtered' if np.isnan(v) else v})
    pd.DataFrame(rows).to_csv(os.path.join(out, 'Spectronaut_Report.tsv'), sep='\t', index=False)


def write_fragpipe(catalog, log2, rng, out):
    """FragPipe combined_protein.tsv (IonQuant MaxLFQ; 0 = not quantified)."""
    rows = [(f'sp|{a}|{g}_HUMAN', a, g, d, log2.loc[a].to_numpy())
            for a, g, d in catalog[['Accession', 'gene_name', 'Description']].itertuples(index=False)]
    for kind, acc, gene, desc, row in extra_rows(rng):
        prefix = 'contam_' if kind == 'contaminant' else 'rev_'
        rows.append((f'{prefix}sp|{acc}|{gene}_HUMAN', acc, gene, desc, row))
    n = len(rows)
    spectra = rng.integers(5, 300, n)
    table = pd.DataFrame({
        'Protein': [r[0] for r in rows], 'Protein ID': [r[1] for r in rows],
        'Entry Name': [r[0].split('|')[-1] for r in rows], 'Gene': [r[2] for r in rows],
        'Protein Length': rng.integers(100, 1800, n), 'Organism': 'Homo sapiens OX=9606',
        'Protein Existence': '1:Experimental evidence at protein level',
        'Description': [protein_name(r[3]) for r in rows],
        'Protein Probability': 1.0, 'Top Peptide Probability': 0.999,
        'Combined Total Peptides': peptides(rng, n), 'Combined Spectral Count': spectra,
        'Combined Unique Spectral Count': spectra, 'Combined Total Spectral Count': spectra})
    intensity = np.nan_to_num(2 ** np.vstack([r[4] for r in rows]))
    for j, s in enumerate(SAMPLES):
        table[f'{s} Spectral Count'] = np.where(intensity[:, j] > 0, spectra // 8 + 1, 0)
        table[f'{s} Intensity'] = np.round(intensity[:, j] * rng.uniform(0.9, 1.1))
        table[f'{s} MaxLFQ Intensity'] = np.round(intensity[:, j])
    table['Indistinguishable Proteins'] = ''
    table.to_csv(os.path.join(out, 'FragPipe_combined_protein.tsv'), sep='\t', index=False)


def write_proteome_discoverer(catalog, log2, rng, out):
    """Proteome Discoverer protein table (label-free), exported to Excel."""
    rows = [(a, g, d, log2.loc[a].to_numpy(), False)
            for a, g, d in catalog[['Accession', 'gene_name', 'Description']].itertuples(index=False)]
    rows += [(acc, gene, desc, row, True) for kind, acc, gene, desc, row in extra_rows(rng)
             if kind == 'contaminant']
    n = len(rows)
    table = pd.DataFrame({
        'Checked': False, 'Protein FDR Confidence: Combined': 'High', 'Master': 'Master Protein',
        'Accession': [r[0] for r in rows], 'Description': [r[2] for r in rows],
        'Contaminant': [r[4] for r in rows],
        'Sum PEP Score': np.round(rng.uniform(5, 500, n), 3), 'Coverage [%]': rng.integers(3, 70, n),
        '# Peptides': peptides(rng, n), '# PSMs': rng.integers(5, 400, n),
        '# Unique Peptides': peptides(rng, n), '# AAs': rng.integers(100, 1800, n),
        'MW [kDa]': np.round(rng.uniform(10, 200, n), 1), 'Gene Symbol': [r[1] for r in rows]})
    intensity = 2 ** np.vstack([r[3] for r in rows])
    for j, s in enumerate(SAMPLES):
        table[f'Abundances (Normalized): F{j + 1}: Sample, {s}, {s.split("_")[0]}'] = intensity[:, j]
    table.to_excel(os.path.join(out, 'ProteomeDiscoverer_proteins.xlsx'), index=False)


def write_rpkm(catalog, log2, rng, out):
    """RPKM matrix for transcriptomics: gene_name, Accession (RefSeq) and one column per sample."""
    # transcripts span a lower range than protein intensities; same planted changes
    values = np.round(2 ** (log2 - 16) - 1, 3).clip(lower=0).fillna(0)
    table = pd.concat([pd.DataFrame({'gene_name': catalog['gene_name'].to_numpy(),
                                     'Accession': catalog['RefSeq'].to_numpy()}),
                       values.reset_index(drop=True)], axis=1)
    table.to_csv(os.path.join(out, 'RPKM_matrix.tsv'), sep='\t', index=False)


def main(out=HERE):
    catalog = pd.read_csv(os.path.join(HERE, 'protein_catalog.csv'))
    rng = np.random.default_rng(SEED)
    log2 = simulate(catalog, rng)
    write_general(catalog, log2, out)
    write_snapshot(catalog, log2, out)
    write_progenesis(catalog, log2, rng, out)
    write_patternlab(catalog, log2, rng, out)
    write_maxquant(catalog, log2, rng, out)
    write_diann(catalog, log2, rng, out)
    write_spectronaut(catalog, log2, out)
    write_fragpipe(catalog, log2, rng, out)
    write_proteome_discoverer(catalog, log2, rng, out)
    write_rpkm(catalog, log2, rng, out)
    pdata().to_excel(os.path.join(out, 'pdata.xlsx'), index=False)
    expected = catalog[['Accession', 'gene_name', 'RefSeq', 'planted_log2FC', 'planted_change']]
    expected.to_csv(os.path.join(out, 'expected_results.csv'), index=False)


if __name__ == '__main__':
    main()
