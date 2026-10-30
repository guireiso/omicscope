# Example inputs for every import method

One small experiment, written in the output format of each software that OmicScope reads. Open any
file to see what OmicScope expects, run it, and compare with the known answer.

## The experiment

| | |
|---|---|
| Conditions | **Control** and **Treated** |
| Replicates | 4 per condition: `Control_1` … `Control_4`, `Treated_1` … `Treated_4` |
| Proteins | 200 real human proteins (UniProt accessions, gene names and descriptions; RefSeq IDs for the RPKM file) |
| Up in Treated | 13 acute-phase / complement proteins (HP, FGA, FGB, FGG, C3, CP, HPX, A2M, AGT, AMBP, CLU, C4BPA, C1QB), log2 fold change from 3.0 down to 1.0 |
| Down in Treated | 13 mitochondrial respiratory-chain proteins (ATP5F1A, ATP5F1B, ATP5F1C, ATP5PO, ATP5PB, COX4I1, COX5A, COX5B, COX6B1, COX7C, CYC1, CYCS, ETFB), log2 fold change from −2.5 to −1.0 |
| On/off | **IFIT1** detected only in Treated, **COX7A1** only in Control |
| Unchanged | the other 172 proteins |
| Removed on import | keratin, trypsin and reversed (decoy) entries, in the formats whose software reports them |

`expected_results.csv` lists every protein with its planted change. Abundances, peptide counts, peptide
sequences and UniProt entry names (e.g. `HP_HUMAN`) are simulated placeholders.

## Files

| Import method | File | pdata |
|---|---|---|
| General (OmicScope runs the statistics) | `General_perform_statistics.xlsx` (sheets: assay, rdata, pdata) | inside the file |
| General (statistics imported from the file) | `General_import_statistics.xlsx` (`pvalue`, `pAdjusted` in rdata) | inside the file |
| Snapshot (results only) | `Snapshot.xlsx` | not needed |
| Progenesis QI | `Progenesis_proteins.csv` | not needed (conditions in the header) |
| PatternLab V | `PatternLab.plp` | not needed (classes in the file) |
| MaxQuant | `MaxQuant_proteinGroups.txt` | `pdata.xlsx` |
| DIA-NN | `DIANN_report.pg_matrix.tsv` (recommended) or `DIANN_report.tsv` | `pdata.xlsx` |
| Spectronaut | `Spectronaut_Report_Pivot.tsv` | `pdata.xlsx` |
| Spectronaut | `Spectronaut_Report.tsv` (long report with `R.Condition`) | not needed |
| FragPipe | `FragPipe_combined_protein.tsv` | `pdata.xlsx` |
| Proteome Discoverer | `ProteomeDiscoverer_proteins.xlsx` | not needed (study factors in the column names) |
| RPKM (transcriptomics) | `RPKM_matrix.tsv` | `pdata.xlsx` |

The same `pdata.xlsx` (columns `Sample`, `Condition`, `Biological`) works for every method that needs
one: sample names are matched to the column names of each file, ignoring folders (`C:\Data\...`),
file extensions (`.raw`, `.d`, ...) and letter case.

## What you should see

```python
import omicscope as omics

data = omics.OmicScope('MaxQuant_proteinGroups.txt', Method='MaxQuant',
                       pdata='pdata.xlsx', ControlGroup='Control')
```

- 200 proteins quantified; contaminants and decoys removed.
- All 13 up- and 13 down-regulated proteins significant (adjusted p ≤ 0.05) with the planted direction.
- IFIT1 and COX7A1 labelled *exclusive to Treated* / *exclusive to Control* (presence/absence test).
  In the RPKM file they are 0 RPKM in one condition, i.e. log2(0 + 1) = 0, a measured value, so they come
  out as ordinary large changes.
- With the default settings, about 3 of the 172 unchanged proteins also reach adjusted p ≤ 0.05, with
  |log2 fold change| ≤ 0.33. This is what a 5% false discovery rate allows. With a fold change cutoff of
  0.5 (`FoldChange_cutoff=0.5`, or *Significance → Fold change cutoff* in the app) the list is exactly
  the planted one.
- Enrichment (ORA, Enrichr) of the up-regulated proteins gives *Complement and coagulation cascades*
  (KEGG) and *Complement cascade* (Reactome); of the down-regulated ones, *Respiratory electron transport*
  (Reactome) and *Oxidative phosphorylation* (KEGG).

`tests/test_example_inputs.py` checks all of the above for every file.

## Regenerating

```bash
python make_example_inputs.py
```

The files are rebuilt from `protein_catalog.csv` with a fixed random seed, so the text files come out
identical.
