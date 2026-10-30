
Changelog
=========

2.0.0 (unreleased)
------------------

Statistics

* Moderated (empirical Bayes) t-test and ANOVA by default; Welch/Games-Howell chosen automatically when
  variances differ; vectorised Tukey-Kramer and Games-Howell (much faster).
* Two-stage Benjamini-Krieger-Yekutieli FDR by default (also Benjamini-Hochberg and Storey q-values).
* Presence/absence test for proteins detected in some conditions only, with censored fold changes.
* Effect sizes: mean log2 difference with 95% confidence interval and Hedges' g; eta squared for ANOVA.
* Batch correction (``batch_correction=True``) that keeps the condition effect.
* Longitudinal workflow: correct degrees of freedom, individuals as fixed effects, ``pTime``.
* Every statistical default validated against simulated data with a known answer (``tests/``).

New analyses

* Quality control: missing values, coefficients of variation, sample outliers, batch effects;
  ``exclude_samples`` leaves samples out and recomputes the statistics.
* Protein co-expression map (t-SNE/PCA) and co-expression modules (light WGCNA) with module-trait heatmap.
* PLS-DA with cross-validated accuracy and permutation test; ROC/AUC per protein and cross-validated
  biomarker panels; power analysis.
* GSEA ranked by signed p-value by default, GSEA running-score plot, fold change per term, enrichment per
  cluster, custom gene sets (GMT or dict), redundant-term collapsing.
* Nebula protein meta-analysis (Stouffer, per-study z capped).
* ``methods_text()`` and ``reproducible_script()`` for reporting.

Inputs

* Example inputs for every import method (``docs/example_inputs``): one Control vs Treated experiment with
  planted changes, a shared pdata and the expected results; ``tests/test_example_inputs.py`` checks them.
* RPKM: values were log-transformed twice (log2(RPKM + 1) by the importer, then log2 again by
  ``logTransform``) and fold changes were ratios of log values. The importer now adds 1 and the statistics
  apply log2 once, so fold changes are on the log2(RPKM + 1) scale.
* Normalization and imputation no longer treat an abundance of exactly 1 (log2 = 0) as missing.
* Proteome Discoverer: a condition named like a PD sample type (e.g. "Control") is no longer dropped.
* Snapshot: condition names no longer start with a space.

Other

* Python 3.12+, uv project, pandas 3 / numpy 2 / matplotlib 3.10 compatibility.
* Sample names no longer depend on hash randomisation when pdata has extra columns.
* New logo.

1.0.0 (2022-11-10)
------------------

* First release on PyPI.
