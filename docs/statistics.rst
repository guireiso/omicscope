=====================
Statistical workflow
=====================

This page describes what OmicScope does between reading your table and reporting differentially
abundant proteins, which defaults it uses and why. Every choice is recorded in ``object.Params``, and
``object.methods_text()`` turns it into a methods paragraph with references.

Each default was checked against simulated data sets with a known answer (``tests/data``). The numbers
quoted below come from those tests.

Preprocessing
=============

1. **Technical replicates** of the same biological sample are averaged. Missing technical values are
   ignored rather than counted as zero.
2. **Filtering.** A protein is tested on intensity when it is quantified in at least two biological
   replicates of every condition. Proteins quantified in some conditions only go to the presence/absence
   layer (below) or are listed in ``object.excluded_values``.
3. **Contaminants** (Frankenfield et al., 2022) are removed when ``ExcludeContaminants=True``.
4. **Normalization** (``normalization_method``): ``'median'``, ``'average'`` or ``'quantile'``. Check the
   result with ``object.normalization_boxplot()``.
5. **Imputation** (``imputation_method``): ``'mean'``, ``'median'``, ``'knn'`` or ``'minprob'``. Missing values
   in proteomics are mostly below the detection limit; ``object.missing_values()`` shows whether that is
   the case in your data. If it is, prefer ``'minprob'`` or the presence/absence layer.
6. **log2 transformation** (``logTransform=True``).
7. **Batch correction** (``batch_correction=True``, column ``Batch`` in pdata): batch terms are removed with
   a linear model that keeps the condition effect. A design where batch and condition are confounded
   cannot be corrected and raises an error. ``object.batch_plot()`` shows the samples before and after.

Two conditions
==============

``moderated=True`` (default)
   t-test with empirical Bayes variance moderation (limma's ``squeezeVar``, Smyth 2004). Each protein's
   variance is shrunk towards a prior estimated from all proteins, which stabilises tests with few
   replicates. With three replicates per group it found 96% of the true changes, against 70% for the
   ordinary t-test, at the same false discovery rate.

``welch='auto'`` (default)
   Welch's t-test is used instead when the typical variance differs more than two-fold between
   conditions and every condition has at least four replicates. ``welch=True`` or ``False`` forces the
   choice. Welch tests are not moderated.

``independent_ttest=False``
   Paired t-test for related samples (same individual in both conditions).

Effect sizes
   ``quant_data`` gains ``Mean difference`` (difference of log2 means), ``CI95 low``/``CI95 high`` (its
   95% confidence interval, using the same standard error as the test) and ``Hedges' g``. In simulations
   the moderated interval covered the true fold change 95.4% of the time while being 23% narrower than
   the ordinary t interval. ``object.effect_plot()`` draws them as a forest plot.

More than two conditions
========================

One-way ANOVA (moderated by default, Welch ANOVA under the same rule as above), followed by post-hoc
tests on the proteins significant after correction: Tukey-Kramer, or Games-Howell when Welch tests are
used. ``pTukey`` holds the post-hoc p-value of the pair that differs most. ``Eta squared`` gives the share
of each protein's variance explained by the conditions.

Multiple testing
================

``fdr='fdr_tsbh'`` (default)
   Two-stage Benjamini-Krieger-Yekutieli procedure. It estimates the share of unchanged proteins and is
   less conservative than Benjamini-Hochberg when many proteins change, while controlling the FDR.

``fdr='fdr_bh'``
   Benjamini-Hochberg.

``fdr='qvalue'``
   Storey's q-values.

Proteins detected in some conditions only
=========================================

A protein quantified in one condition and absent from another cannot be tested on intensity, but its
absence is information. With ``detection_layer=True`` (default) such proteins, when detected in at least
``min_detected`` replicates of some condition, are tested on detection frequency with an exact
Fisher-Freeman-Halton test. Their FDR is controlled separately, they are labelled
``exclusive to <condition>`` in ``Regulation``, and their fold change is bounded by the detection limit
(``fc_censored``). In simulations at least 90% of the on/off proteins were recovered with an empirical
FDR below 7%.

Longitudinal designs
====================

``ExperimentalDesign='longitudinal'`` fits natural cubic splines over time (``degrees_of_freedom``,
default 2) and compares a model with condition-specific trajectories against a model with a common
trajectory (Storey et al., 2005). When the same individuals are measured repeatedly (``Individual``
column in pdata), individuals enter the model as fixed effects. ``pTime`` tests whether a protein changes
over time at all. In simulations the false discovery rate of the default workflow stayed at the
nominal level (3% at FDR 5%).

Imported statistics
===================

If your table already has a ``pvalue`` column (``General`` method) or you use ``Snapshot``, OmicScope keeps
your statistics and only computes what is missing.
