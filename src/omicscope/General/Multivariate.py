""" Supervised multivariate analysis and biomarker evaluation of an OmicScope object.

plsda: partial least squares discriminant analysis (PLS-DA) of the samples, with variable importance in
    projection (VIP) per protein. A supervised projection always separates the groups of the samples it
    was fitted on, even for pure noise, so the separation is judged by cross-validated accuracy and a
    label-permutation p-value, never by the score plot alone.
roc_table: per-protein receiver operating characteristic (ROC) area under the curve (AUC) between two
    conditions, with a 95% confidence interval.
biomarker_panel: cross-validated logistic-regression panel of a few proteins; feature selection runs
    inside each training fold, so the reported AUC is not inflated by information leakage.

Figures: plsda_plot (scores + top VIP proteins) and roc_plot (single proteins + cross-validated panel).

Missing values are treated as below the detection limit and filled, per protein, with the lowest value
of that protein (as in Exploration.protein_embedding).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .Exploration import _condition_of, _log_matrix, _sample_matrix


def _filled_log(self, data):
    """log2 matrix (proteins x samples) with missing values set to each protein's minimum."""
    mat, genes = _sample_matrix(self, data)
    logm = _log_matrix(self, mat)
    logm = logm.apply(lambda r: r.fillna(r.min()), axis=1)   # missing = below detection
    logm = logm.dropna()
    return logm[logm.std(axis=1) > 0], genes


def _two_groups(self, control, case):
    """log2 matrix (samples x proteins) of two conditions and y (1 = case)."""
    control = self.Conditions[0] if control is None else control
    case = self.Conditions[1] if case is None else case
    if control == case or control not in self.Conditions or case not in self.Conditions:
        raise ValueError(f'control and case must be two different conditions among {self.Conditions}')
    cond = _condition_of(self)
    mat, genes = _sample_matrix(self, self.quant_data)
    samples = [s for s in mat.columns if cond[s] in (control, case)]
    logm, genes = _filled_log(self, self.quant_data[['Accession', 'gene_name'] + samples])
    y = np.array([cond[s] == case for s in logm.columns], dtype=int)
    return logm.T, y, genes, control, case


# ---------------------------------------------------------------------------- PLS-DA


def _vip(pls):
    """Variable importance in projection (Wold et al. 2001; Chong & Jun 2005)."""
    t, w, q = pls.x_scores_, pls.x_weights_, pls.y_loadings_
    ss = np.sum(q ** 2, axis=0) * np.sum(t ** 2, axis=0)       # y variance explained by each component
    wn = w / np.linalg.norm(w, axis=0)
    return np.sqrt(w.shape[0] * (wn ** 2 @ ss) / ss.sum())


def _plsda_cv_accuracy(X, y, n_components, splits):
    """Accuracy of out-of-fold class predictions (argmax of the predicted one-hot matrix)."""
    from sklearn.cross_decomposition import PLSRegression

    Y = np.eye(y.max() + 1)[y]
    pred = np.empty_like(y)
    for train, test in splits:
        k = min(n_components, len(train) - 1)
        pls = PLSRegression(n_components=k, scale=True).fit(X[train], Y[train])   # scaling learnt on train only
        pred[test] = pls.predict(X[test]).argmax(axis=1)
    return float(np.mean(pred == y))


def plsda(self, pvalue=1.0, n_components=2, cv_folds=5, permutations=100, max_features=5000, random_state=0):
    """PLS-DA of the samples, VIP scores, cross-validated accuracy and permutation p-value.

    PLS-DA finds the directions of protein space that best separate the conditions. With thousands of
    proteins and tens of samples it ALWAYS separates the groups of the samples it was fitted on, even
    when there is no biological difference, so the score plot alone proves nothing. The separation is
    therefore judged by:

    * cross-validated accuracy: stratified K-fold; the model (including z-scoring) is fitted on the
      training samples and classes of the held-out samples are predicted by the argmax of the
      predicted one-hot matrix;
    * permutation p-value: the whole cross-validation is repeated with shuffled condition labels;
      p = (1 + number of permuted accuracies >= observed) / (1 + permutations).

    Scores and VIP come from the model fitted on all samples and proteins passing the p-value cutoff.
    To keep the cross-validation and permutations fast, they use the `max_features` proteins with the
    highest variance across samples (a selection that ignores the labels, so it does not leak them).

    Args:
        pvalue (float, optional): keep proteins with self.pvalue <= pvalue; 1.0 keeps all. Defaults to 1.0.
            Selecting proteins by a p-value computed on the same labels makes the cross-validated accuracy
            optimistic; keep 1.0 to assess the separation.
        n_components (int, optional): PLS components. Defaults to 2.
        cv_folds (int, optional): folds of the stratified cross-validation. Defaults to 5.
        permutations (int, optional): label permutations; 0 skips the test. Defaults to 100.
        max_features (int, optional): proteins (highest variance) used in cross-validation and
            permutations; None uses all. Defaults to 5000.
        random_state (int, optional): seed for folds and permutations. Defaults to 0.

    Returns:
        dict: scores [Sample, Condition, comp1, comp2, ...], vip [Accession, gene_name, VIP] (descending),
        cv_accuracy, chance (accuracy of always predicting the largest class), permutation_p,
        permuted_accuracy (array), n_components.
    """
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.model_selection import StratifiedKFold

    qd = self.quant_data
    data = qd if pvalue >= 1 else qd[qd[self.pvalue] <= pvalue]
    logm, genes = _filled_log(self, data)
    if len(logm) < 2:
        raise ValueError('At least 2 proteins are needed for PLS-DA; relax the p-value cutoff.')
    cond = _condition_of(self)
    labels = np.array([cond[s] for s in logm.columns])
    levels = [c for c in self.Conditions if c in set(labels)]
    y = np.array([levels.index(c) for c in labels])
    X = logm.T.to_numpy(float)
    n_components = int(min(n_components, X.shape[0] - 1, X.shape[1]))

    pls = PLSRegression(n_components=n_components, scale=True).fit(X, np.eye(len(levels))[y])
    scores = pd.DataFrame(pls.transform(X), columns=[f'comp{i + 1}' for i in range(n_components)])
    scores.insert(0, 'Condition', labels)
    scores.insert(0, 'Sample', logm.columns)
    vip = pd.DataFrame({'Accession': logm.index, 'gene_name': genes.reindex(logm.index).to_numpy(),
                        'VIP': _vip(pls)}).sort_values('VIP', ascending=False, ignore_index=True)

    folds = int(min(cv_folds, np.bincount(y).min()))
    if folds < 2:
        raise ValueError('Each condition needs at least 2 samples for cross-validation.')
    Xcv = X
    if max_features is not None and X.shape[1] > max_features:
        Xcv = X[:, np.argsort(X.var(axis=0))[::-1][:max_features]]
    cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=random_state)
    accuracy = _plsda_cv_accuracy(Xcv, y, n_components, list(cv.split(Xcv, y)))
    rng = np.random.default_rng(random_state)
    permuted = np.empty(permutations)
    for b in range(permutations):
        yp = rng.permutation(y)
        permuted[b] = _plsda_cv_accuracy(Xcv, yp, n_components, list(cv.split(Xcv, yp)))
    p = (1 + np.sum(permuted >= accuracy)) / (1 + permutations) if permutations else np.nan
    return {'scores': scores, 'vip': vip, 'cv_accuracy': accuracy, 'chance': np.bincount(y).max() / len(y),
            'permutation_p': float(p), 'permuted_accuracy': permuted, 'n_components': n_components}


def plsda_plot(self, top=20, pvalue=1.0, n_components=2, cv_folds=5, permutations=100, max_features=5000,
               random_state=0, palette='Dark2', save=None, dpi=300, vector=True):
    """PLS-DA score plot and the proteins with the highest VIP.

    The score plot of a PLS-DA always shows separated groups, because the projection is chosen to
    separate them (with many proteins, even pure noise is separated). Read the title instead: the
    cross-validated accuracy should be well above chance and the permutation p small before the
    separation, or the VIP ranking, is interpreted. VIP > 1 is the usual threshold for proteins
    that contribute more than average to the projection.

    Args:
        top (int, optional): number of VIP proteins to show. Defaults to 20.
        pvalue (float, optional): p-value cutoff for proteins. Defaults to 1.0 (all).
        n_components (int, optional): PLS components. Defaults to 2.
        cv_folds (int, optional): cross-validation folds. Defaults to 5.
        permutations (int, optional): label permutations. Defaults to 100.
        max_features (int, optional): proteins used in cross-validation. Defaults to 5000.
        random_state (int, optional): seed. Defaults to 0.
        palette (str, optional): Palette for conditions. Defaults to 'Dark2'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.
    """
    res = plsda(self, pvalue=pvalue, n_components=max(2, n_components), cv_folds=cv_folds,
                permutations=permutations, max_features=max_features, random_state=random_state)
    scores, vip = res['scores'], res['vip'].head(top)
    colors = dict(zip(self.Conditions, sns.color_palette(palette, len(self.Conditions)).as_hex()))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(5, 0.28 * len(vip) + 1.5)),
                                   gridspec_kw={'width_ratios': [1.3, 1]}, dpi=dpi)
    sns.scatterplot(data=scores, x='comp1', y='comp2', hue='Condition', palette=colors, s=70, ax=ax1)
    p = res['permutation_p']
    ptxt = f', permutation p = {p:.3g}' if np.isfinite(p) else ''
    ax1.set(xlabel='comp 1', ylabel='comp 2',
            title=f"PLS-DA scores\nCV accuracy = {res['cv_accuracy']:.2f} (chance {res['chance']:.2f}){ptxt}")
    ax1.legend(frameon=False, bbox_to_anchor=(1, 1), loc='upper left')
    ax2.barh(range(len(vip))[::-1], vip['VIP'], color='#555')
    ax2.set_yticks(range(len(vip))[::-1], vip['gene_name'], fontsize=9)
    ax2.axvline(1, color='#d95f02', ls='--', lw=1)
    ax2.set(xlabel='VIP', title=f'Top {len(vip)} proteins by VIP')
    sns.despine(fig)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save + 'plsda.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return res


# ---------------------------------------------------------------------------- biomarkers


def _auc(X, y):
    """AUC of each column of X for y (1 = case): P(case > control) + 0.5 P(tie), via Mann-Whitney U."""
    from scipy.stats import rankdata

    ranks = rankdata(X, axis=0)
    n1, n0 = y.sum(), (1 - y).sum()
    u = ranks[y == 1].sum(axis=0) - n1 * (n1 + 1) / 2
    return u / (n1 * n0)


def roc_table(self, control=None, case=None):
    """Per-protein ROC AUC between two conditions, with 95% confidence interval.

    The AUC is the probability that a random `case` sample has a higher abundance than a random
    `control` sample (Mann-Whitney U / (n1 n0)). It is direction-aware: proteins lower in `case`
    get AUC = 1 - AUC and direction 'down', so AUC >= 0.5 and 0.5 means no discrimination.
    The confidence interval uses the Hanley & McNeil (1982) standard error; with few samples per
    group AUCs of 1 are common and the interval is wide.

    Args:
        control (str, optional): reference condition. Defaults to self.Conditions[0].
        case (str, optional): compared condition. Defaults to self.Conditions[1].

    Returns:
        pd.DataFrame: [Accession, gene_name, AUC, CI_low, CI_high, direction, p] sorted by AUC, with p the
        two-sided Mann-Whitney p-value.
    """
    from scipy.stats import mannwhitneyu

    X, y, genes, control, case = _two_groups(self, control, case)
    auc = _auc(X.to_numpy(float), y)
    a = np.maximum(auc, 1 - auc)
    n1, n0 = y.sum(), (1 - y).sum()
    q1, q2 = a / (2 - a), 2 * a ** 2 / (1 + a)
    se = np.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a ** 2) + (n0 - 1) * (q2 - a ** 2)) / (n1 * n0))
    p = mannwhitneyu(X.to_numpy(float)[y == 1], X.to_numpy(float)[y == 0], axis=0).pvalue
    out = pd.DataFrame({'Accession': X.columns, 'gene_name': genes.reindex(X.columns).to_numpy(),
                        'AUC': a, 'CI_low': np.clip(a - 1.96 * se, 0, 1), 'CI_high': np.clip(a + 1.96 * se, 0, 1),
                        'direction': np.where(auc >= 0.5, 'up', 'down'), 'p': p})
    return out.sort_values(['AUC', 'p'], ascending=[False, True], ignore_index=True)


def _t_abs(X, y):
    """absolute Welch t of each column of X between y == 1 and y == 0."""
    a, b = X[y == 1], X[y == 0]
    with np.errstate(invalid='ignore', divide='ignore'):
        t = (a.mean(0) - b.mean(0)) / np.sqrt(a.var(0, ddof=1) / len(a) + b.var(0, ddof=1) / len(b))
    return np.nan_to_num(np.abs(t))


def biomarker_panel(self, control=None, case=None, n_features=5, cv_folds=5, repeats=10, random_state=0):
    """Cross-validated logistic-regression panel of the most discriminant proteins.

    In each training fold the `n_features` proteins with the largest absolute Welch t are selected, z-scored
    and combined by an L2-penalised logistic regression; the held-out samples are then scored. Because
    the selection is repeated inside every fold, the held-out samples never influence which proteins
    are chosen: selecting the proteins on all samples first and cross-validating afterwards gives high
    AUCs even for pure noise. The cross-validation is repeated `repeats` times with different folds.

    Args:
        control (str, optional): reference condition. Defaults to self.Conditions[0].
        case (str, optional): compared condition. Defaults to self.Conditions[1].
        n_features (int, optional): proteins in the panel. Defaults to 5.
        cv_folds (int, optional): folds (capped at the size of the smallest group). Defaults to 5.
        repeats (int, optional): repetitions of the cross-validation. Defaults to 10.
        random_state (int, optional): seed. Defaults to 0.

    Returns:
        dict: fpr, tpr (ROC of the out-of-fold predictions of all repeats pooled), auc (pooled),
        auc_mean, auc_sd (across repeats), selection [Accession, gene_name, frequency] (share of the
        folds in which each protein entered the panel), features (panel chosen on all samples, for
        reporting) and control, case.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score, roc_curve
    from sklearn.model_selection import StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    Xdf, y, genes, control, case = _two_groups(self, control, case)
    X = Xdf.to_numpy(float)
    folds = int(min(cv_folds, np.bincount(y).min()))
    if folds < 2:
        raise ValueError('Each condition needs at least 2 samples for cross-validation.')
    n_features = int(min(n_features, X.shape[1]))
    counts = np.zeros(X.shape[1])
    pooled_y, pooled_p, aucs = [], [], []
    for r in range(repeats):
        prob = np.empty(len(y))
        for train, test in StratifiedKFold(folds, shuffle=True, random_state=random_state + r).split(X, y):
            keep = np.argsort(_t_abs(X[train], y[train]))[::-1][:n_features]   # selection on the training fold only
            counts[keep] += 1
            model = make_pipeline(StandardScaler(), LogisticRegression(C=1.0)).fit(X[train][:, keep], y[train])
            prob[test] = model.predict_proba(X[test][:, keep])[:, 1]
        aucs.append(roc_auc_score(y, prob))
        pooled_y.append(y)
        pooled_p.append(prob)
    pooled_y, pooled_p = np.concatenate(pooled_y), np.concatenate(pooled_p)
    fpr, tpr, _ = roc_curve(pooled_y, pooled_p)
    selection = pd.DataFrame({'Accession': Xdf.columns, 'gene_name': genes.reindex(Xdf.columns).to_numpy(),
                              'frequency': counts / (folds * repeats)})
    selection = selection[selection['frequency'] > 0].sort_values('frequency', ascending=False, ignore_index=True)
    final = np.argsort(_t_abs(X, y))[::-1][:n_features]
    return {'fpr': fpr, 'tpr': tpr, 'auc': float(roc_auc_score(pooled_y, pooled_p)),
            'auc_mean': float(np.mean(aucs)), 'auc_sd': float(np.std(aucs, ddof=1)) if repeats > 1 else 0.0,
            'selection': selection, 'features': list(genes.reindex(Xdf.columns[final])),
            'control': control, 'case': case}


def roc_plot(self, *Proteins, panel=True, control=None, case=None, top=5, n_features=5, cv_folds=5, repeats=10,
             random_state=0, palette='Dark2', save=None, dpi=300, vector=True):
    """ROC curves of single proteins and of a cross-validated multi-protein panel.

    Single-protein curves are computed on all samples (the protein was not chosen by the curve, but if
    it was picked as the best of thousands its AUC is optimistic). The panel curve is cross-validated
    with feature selection inside each fold (see biomarker_panel), so it estimates how the panel would
    perform on new samples.

    Args:
        *Proteins (str): gene names to plot; none plots the `top` proteins by AUC.
        panel (bool, optional): add the cross-validated panel. Defaults to True.
        control (str, optional): reference condition. Defaults to self.Conditions[0].
        case (str, optional): compared condition. Defaults to self.Conditions[1].
        top (int, optional): proteins plotted when none are given. Defaults to 5.
        n_features, cv_folds, repeats, random_state: see biomarker_panel.
        palette (str, optional): Palette for proteins. Defaults to 'Dark2'.
        save (str, optional): Path to save figure. Defaults to None.
        dpi (int, optional): figure resolution. Defaults to 300.
        vector (bool, optional): Save figure as vector (.svg). Defaults to True.

    Returns:
        (pd.DataFrame, dict or None): roc_table and biomarker_panel results.
    """
    from sklearn.metrics import roc_curve

    table = roc_table(self, control=control, case=case)
    X, y, genes, control, case = _two_groups(self, control, case)
    rows = table[table['gene_name'].isin(Proteins)] if Proteins else table.head(top)
    fig, ax = plt.subplots(figsize=(6.5, 6), dpi=dpi)
    for (_, r), color in zip(rows.iterrows(), sns.color_palette(palette, max(len(rows), 1))):
        score = X[r['Accession']].to_numpy(float) * (1 if r['direction'] == 'up' else -1)
        fpr, tpr, _ = roc_curve(y, score)
        ax.plot(fpr, tpr, color=color, lw=1.5, label=f"{r['gene_name']} ({r['direction']}, AUC {r['AUC']:.2f})")
    res = None
    if panel:
        res = biomarker_panel(self, control=control, case=case, n_features=n_features, cv_folds=cv_folds,
                              repeats=repeats, random_state=random_state)
        ax.plot(res['fpr'], res['tpr'], color='black', lw=2.5,
                label=f"CV panel of {n_features} (AUC {res['auc_mean']:.2f} ± {res['auc_sd']:.2f})")
    ax.plot([0, 1], [0, 1], color='#999', ls='--', lw=1)
    ax.set(xlabel='False positive rate (1 - specificity)', ylabel='True positive rate (sensitivity)',
           xlim=(0, 1), ylim=(0, 1.01), title=f'ROC: {case} vs {control}')
    ax.legend(frameon=False, loc='lower right', fontsize=9)
    sns.despine()
    if save is not None:
        plt.savefig(save + 'roc_plot.' + ('svg' if vector else 'png'), dpi=dpi, bbox_inches='tight')
    plt.show()
    return table, res
