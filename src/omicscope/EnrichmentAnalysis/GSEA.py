""" Module for Functional Enrichment Analysis.

This OmicScope module is totally dedicated to Functional Enrichment
Analysis of Omics Experiment.

In this workflow we used GSEApy algorithm to
get two main analysis: Over-Representation Analysis and Gene-Set Enrichment
Analysis. For both cases, GSEApy runs the analysis against Enrichr databases.

Additionally to prepare data to perform enrichment analysis,
OmicScope perform several analysis associated to enrichment analysis and
also can produces ready-to-publish figures.

@author: Reis-de-Oliveira G <guioliveirareis@gmail.com>
"""
import copy
import warnings

import numpy as np
import pandas as pd

RANK_METRICS = ('signed_p', 'log2fc', 'pi_score')


def rank_entities(quant_data, rank_by='signed_p', pvalue='pvalue'):
    """Ranked list used by pre-ranked GSEA.

    signed_p: sign(log2FC) * -log10(p). With a moderated t-test every protein shares the same degrees of
        freedom, so this is the moderated t statistic on a monotone scale: precise and large changes go to
        the extremes, noisy ones stay near the middle. Recommended.
    log2fc: fold change only (the former OmicScope behaviour); ignores variability.
    pi_score: log2FC * -log10(p) (Xiao et al. 2014).

    Proteins without gene name or statistics are dropped and duplicated gene names keep the entry with the
    largest absolute score, as GSEA needs one value per gene.
    """
    if rank_by not in RANK_METRICS:
        raise ValueError(f'rank_by must be one of {RANK_METRICS}')
    if pvalue not in quant_data.columns:
        pvalue = 'pAdjusted' if 'pAdjusted' in quant_data.columns else None
    df = quant_data[['gene_name', 'log2(fc)'] + ([pvalue] if pvalue else [])].copy()
    df['gene_name'] = df['gene_name'].astype('string').str.strip()
    df = df[df['gene_name'].notna() & (df['gene_name'] != '') & (df['gene_name'].str.lower() != 'nan')]
    fc = pd.to_numeric(df['log2(fc)'], errors='coerce')
    if rank_by == 'log2fc' or pvalue is None:
        score = fc
    else:
        p = pd.to_numeric(df[pvalue], errors='coerce').clip(lower=1e-300, upper=1)
        score = np.sign(fc) * -np.log10(p) if rank_by == 'signed_p' else fc * -np.log10(p)
        # ties at p = 1 are broken by the fold change, so the order inside the middle of the list is still meaningful
        score = score + fc * 1e-6
    df['score'] = score
    df = df[np.isfinite(df['score'])]
    df = df.loc[df['score'].abs().groupby(df['gene_name']).idxmax()]
    return df.set_index('gene_name')['score'].sort_values(ascending=False)


class Enrichment:
    def __init__(self, OmicScope, dbs,
                 organism='Human', padjust_cutoff=0.05, background=None,
                 rank_by='signed_p', min_size=10, max_size=500, permutation_num=1000, seed=123, threads=4):

        self.Analysis = 'GSEA'
        self.dbs = copy.copy(dbs)
        self.OmicScope = copy.copy(OmicScope)
        self.organism = str(organism).lower()   # gseapy only accepts lower-case organism names
        self.padjust_cutoff = copy.copy(padjust_cutoff)
        self.params = dict(min_size=min_size, max_size=max_size, permutation_num=permutation_num,
                           seed=seed, threads=threads)
        conditions = getattr(OmicScope, 'Conditions', None)
        if conditions is not None and len(conditions) > 2:
            warnings.warn('GSEA ranks proteins by a single comparison, but this dataset has more than two '
                          'conditions: log2(fc) of each protein comes from its own post-hoc pair. Prefer a '
                          'two-group OmicScope object (e.g. drop the other groups) for GSEA.')
        self.ranking = rank_entities(OmicScope.quant_data, rank_by=rank_by)
        self.rank_by = rank_by
        self.curves = {}
        dfs = [self.gsea(db=i) for i in dbs]
        results = pd.concat(dfs)
        self.results = results.sort_values('Adjusted P-value', ignore_index=True)

    def gsea(self, db):
        """GSEA workflow according to gseapy
        """
        from gseapy import prerank

        from .ClusterEnrichment import gene_set_library
        db, name = gene_set_library(db)   # .gmt files are read here; dicts are named 'custom'
        foldchange = dict(zip(self.OmicScope.quant_data['gene_name'].astype(str).str.upper(),
                              pd.to_numeric(self.OmicScope.quant_data['log2(fc)'], errors='coerce')))
        gsea_result = prerank(self.ranking, gene_sets=db, organism=self.organism,
                              no_plot=True, outdir=None, **self.params)
        # running enrichment score and hit positions of each term, for the classic GSEA plot
        self.curves.update({(name, k): {'hits': np.asarray(v['hits'], int),
                                        'RES': np.asarray(v['RES'], np.float32)}
                            for k, v in gsea_result.results.items()
                            if float(v['fdr']) <= self.padjust_cutoff})
        gsea_result = pd.DataFrame.from_dict(gsea_result.results).T
        if gsea_result.empty:
            return pd.DataFrame(columns=['Term', 'ES', 'NES', 'P-value', 'Adjusted P-value', 'Genes', 'Gene_set'])
        gsea_result = gsea_result.reset_index()
        gsea_result = gsea_result.rename(columns={
                                'index': 'Term',
                                'es': 'ES',
                                'nes': 'NES',
                                'pval': 'P-value',
                                'fdr': 'Adjusted P-value',
                                'fwerp': 'FWER p-val',
                                'tag %': 'Overlap',
                                'gene %': 'Gene %',
                                'lead_genes': 'Genes',
                                })
        for col in ['ES', 'NES', 'P-value', 'Adjusted P-value', 'FWER p-val']:
            gsea_result[col] = pd.to_numeric(gsea_result[col], errors='coerce')
        gsea_result = gsea_result.sort_values('Adjusted P-value').reset_index(drop=True)
        gsea_result['Gene_set'] = name
        gsea_result = gsea_result[gsea_result['Adjusted P-value'] <= self.padjust_cutoff].copy()
        gsea_result['Genes'] = gsea_result['Genes'].astype(str).str.split(';').apply(lambda x: [g for g in x if g])
        # one value per gene (NaN if the symbol is not found), so Genes and regulation stay aligned
        gsea_result['regulation'] = gsea_result['Genes'].apply(lambda x: [foldchange.get(g.upper(), np.nan) for g in x])
        gsea_result['down-regulated'] = gsea_result['regulation'].apply(lambda x: sum(1 for i in x if i < 0))
        gsea_result['up-regulated'] = gsea_result['regulation'].apply(lambda x: sum(1 for i in x if i > 0))
        gsea_result['N_Proteins'] = gsea_result['down-regulated'] + gsea_result['up-regulated']
        # permutation FDR can be exactly 0: report it as below the resolution of the permutations
        floor = 1 / (self.params['permutation_num'] + 1)
        gsea_result['-log10(pAdj)'] = -np.log10(gsea_result['Adjusted P-value'].clip(lower=floor))
        return gsea_result
