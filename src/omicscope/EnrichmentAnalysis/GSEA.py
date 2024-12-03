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

import numpy as np
import pandas as pd


class Enrichment:
    def __init__(self, OmicScope, dbs,
                 organism='Human', padjust_cutoff=0.05, background=None):

        self.Analysis = 'GSEA'
        self.dbs = copy.copy(dbs)
        self.OmicScope = copy.copy(OmicScope)
        self.organism = copy.copy(organism)
        self.padjust_cutoff = copy.copy(padjust_cutoff)
        dfs = []
        for i in dbs:
            try:
                df = self.gsea(db=i)
            except NotADirectoryError:
                df = self.gsea(db=i)
            dfs.append(df)
        results = pd.concat(dfs)
        self.results = results.sort_values('Adjusted P-value', ignore_index=True)

    def gsea(self, db):
        """GSEA workflow according to gseapy
        """
        from gseapy import prerank
        omics = copy.copy(self.OmicScope.quant_data)
        omics = omics[['gene_name', "log2(fc)"]]
        omics = omics.sort_values('log2(fc)',ascending=False).reset_index(drop=True)
        foldchange = dict(zip(omics.gene_name.str.upper(), omics['log2(fc)']))
        gsea_result = prerank(omics, gene_sets=db,
                              no_plot=True, outdir=None,)
        gsea_result = pd.DataFrame.from_dict(gsea_result.results).T
        gsea_result = gsea_result.sort_values('fdr')
        gsea_result = gsea_result.reset_index()
        gsea_result = gsea_result.rename(columns={
                                'index':'Term',
                                'es':'ES',
                                'nes':'NES',
                                'pval':'P-value',
                                'fdr':'Adjusted P-value',
                                'fwerp':'FWER p-val',
                                'tag %': 'Overlap',
                                'gene %':'Gene %',
                                'lead_genes':'Genes',                        
                                })
        gsea_result = gsea_result.reset_index()
        gsea_result['Gene_set'] = db
        gsea_result = gsea_result[gsea_result['Adjusted P-value'] <= self.padjust_cutoff]
        gsea_result['Genes'] = gsea_result['Genes'].str.split(';')
        gsea_result['regulation'] = gsea_result['Genes'].apply(lambda x: [foldchange[i] for i in x if x != ['']])
        gsea_result['down-regulated'] = gsea_result['regulation'].apply(lambda x: len([i for i in x if i < 0]))
        gsea_result['up-regulated'] = gsea_result['regulation'].apply(lambda x: len([i for i in x if i > 0]))
        gsea_result['N_Proteins'] = gsea_result['down-regulated'] + gsea_result['up-regulated']
        gsea_result['-log10(pAdj)'] = gsea_result['Adjusted P-value'].apply(lambda x: -np.log10(x) if type(x) is float else 0)
        return gsea_result
