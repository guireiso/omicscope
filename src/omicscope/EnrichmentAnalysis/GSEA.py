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
        from gseapy import gsea
        # pvalue_cutoff = copy.copy(self.OmicScope.PValue_cutoff)
        # foldchange_cutoff = copy.copy(self.OmicScope.FoldChange_cutoff)
        omics = copy.copy(self.OmicScope.quant_data)
        conditions = [copy.copy(self.OmicScope.ControlGroup)]
        conditions.extend(self.OmicScope.experimental)
        conditions = conditions[::-1]
        # GSEA DO NOT REQUIRE FILTERING STEPS TO ANALYSE DATA
        # omics = omics.loc[(omics['log2(fc)'] <= -foldchange_cutoff) | (omics['log2(fc)'] >= foldchange_cutoff)]
        # omics = omics[omics[self.OmicScope.pvalue] <= pvalue_cutoff]
        foldchange = dict(zip(omics.gene_name.str.upper(), omics['log2(fc)']))
        gsea_input = omics.set_index('gene_name')
        gsea_input = gsea_input.iloc[:, gsea_input.columns.str.contains('.', regex=False)]
        gsea_input = np.log2(gsea_input)
        gsea_input = gsea_input.replace(-np.inf, 0)
        gsea_input = gsea_input.replace(np.nan, 0)
        sorting_columns = [column for cond in conditions for column in gsea_input.columns if column.endswith(cond)]
        gsea_input = gsea_input.reindex(columns=sorting_columns)
        groups = gsea_input.columns.str.split('.', regex=False).str[1]
        gsea_input = gsea_input.reset_index()
        gsea_input['gene_name'] = gsea_input['gene_name'].str.upper()
        gsea_input = gsea_input[gsea_input['gene_name'].notna()]
        gsea_result = gsea(gsea_input, gene_sets=db,
                           cls=groups, no_plot=True, outdir=None)
        gsea_result = gsea_result.res2d
        gsea_result = gsea_result.reset_index()
        gsea_result['DB'] = db
        gsea_result = gsea_result.rename({'DB': 'Gene_set',
                                          'Lead_genes': 'Genes',
                                          'FDR q-val': 'Adjusted P-value',
                                          'NOM p-val': 'P-value',
                                          'Tag %': 'Overlap'},
                                         axis='columns')

        gsea_result = gsea_result[gsea_result['Adjusted P-value'] <= self.padjust_cutoff]
        gsea_result['Genes'] = gsea_result['Genes'].str.split(';')
        gsea_result['regulation'] = gsea_result['Genes'].apply(lambda x: [foldchange[i] for i in x if x != ['']])
        gsea_result['down-regulated'] = gsea_result['regulation'].apply(lambda x: len([i for i in x if i < 0]))
        gsea_result['up-regulated'] = gsea_result['regulation'].apply(lambda x: len([i for i in x if i > 0]))
        gsea_result['N_Proteins'] = gsea_result['down-regulated'] + gsea_result['up-regulated']
        # try:
        gsea_result['-log10(pAdj)'] = gsea_result['Adjusted P-value'].apply(lambda x: -np.log10(x) if type(x) is float else 0)
        # except TypeError:
        #     gsea_result['-log10(pAdj)'] = 1
        gsea_result = gsea_result[['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value',
                                   'NES', 'ES', 'Genes', 'N_Proteins', '-log10(pAdj)', 'regulation',
                                   'down-regulated', 'up-regulated']]
        return gsea_result
