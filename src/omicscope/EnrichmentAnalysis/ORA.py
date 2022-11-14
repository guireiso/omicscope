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
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

class Enrichment:
    """Over-Representation Analysis Results
    """
    def __init__(self, OmicScope, dbs=['KEGG_2021_Human'],
                 padjust_cutoff = 0.05,
                 organism='Human', save=''):
        """Over-Representation Analysis (also known as Enrichment Analysis)

        Args:
            OmicScope (OmicScope object): OmicScope Experiment
            dbs (list, optional): Databases to search against . Defaults to ['KEGG_2021_Human'].
            organism (str, optional): Organism to search in dbs. Defaults to 'Human'.
            save (str, optional): File path to save OmicScope Experiment.
            Defaults to ''.
        """
        self.Analysis = 'ORA'
        self.OmicScope = copy.copy(OmicScope)
        self.padjust_cutoff = padjust_cutoff,
        self.dbs = copy.copy(dbs)
        self.organism = copy.copy(organism)
        # Perform Over Representation Analysis
        try:
            self.results = copy.copy(self.enrichr())
        except NotADirectoryError:
            self.results = copy.copy(self.enrichr())
        # Saving
        if save != '':
            self.output(save=save)


    def enrichr(self):
        """Enrichr workflow according to gseapy
        """
        from gseapy import enrichr
        organism = self.organism
        pvalue_cutoff = self.OmicScope.PValue_cutoff
        foldchange_cutoff = self.OmicScope.FoldChange_cutoff
        gene_set = self.dbs
        omics = self.OmicScope.quant_data
        # Filtering data based on Fold Change and P-value
        omics = omics.loc[(omics['log2(fc)']<=-foldchange_cutoff) | (omics['log2(fc)']>=foldchange_cutoff)]
        genes = list(omics[omics['pvalue'] < pvalue_cutoff]['gene_name'].dropna())
        # Enrichment
        enr = enrichr(gene_list=genes, #proteins diferrentially regulated
                         gene_sets=gene_set, #dbs
                         organism=organism,
                         cutoff=1, outdir=None, no_plot=True)
        # Extracting just results file from analysis
        df = enr.results
        # Filtering results Adjusted P-value 
        df = df[df['Adjusted P-value'] < self.padjust_cutoff]
        df['N_Proteins'] = df['Overlap'].str.split('/').str[0].astype(int)
        df['-log10(pAdj)'] = -np.log10(df['Adjusted P-value'])
        foldchange = dict(zip(omics.gene_name.str.upper(), omics['log2(fc)']))
        df.Genes = df.Genes.str.split(';')
        df['regulation'] = df['Genes'].apply(lambda x: [foldchange[i] for i in x])
        df['down-regulated'] = df['regulation'].apply(lambda x: len([i for i in x if i < 0]))
        df['up-regulated'] = df['regulation'].apply(lambda x: len([i for i in x if i > 0]))
        return df