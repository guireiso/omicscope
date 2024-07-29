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


class Enrichment:
    """Over-Representation Analysis Results
    """

    def __init__(self, OmicScope, dbs,
                 padjust_cutoff=0.05,
                 organism='Human', background=None, save=''):
        """Over-Representation Analysis (also known as Enrichment Analysis)

        Args:
            OmicScope (OmicScope object): OmicScope Experiment
            dbs (list, optional): Databases to search against . Defaults to ['KEGG_2021_Human'].
            organism (str, optional): Organism to search in dbs. Defaults to 'Human'.
            background (int, list, str, bool): Background genes. By default, all genes listed in the `gene_sets` input will be used 
                as background. Alternatively, user can use all genes evaluated in study (Recommended, background = True). Still,
                user can insert a specific number (integer) to use as background (Not recommended), such as number of reviewed proteins 
                in the target organism on Uniprot.
            save (str, optional): File path to save OmicScope Experiment.
            Defaults to ''.
        """
        self.Analysis = 'ORA'
        self.OmicScope = copy.copy(OmicScope)
        self.padjust_cutoff = padjust_cutoff,
        self.dbs = copy.copy(dbs)
        self.organism = copy.copy(organism)
        self.background = background
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
        background = self.background
        if background is True:
            background = omics.gene_name.dropna().drop_duplicates()
        # Filtering data based on Fold Change and P-value
        omics = omics.loc[(omics['log2(fc)'] <= -foldchange_cutoff) | (omics['log2(fc)'] >= foldchange_cutoff)]
        genes = list(omics[omics[self.OmicScope.pvalue] <= pvalue_cutoff]['gene_name'].dropna())
        # Enrichment
        enr = enrichr(gene_list=genes,  # proteins diferrentially regulated
                      gene_sets=gene_set,
                      outdir=None,  # dbs
                      organism=organism,
                      background=background,
                      cutoff=1, no_plot=True,
                      verbose=False)
        # Extracting just results file from analysis
        df = enr.results.copy()
        # Filtering results Adjusted P-value
        df = df[df['Adjusted P-value'] <= self.padjust_cutoff]
        # GSEApy package presents a threshold to consider pvalue = 0
        try:
            minimum_value = df['Adjusted P-value'].drop_duplicates().sort_values()
            minimum_value = minimum_value[minimum_value != 0].reset_index(drop=True)
            minimum_value = minimum_value[0]
            df['Adjusted P-value'] = df['Adjusted P-value'].replace(0, minimum_value)
        except KeyError:
            pass
        # Adjusting table
        df['-log10(pAdj)'] = -np.log10(df['Adjusted P-value'])
        foldchange = dict(zip(omics.gene_name.str.upper(), omics['log2(fc)']))
        df.Genes = df.Genes.str.split(';')
        df['N_Proteins'] = df['Genes'].apply(lambda x: len(x))
        df['regulation'] = df['Genes'].apply(lambda x: [foldchange[i.upper()] for i in x])
        df['down-regulated'] = df['regulation'].apply(lambda x: len([i for i in x if i < 0]))
        df['up-regulated'] = df['regulation'].apply(lambda x: len([i for i in x if i > 0]))
        df = df.sort_values(['Adjusted P-value', 'Combined Score'], 
                            ascending=[True, False])
        df = df.reset_index()
        return df
