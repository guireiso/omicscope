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

ENRICHR_COLUMNS = ['Gene_set', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Odds Ratio', 'Combined Score',
                   'Genes']


class Enrichment:
    """Over-Representation Analysis Results
    """

    def __init__(self, OmicScope, dbs,
                 padjust_cutoff=0.05,
                 organism='Human', background=None, save=''):
        """Over-Representation Analysis (also known as Enrichment Analysis)

        Args:
            OmicScope (OmicScope object): OmicScope Experiment
            dbs (list, optional): Databases to search against: Enrichr libraries, gene-set dicts
                (term -> genes) or paths to .gmt files. Defaults to ['KEGG_2021_Human'].
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
        self.padjust_cutoff = padjust_cutoff
        self.dbs = copy.copy(dbs)
        self.organism = str(organism).lower()   # gseapy only accepts lower-case organism names
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
        pvalue_cutoff = self.OmicScope.PValue_cutoff
        foldchange_cutoff = self.OmicScope.FoldChange_cutoff
        omics = self.OmicScope.quant_data
        background = self.background
        if background is True:
            background = omics.gene_name.dropna().drop_duplicates()
        # Filtering data based on Fold Change and P-value
        omics = omics.loc[(omics['log2(fc)'] <= -foldchange_cutoff) | (omics['log2(fc)'] >= foldchange_cutoff)]
        genes = list(omics[omics[self.OmicScope.pvalue] <= pvalue_cutoff]['gene_name'].dropna())
        foldchange = dict(zip(omics.gene_name.str.upper(), omics['log2(fc)']))
        return run_enrichr(genes, self.dbs, organism=self.organism, background=background,
                           padjust_cutoff=self.padjust_cutoff, foldchange=foldchange)


def run_enrichr(genes, dbs, organism='human', background=None, padjust_cutoff=0.05, foldchange=None):
    """Over-representation of a gene list with gseapy enrichr.

    Args:
        genes (list): gene names to test.
        dbs (list): Enrichr libraries, gene-set dicts (term -> genes) or paths to .gmt files. Dicts and
            .gmt files run offline; their Gene_set is 'custom' and the file stem, respectively.
        organism (str, optional): organism of the Enrichr libraries. Defaults to 'human'.
        background (list, int, optional): background genes (see Enrichment). Defaults to None.
        padjust_cutoff (float, optional): adjusted p-value cutoff. Defaults to 0.05.
        foldchange (dict, optional): upper-case gene name -> log2(fc), for the regulation columns.

    Returns:
        pd.DataFrame: enrichr results with Genes as lists, -log10(pAdj), N_Proteins, regulation,
        down-regulated and up-regulated.
    """
    from gseapy import enrichr

    from .ClusterEnrichment import gene_set_library
    dbs = [dbs] if isinstance(dbs, (str, dict)) else list(dbs)
    libraries = [gene_set_library(db) for db in dbs]
    foldchange = {} if foldchange is None else foldchange
    # Enrichment
    enr = enrichr(gene_list=list(genes),  # proteins diferrentially regulated
                  gene_sets=[lib for lib, _ in libraries],
                  outdir=None,  # dbs
                  organism=str(organism).lower(),
                  background=background,
                  cutoff=1, no_plot=True,
                  verbose=False)
    # Extracting just results file from analysis
    df = enr.results.copy() if isinstance(enr.results, pd.DataFrame) else pd.DataFrame(columns=ENRICHR_COLUMNS)
    # gseapy names gene-set dicts after their position in the list
    df['Gene_set'] = df['Gene_set'].replace({f'gs_ind_{i}': name for i, (_, name) in enumerate(libraries)})
    # Filtering results Adjusted P-value
    df = df[df['Adjusted P-value'] <= padjust_cutoff].copy()
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
    df['Genes'] = df['Genes'].str.split(';')
    df['N_Proteins'] = df['Genes'].apply(lambda x: len(x))
    df['regulation'] = df['Genes'].apply(lambda x: [foldchange.get(i.upper(), np.nan) for i in x])
    df['down-regulated'] = df['regulation'].apply(lambda x: len([i for i in x if i < 0]))
    df['up-regulated'] = df['regulation'].apply(lambda x: len([i for i in x if i > 0]))
    df = df.sort_values(['Adjusted P-value', 'Combined Score'],
                        ascending=[True, False])
    df = df.reset_index()
    return df
