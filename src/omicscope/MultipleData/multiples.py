# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 15:21:23 2022

@author: gui_d
"""
import os
import pandas as pd
import glob
import __main__
import pandas as pd
from itertools import repeat
from copy import copy


class multiples:
    def __init__(self, folder):
        self.original_path = os.getcwd()
        self.read_omics(folder)
        self.group_data = []
        for o, g, l in zip(self.original, self.groups, self.labels):
            self.group_data.append(self.importing(o, g, l))

        print(f'''You imported your data successfully!
        Data description:
        1. N groups imported: {len(self.groups)} 
        2. Groups: {','.join(self.groups)} 
        3. N groups with enchment data: {sum(x is not None for x in self.enrichment)}
        ''')

    def read_omics(self, folder):
        path = folder  # use your path
        all_files = glob.glob(path + "/*.omics")
        groups = []
        labels = []
        original = []
        enrichment = []
        for i in all_files:
            positions = []
            with open(i, 'r') as archive:
                lines = archive.readlines()
                for n, line in enumerate(lines):
                    if line == '-------\n':
                        positions.append(n)
                    if line.startswith('Experimental'):
                        groups.append(line.split('\t')[-1].split('\n')[0])
                        labels.append(line.split('\t')[-1].split('\n')[0])
                    if line.startswith('Statistics'):
                        self.pvalue = line.split('\t')[-1].split('\n')[0]
                if len(positions) == 1:
                    original.append(pd.read_csv(i, header=positions[0]+1, sep='\t'))
                    enrichment.append(None)
                else:
                    original.append(pd.read_csv(i, header=positions[0]+1, sep='\t',
                                                nrows=int((positions[1]/2)-5)))
                    enrichment.append(pd.read_csv(i, header=int((positions[1]+1)/2)+4,
                                                  sep='\t'))
                archive.close()
            self.groups = groups
            self.labels = labels
            self.original = original
            self.enrichment = enrichment

    def importing(self, original, group, label):
        df = original
        df = df[df[self.pvalue] < 0.05][['gene_name', 'log2(fc)']].sort_values('log2(fc)', ignore_index=True)
        df['group'] = label
        df['color'] = df['log2(fc)'].round()
        return (df)

    from .MultipleVisualization import (barplot, protein_overlap,
                                        enrichment_overlap, overlap_pearson,
                                        Differentially_Regulated, network,
                                        overlap_stat, group_network)
    from .circos import circos_plot
