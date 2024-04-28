import glob
import os

import pandas as pd
import seaborn as sns


class nebula:
    from .MultipleVisualization import barplot
    from .MultipleVisualization import circos_plot
    from .MultipleVisualization import circular_term
    from .MultipleVisualization import diff_reg
    from .MultipleVisualization import dotplot_enrichment
    from .MultipleVisualization import enrichment_overlap
    from .MultipleVisualization import stat_heatmap
    from .MultipleVisualization import stat_network
    from .MultipleVisualization import protein_overlap
    from .MultipleVisualization import similarity_heatmap
    from .MultipleVisualization import similarity_network
    from .MultipleVisualization import whole_network

    def __init__(self, folder, palette='Dark2', pvalue_cutoff=0.05):
        self.original_path = os.getcwd()
        self.read_omics(folder, palette=palette)
        self.group_data = []
        for o, g, l in zip(self.original, self.groups, self.labels):
            self.group_data.append(self.importing(o, g, l, pvalue_cutoff=pvalue_cutoff))
        if len(self.groups) != len(set(self.groups)):
            # Create a dictionary to keep track of counts for each string
            string_count = {}
            group_without_duplicates = []
            # Loop through the original list and modify duplicates
            for item in self.groups:
                # If the string is already in the dictionary (a duplicate)
                if item in string_count:
                    # Increment the count
                    string_count[item] += 1
                    # Append the modified string with count to the second list
                    modified_item = f"{item}_{string_count[item]}"
                    group_without_duplicates.append(modified_item)
                else:
                    # If it's the first occurrence, add it as is to the second list
                    group_without_duplicates.append(item)
                    # Add the string to the dictionary with a count of 1
                    string_count[item] = 1
            self.groups = group_without_duplicates
            self.labels = group_without_duplicates
        self.enrichment = self.remap()
        print(f'''You imported your data successfully!
        Data description:
        1. N groups imported: {len(self.groups)}
        2. Groups: {','.join(self.groups)}
        3. N groups with enchment data: {sum(x is not None for x in self.enrichment)}
        ''')

    def read_omics(self, folder, palette):
        path = folder  # use your path
        all_files = glob.glob(path + "/*.omics")
        if len(all_files) == 0:
            raise ValueError('Nebula cannot import .omics file.' +
                             ' Check if the selected folder contains respective files.')
        groups = []
        labels = []
        original = []
        enrichment = []
        for i in all_files:
            positions = []
            with open(i, 'r') as archive:
                lines = archive.readlines()
                for n, line in enumerate(lines):
                    line = line.replace('\r\n', '\n')
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
                                                nrows=int(positions[1]-10)))
                    enrichment_original = pd.read_csv(i, header=int(positions[1]),
                                                      sep='\t')
                    enrichment_original['Genes'] = enrichment_original['Genes'].str.replace("'", '', regex=False)
                    enrichment_original['Genes'] = enrichment_original['Genes'].str.replace("[", '', regex=False)
                    enrichment_original['Genes'] = enrichment_original['Genes'].str.replace("]", '', regex=False)
                    enrichment_original['Genes'] = enrichment_original['Genes'].str.split(', ', regex=False)
                    enrichment.append(enrichment_original)
                archive.close()

            self.groups = groups
            self.colors = sns.color_palette(palette, as_cmap=False, n_colors=len(groups)).as_hex()
            self.labels = labels
            self.original = original
            self.enrichment = enrichment

    def importing(self, original, group, label, pvalue_cutoff):
        df = original
        df = df[df[self.pvalue] < pvalue_cutoff][['gene_name', 'log2(fc)']].sort_values('log2(fc)', ignore_index=True)
        df['group'] = label
        df['color'] = df['log2(fc)'].round()
        return (df)

    def remap(self):
        # Remap all genes from enrichment into respective gene annotation
        dic = pd.concat(self.original)[['gene_name']].drop_duplicates().reset_index(drop=True)
        dic.index = dic['gene_name'].str.upper()
        dic = dic.to_dict()
        dic = dic['gene_name']
        remap_enrichment = []
        original_enrichment = self.enrichment
        for i in original_enrichment:
            try:
                i['Genes'] = i['Genes'].apply(lambda x: [dic.get(i, i) for i in x])
                remap_enrichment.append(i)
            except TypeError:
                remap_enrichment.append(None)
        return remap_enrichment

    __all__ = [
        'circos_plot',
        'barplot',
        'circular_term',
        'similarity_heatmap',
        'similarity_network',
        'diff_reg',
        'dotplot_enrichment',
        'enrichment_overlap',
        'stat_heatmap',
        'stat_network',
        'whole_network',
        'protein_overlap',
    ]
