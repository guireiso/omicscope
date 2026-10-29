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
    from .Meta import meta_analysis
    from .Meta import meta_plot
    from .Meta import forest_plot
    from .Meta import has_effect_sizes
    from .Integration import shared_universe
    from .Integration import overlap_test
    from .Integration import overlap_matrix
    from .Integration import concordance
    from .Integration import concordance_table
    from .Integration import concordance_plot
    from .Integration import study_diagnostics
    from .Integration import pvalue_histograms
    from .Sankey import sankey_plot
    from .Sankey import sankey_table

    def __init__(self, folder, palette='Dark2', pvalue_cutoff=0.05):
        self.original_path = os.getcwd()
        self.pvalue_cutoff = pvalue_cutoff
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
        from .omics_file import parse_omics
        all_files = sorted(glob.glob(folder + "/*.omics"))
        if len(all_files) == 0:
            raise ValueError('Nebula cannot import .omics file.' +
                             ' Check if the selected folder contains respective files.')
        parsed = []
        for path in all_files:
            with open(path, encoding='utf-8') as handle:
                parsed.append(parse_omics(handle.read()))
        self.pvalue = parsed[-1]['pvalue']
        self.groups = [p['group'] for p in parsed]
        self.labels = [p['label'] for p in parsed]
        self.original = [p['original'] for p in parsed]
        self.enrichment = [p['enrichment'] for p in parsed]
        self.metadata = [p['metadata'] for p in parsed]
        self.colors = sns.color_palette(palette, as_cmap=False, n_colors=len(self.groups)).as_hex()

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
        'meta_analysis',
        'meta_plot',
        'sankey_plot',
        'sankey_table',
        'forest_plot',
        'has_effect_sizes',
        'shared_universe',
        'overlap_test',
        'overlap_matrix',
        'concordance',
        'concordance_table',
        'concordance_plot',
        'study_diagnostics',
        'pvalue_histograms',
    ]
