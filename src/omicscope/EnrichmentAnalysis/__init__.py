import copy
import warnings

warnings.filterwarnings("ignore")


class Enrichmentscope():
    from .EnrichmentVisualization import dotplot
    from .EnrichmentVisualization import enrichment_map
    from .EnrichmentVisualization import enrichment_network
    from .EnrichmentVisualization import gsea_heatmap
    from .EnrichmentVisualization import heatmap
    from .EnrichmentVisualization import number_deps

    def __init__(self, OmicScope, Analysis, dbs,
                 padjust_cutoff=0.05, organism='Human', background=None):
        """EnrichmentScope is the module designed to perform over-representation
            and Gene-Set Enrichment Analyses of proteins and genes.
            Args:
                OmicScope (Omicscope): Omicscope object
                Analysis (str): Over-representation Analysis (ORA) or Gene-Set Enrichment Analysis (GSEA).
                Defaults to 'ORA'.
                dbs (List[str]): List of enrichment databases to perform the enrichment analysis.
                Defaults to ['KEGG_2021_Human'].
                padjust_cutoff (float, optional): statistical cutoff. Defaults to 0.05.
                organism (str, optional): Organism from which entities belong. Defaults to 'human'.
                background (int, list, str, bool): Background genes. By default, all genes listed in the `gene_sets` input will be used 
                    as background. Alternatively, user can use all genes evaluated in study (Recommended, background = True). Still,
                    user can insert a specific number (integer) to use as background (Not recommended), such as number of reviewed proteins 
                    in the target organism on Uniprot.

                Raises:
                    ValueError: Users choose an invalid analysis.
                    Exception: Some functions just can be run in GSEA.
            Returns:
                Enrichmentscope: Return a EnrichmentScope obj. The results is stored to obj.results.

        """

        if Analysis == 'ORA':
            from .ORA import Enrichment
        elif Analysis == 'GSEA':
            from .GSEA import Enrichment
        else:
            raise ValueError('You must choose between ORA or GSEA workflow')
        self.OmicScope = copy.copy(OmicScope)
        self.Analysis = Analysis
        self.dbs = copy.copy(dbs)
        self.organism = copy.copy(organism)
        self.padjust_cutoff = copy.copy(padjust_cutoff)
        enrichment = Enrichment(OmicScope=OmicScope, dbs=dbs,
                                padjust_cutoff=padjust_cutoff, organism=organism,
                                background=background)
        self.results = enrichment.results

    def libraries(self):
        """Get Libraries from Enrichr
        """
        from gseapy import get_library_name
        libraries = get_library_name()
        return (libraries)

    def savefile(self, Path: str):
        from copy import copy
        import os
        save = Path
        data = copy(self)
        string = '-'.join(data.OmicScope.Conditions)
        with open(save + '/' + string + '.omics', 'w') as f:
            expression = data.OmicScope.quant_data[['gene_name', 'Accession',
                                                    data.OmicScope.pvalue, 'log2(fc)', 'TotalMean']].to_csv(sep='\t', index=False).replace(
                                                        '\r', "")
            if self.Analysis == 'ORA':
                enrichment = data.results[['Gene_set', 'Term', 'Overlap', 'Adjusted P-value', 'Genes']
                                          ].to_csv(sep='\t', index=False).replace('\r', "")
            elif self.Analysis == 'GSEA':
                enrichment = data.results[['Gene_set', 'Term', 'NES', 'Overlap',
                                           'Adjusted P-value', 'Genes']].to_csv(sep='\t', index=False).replace('\r', "")
            file = str("OmicScope" + "\n" +
                    "This file is the output performed by OmicScope pipeline and can be used as input" +
                    " for group comparisons having the controlling group used as used according to OmicScope." +
                    "Please, cite: Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: an Comprehensive Python " +
                    "library for Systems Biology Visualization" +
                    '\nControlGroup:' + '\t' + data.OmicScope.ctrl + '\n' +
                    'Experimental:' + '\t' + '\t'.join(data.OmicScope.experimental) + '\n' +
                    'Statistics:' + '\t' + data.OmicScope.pvalue + '\n' +
                    'Expression:\n' + '-------\n' +
                    expression + '\n' +
                    'Enrichment Analysis:\n' + '-------\n' +
                    enrichment)
            f.write(file.replace(os.linesep, '\n'))
    __all__ = ['dotplot',
               'heatmap',
               'number_deps',
               'enrichment_network',
               'enrichment_map',
               'gsea_heatmap']
