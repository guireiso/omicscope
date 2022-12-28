import copy
import pandas as pd


class EnrichmentScope():
    def __init__(self, OmicScope, Analysis, dbs = ['KEGG_2021_Human'], 
                 padjust_cutoff = 0.05, organism = 'Human'):
        """_summary_

        Args:
            OmicScope (_type_): _description_
            Analysis (_type_): _description_
            dbs (list, optional): _description_. Defaults to ['KEGG_2021_Human'].
            padjust_cutoff (float, optional): _description_. Defaults to 0.05.
            organism (str, optional): _description_. Defaults to 'Human'.
            save (str, optional): _description_. Defaults to ''.
        """
        if Analysis == 'ORA':
            from .ORA import Enrichment
            from .EnrichmentVisualization import (dotplot,
                                                heatmap,number_deps,
                                                enrichment_network,
                                                enrichment_map)
        elif Analysis == 'GSEA':
            from .GSEA import Enrichment
        else:
            raise ValueError('You must choose between ORA or GSEA workflow')
        self.OmicScope = copy.copy(OmicScope)
        self.Analysis = Analysis
        self.dbs = copy.copy(dbs)
        self.organism = copy.copy(organism)
        self.padjust_cutoff = copy.copy(padjust_cutoff)
        enrichment = Enrichment(OmicScope = OmicScope, dbs = dbs,
                     padjust_cutoff = padjust_cutoff, organism = organism)
        self.results = enrichment.results
    
    from .EnrichmentVisualization import (dotplot,
                    heatmap,number_deps,
                    enrichment_network,
                    enrichment_map, gsea_heatmap)

    def libraries(self):
        """Get Libraries from Enrichr
        """
        from gseapy import get_library_name
        libraries = get_library_name()
        return(libraries)

    def savefile(self, Path: str):
        from copy import copy
        save = Path
        data = copy(self)
        string = '-'.join(data.OmicScope.Conditions)
        with open(save + '/' + string + '.omics', 'w') as f:
            expression = data.OmicScope.quant_data[['gene_name', 'Accession', data.OmicScope.pvalue, 'log2(fc)']].to_csv(sep='\t', index=False)
            if self.Analysis == 'ORA':
                enrichment = data.results[['Gene_set', 'Term', 'Overlap', 'Adjusted P-value', 'Genes']].to_csv(sep='\t', index=False)
            elif self.Analysis == 'GSEA':
                enrichment = data.results[['Gene_set', 'Term', 'NES', 'Overlap', 'Adjusted P-value', 'Genes']].to_csv(sep='\t', index=False)
            
            f.write("OmicScope v1.0.0" + "\n" +
                    "This file is the output performed by OmicScope pipeline and can be used as input" +
                    " for group comparisons having the controling group used as used according to OmicScope." +
                    "Please, cite: Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: an Comprehensive Python library for Systems Biology Visualization" +
                    '\nControlGroup:' + '\t' + data.OmicScope.ctrl + '\n' +
                    'Experimental:' + '\t' + '\t'.join(data.OmicScope.experimental) + '\n' +
                    'Statistics:' +'\t'+ data.OmicScope.pvalue + '\n' +
                    'Expression:\n' + '-------\n' +
                    expression + '\n' +
                    'Enrichment Analysis:\n' + '-------\n' +
                    enrichment)
