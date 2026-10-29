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
    from .EnrichmentExtras import gsea_plot
    from .ClusterEnrichment import collapse_terms
    from .EnrichmentExtras import term_foldchange
    from .EnrichmentExtras import term_foldchange_table

    def __init__(self, OmicScope, Analysis, dbs,
                 padjust_cutoff=0.05, organism='Human', background=None, **gsea_options):
        """EnrichmentScope is the module designed to perform over-representation
            and Gene-Set Enrichment Analyses of proteins and genes.
            Args:
                OmicScope (Omicscope): Omicscope object
                Analysis (str): Over-representation Analysis (ORA) or Gene-Set Enrichment Analysis (GSEA).
                Defaults to 'ORA'.
                dbs (List[str]): List of enrichment databases to perform the enrichment analysis:
                Enrichr library names, gene-set dicts (term -> genes; Gene_set 'custom') or paths to .gmt
                files (Gene_set = file name without extension). Defaults to ['KEGG_2021_Human'].
                padjust_cutoff (float, optional): statistical cutoff. Defaults to 0.05.
                organism (str, optional): Organism from which entities belong. Defaults to 'human'.
                background (int, list, str, bool): Background genes. By default, all genes listed in the `gene_sets` input will be used 
                    as background. Alternatively, user can use all genes evaluated in study (Recommended, background = True). Still,
                    user can insert a specific number (integer) to use as background (Not recommended), such as number of reviewed proteins 
                    in the target organism on Uniprot.
                **gsea_options: GSEA only. rank_by ('signed_p' (default): sign(log2FC) * -log10(p); 'log2fc';
                    'pi_score'), min_size (default 10: proteomics covers a fraction of each gene set),
                    max_size (500), permutation_num (1000), seed (123) and threads (4).

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
        if isinstance(dbs, dict) or (isinstance(dbs, str) and dbs.lower().endswith('.gmt')):
            dbs = [dbs]
        self.OmicScope = copy.copy(OmicScope)
        self.Analysis = Analysis
        self.dbs = copy.copy(dbs)
        self.organism = copy.copy(organism)
        self.padjust_cutoff = copy.copy(padjust_cutoff)
        if gsea_options and Analysis != 'GSEA':
            raise ValueError(f'Options {sorted(gsea_options)} only apply to GSEA')
        enrichment = Enrichment(OmicScope=OmicScope, dbs=dbs,
                                padjust_cutoff=padjust_cutoff, organism=organism,
                                background=background, **gsea_options)
        self.results = enrichment.results
        if Analysis == 'GSEA':
            self.ranking = enrichment.ranking   # the ranked list GSEA used, for running-score plots
            self.curves = enrichment.curves
            self.rank_by = enrichment.rank_by

    def libraries(self):
        """Get Libraries from Enrichr
        """
        from gseapy import get_library_name
        libraries = get_library_name()
        return (libraries)

    def savefile(self, Path: str):
        """Save the .omics file with the enrichment results (for Nebula) and a pickle of the object."""
        import pickle

        from ..MultipleData.omics_file import omics_text
        string = '-'.join(self.OmicScope.Conditions)
        with open(Path + '/' + string + '.omics', 'w', encoding='utf-8', newline='\n') as f:
            f.write(omics_text(self.OmicScope, enrichment=self))
        with open(Path + '/' + string + '.results', 'wb') as bi:
            pickle.dump(self, bi)

    __all__ = ['dotplot',
               'heatmap',
               'number_deps',
               'enrichment_network',
               'enrichment_map',
               'gsea_heatmap',
               'gsea_plot',
               'collapse_terms',
               'term_foldchange',
               'term_foldchange_table']
