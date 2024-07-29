__version__ = '1.4.2'
from typing import List
from typing import Optional

from .EnrichmentAnalysis import Enrichmentscope
from .General.Omicscope import Omicscope
from .General.Snapshot import Omicscope_Snapshot
from .MultipleData.Nebula import nebula

if __name__ != "__main__":
    print("OmicScope v " + __version__ + " For help: https://omicscope.readthedocs.io/en/latest/ or https://omicscope.ib.unicamp.br " +
          "If you use  in published research, please cite:\n" +
          "'Reis-de-Oliveira, G., et al (2024). OmicScope unravels systems-level insights from quantitative proteomics data'.")


def OmicScope(Table: str,
              Method: str,
              ControlGroup: Optional[str] = None,
              ExperimentalDesign: str = 'static',
              pvalue: str = 'pAdjusted',
              PValue_cutoff: float = 0.05,
              normalization_method: Optional[str] = None,
              imputation_method: Optional[str] = None,
              FoldChange_cutoff: float = 0.0,
              logTransform: bool = True,
              ExcludeContaminants: bool = True,
              degrees_of_freedom: int = 2,
              independent_ttest=True,
              ** kwargs) -> Omicscope:
    """OmicScope - Differential Proteomics

        OmicScope was designed to be compatible with several Proteomics software,
        such as Progenesis Qi for Proteomics, PatternLab V, MaxQuant, 
        DIA-NN, ProteomeDiscoverer, and FragPipe.

        Additionally, users can also input data from other Omics sources (e.g.Transcriptomics),
         using `General` and `Snapshot` methods. In General, users can analyse data in a pre-specified format using excel workbooks.
         On the other hand, Snapshot enables users to import pre-analyzed data into OmicScope quickly.

        OmicScope are able to perform differential proteomics analysis, returning p-value, adjusted p-value
        (Benjamini-Hochberg approach), and  fold-changes.

    Args:
        Table (str): Quantitative data.
        Method (str): Method used to import data.
        ControlGroup (Optional[str], optional): Control group. Defaults to None.
        ExperimentalDesign (str, optional): Experimental design to perform statistical analysis.
         Options: 'static', 'longitudinal'. Defaults to 'static'.
        pvalue (str, optional): Statistical parameter to consider entities differentially regulated.
          Options: 'pvalue', 'pAdjusted', 'pTukey'. Defaults to 'pAdjusted'.
        PValue_cutoff (float, optional): Statistical cutoff. Defaults to 0.05.
        normalization_method (str, optional): Data normalization can be performed. Options:
          Options: 'average', 'median', 'quantile'. Defaults to None.
        imputation_method (str, optional): Impute values to data instead of NaN. 
          Options: "median", "mean", "knn". Defaults to None.
        FoldChange_cutoff (float, optional): Absolute fold-change cutoff. Defaults to 0.0.
        logTransform (bool, optional): Log-transform protein abundances. Defaults to False.
        ExcludeContaminants (bool, optional): Exclude the list of Contaminant proteins. Defaults to True.
        degrees_of_freedom (int, optional): Degrees of freedom to run longitudinal analysis. Defaults to 2.
        independent_ttest (bool, optional): while running a t-test, the user can specify if data sampling
            is independent (default) or paired (independent_ttest=False). Defaults to True.

    Returns:
        OmicScope: Return a OmicScope obj. The quantitation data is stored to obj.quant_data.
    """
    if Method == 'Snapshot':
        omics = Omicscope_Snapshot(
            Table,
            pvalue,
            PValue_cutoff,
            FoldChange_cutoff)
    else:
        omics = Omicscope(
            Table,
            Method,
            ControlGroup,
            ExperimentalDesign,
            pvalue,
            PValue_cutoff,
            normalization_method,
            imputation_method,
            FoldChange_cutoff,
            logTransform,
            ExcludeContaminants,
            degrees_of_freedom,
            independent_ttest,
            **kwargs)

    return omics


def EnrichmentScope(OmicScope: Omicscope,
                    Analysis: str = 'ORA',
                    dbs: List[str] = ['KEGG_2021_Human'],
                    padjust_cutoff: float = 0.05,
                    organism: str = 'human',
                    background=None) -> Enrichmentscope:
    """EnrichmentScope - Enrichment Analysis

        EnrichmentScope is the module designed to perform over-representation and Gene-Set Enrichment Analyses of proteins and genes.
        In EnrichmentScope, several figures enable user to see enriched terms with their respective proteins. 

    Args:
        OmicScope (Omicscope): Omicscope object
        Analysis (str): Over-representation Analysis (ORA) or Gene-Set Enrichment Analysis (GSEA). Defaults to 'ORA'.
        dbs (List[str]): List of enrichment databases to perform the enrichment analysis. Defaults to ['KEGG_2021_Human'].
        padjust_cutoff (float, optional): P-Adjusted cutoff . Defaults to 0.05.
        organism (str, optional): Organism to perform enrichment analysis. Defaults to 'human'.
        background (int, list, str, bool): Background genes. By default, all genes listed in the `gene_sets` input will be used 
            as background. Alternatively, user can use all genes evaluated in study (Recommended, background = True). Still,
            user can define a specific number (integer) to use as background (Not recommended), such as number of reviewed proteins 
            in the target organism on Uniprot.

    Returns:
        Enrichmentscope: Return a EnrichmentScope obj. The results is stored to obj.results.
    """
    enrichment = Enrichmentscope(
        OmicScope,
        Analysis,
        dbs,
        padjust_cutoff,
        organism,
        background)

    return enrichment


def Nebula(folder: str,
           palette: str = 'Dark2',
           pvalue_cutoff: float = 0.05) -> nebula:
    """Nebula - Multiple group comparison

        Nebula is the module to integrate all data generated by OmicScope and EnrichmentScope pipelines.

    Args:
        folder (str): path to folder that contains all .omics files
        palette (str): Palette to assign colors and discriminate groups
        pvalue_cutoff (float): P-value threshold to consider differentially regulated proteins

    Returns:
        Nebula: Return a Nebula obj.
    """
    multiple_groups = nebula(folder, palette=palette, pvalue_cutoff=pvalue_cutoff)

    return multiple_groups
