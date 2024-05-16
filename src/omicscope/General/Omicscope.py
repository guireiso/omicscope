"""  OmicScope Module

OmicScope allows user to import data and perform statistical analysis
for differential gene/protein expression. This is the key module of
all OmicScope workflow, being the input for all other modules.

"""
import warnings

from .Input import Input
from .Stats.imputation import value_imputation
from .Stats.normalization import normalization

warnings.filterwarnings("ignore")


class Omicscope(Input):
    from .GeneralVisualization import DynamicRange
    from .GeneralVisualization import MAplot
    from .GeneralVisualization import PPInteractions
    from .GeneralVisualization import bar_ident
    from .GeneralVisualization import bar_protein
    from .GeneralVisualization import boxplot_protein
    from .GeneralVisualization import correlation
    from .GeneralVisualization import heatmap
    from .GeneralVisualization import k_trend
    from .GeneralVisualization import pca
    from .GeneralVisualization import volcano
    from .GeneralVisualization import normalization_boxplot

    def __init__(self, Table, Method, ControlGroup, ExperimentalDesign='static',
                 pvalue='pAdjusted', PValue_cutoff=0.05, 
                 normalization_method = None, imputation_method=None,
                 FoldChange_cutoff=0, logTransform=True, ExcludeContaminants=True,
                 degrees_of_freedom=2, independent_ttest=True, **kwargs):
        """  OmicScope was specially designed taking into account the
        proteomic workflow, in which proteins are identified, quantified
        and normalized with several softwares, such as Progenesis Qi for
        Proteomics, PatternLab V and MaxQuant.

        Despite proteomic optimization, users can also input data from other
        Omics sources (e.g. Metabolomics and Transcriptomics), using
        Method = 'General'. To run this approach user need to follow the
        instructions from Input.General directory in a Excel file.

        Additionally, OmicScope also perform differential expression analysis,
        returning p-value (from t-test, Anova, and post-hoc correction),
        q-value and respective conditions Fold-Changes. If user desire to
        import the statistical analysis from other sources, the row data
        (rdata) must have columns "pvalue".

        Args:
            Table (str): Path to quantitative data.
            Method (str): Algorithm/software used to perform quantitative.
            ControlGroup (Optional[str], optional): Control group. Defaults to None.
            ExperimentalDesign (str, optional): Study experimental design in which
            OmicScope must take account for statistical analysis. Defaults to 'static'.
            pvalue (str, optional): Statistical parameter to take into account
            to consider entities differentially regulated. Defaults to 'pAdjusted'.
            pdata (Optional[str], optional): Path to phenotype data of each sample. Defaults to None.
            PValue_cutoff (float, optional): Statistical cutoff. Defaults to 0.05.
            normalization_method (str, optional): Data normalization can be performed. Options:
              Options: 'average', 'median', 'quantile'. Defaults to None.
            imputation_method (str, optional): Impute values to data instead of NaN. 
              Options: "median", "mean", "knn". Defaults to None.
            FoldChange_cutoff (float, optional): Difference cutoff. Defaults to 0.0.
            logTransform (bool, optional): Log-transform protein abundances. Defaults to True.
            ExcludeContaminants (bool, optional): Contaminant proteins is excluded. Defaults to True.
            degrees_of_freedom (int, optional): Degrees of freedom used to run longitudinal analysis. Defaults to 2.
            independent_ttest (bool, optional): If running a t-test, the user can specify if samples
                are independent (default) or paired (independent_ttest=False). Defaults to True.
        """
        import pandas as pd
        self.Params = {'Params': {}}
        super().__init__(Table, Method=Method, **kwargs)
        self.PValue_cutoff = PValue_cutoff
        self.FoldChange_cutoff = FoldChange_cutoff
        self.logTransform = logTransform
        self.normalization_method = normalization_method
        self.imputation_method = imputation_method
        self.ExcludeContaminants = ExcludeContaminants
        self.pvalue = pvalue
        self.ControlGroup = ControlGroup
        self.ind_variables = independent_ttest
        pvalues = ['pvalue', 'pAdjusted', 'pTukey']
        if pvalue not in pvalues:
            raise ValueError(
                "Invalid pvalue specification. Expected one of: %s" % pvalues)
        if 'pdata' in kwargs:
            # If pdata was assigned by user, OmicScope read
            # excel or csv frames.
            try:
                self.pdata = pd.read_excel(kwargs['pdata'])
            except ValueError:
                self.pdata = pd.read_csv(kwargs['pdata'])
                if self.pdata.shape[1] == 1:
                    raise KeyError('OmicScope does not import pdata correctly. Please, confirm if data is xlsx, xls or csv.')

        self.define_conditions()

        self.ctrl = self.ControlGroup
        #  Has user already performed statistical analyses?
        statistics = ['Anova (p)', 'pvalue', 'p-value', 'q-value', 'q Value', 'qvalue',
                      'padj', 'p-adjusted', 'p-adj', 'pAdjusted']
        if True in self.rdata.columns.isin(statistics):
            from .Stats.Performed_Stat import imported_stat
            self.quant_data = imported_stat(self, statistics)
            self.expression = self.assay.copy()
            self.pdata['Sample'] = self.pdata['Sample'] + \
                '.'+self.pdata['Condition']
            print('User already performed statistical analysis')

        elif ExperimentalDesign == 'static':  # OmicScope perform statistics
            from .Stats.Statistic_Module import perform_static_stat
            # Construct pivot-table considering technical and biological replicates
            expression = self.expression_function()
            expression = normalization(self, expression)
            expression = value_imputation(self, expression)
            if len(expression)==0:
                raise ValueError("""
                OmicScope filtered out all proteins in our dataset.
                Please verify if quantitative values were correctly inserted and/or if they differ from zero.""")
            self.expression = expression.copy()
                
            # perform stat
            self.quant_data = perform_static_stat(self)
            print('OmicScope performed statistical analysis (Static workflow)')
        
        elif ExperimentalDesign == 'longitudinal':
            from .Stats.Statistic_Module import perform_longitudinal_stat
            expression = self.expression_function()
            expression = normalization(self, expression)
            expression = value_imputation(self, expression)
            self.expression = expression.copy()

            # perform stat
            self.degrees_of_freedom = degrees_of_freedom
            self.quant_data = perform_longitudinal_stat(self)
            print('OmicScope performed statistical analysis (Longitudinal workflow)')
        else:
            raise ValueError("Invalid statistical specification. " +
                             "User did not add statistical column into rdata, " +
                             "neither specify 'static' or 'longitudinal' experimental design. "
                             )
        self.deps = self.deps()
        if len(self.deps) == 0:
            print('ATTENTION: There is no differential regulation in your dataset')
        else:
            print('OmicScope identifies: ' +
                  str(len(self.deps)) + ' deregulations')

    def define_conditions(self):
        """Determine conditions for statistical analysis
        """
        from copy import copy
        if self.ControlGroup is None:
            Control = list(self.pdata.Condition.drop_duplicates())
            Control = sorted(Control)
            self.ControlGroup = Control[0]
            Control.remove(self.ControlGroup)
            self.experimental = Control
        else:
            Conditions = list(self.pdata.Condition.drop_duplicates())
            try:
                Conditions.remove(self.ControlGroup)
                self.experimental = Conditions
            except ValueError:
                raise ValueError(
                    'The user-defined ControlGroup is not present in pdata.')
        self.Conditions = copy([self.ControlGroup]) + copy(self.experimental)

    def expression_function(self):
        """Joins the technical replicates and organizes biological
        conditions.
        """
        from copy import copy

        import numpy as np
        pdata_original = self.pdata.copy()
        assay = self.assay.copy()
        # Create a custom sorting key using the template
        sorting_samples = list(pdata_original.Sample)
        assay = assay[sorting_samples]
        pdata = []
        for i in pdata_original.columns:
            pdata.append(pdata_original[i])
        expression = assay.T
        expression = expression.set_index(pdata).T

        rdata = []
        for i in self.rdata.columns:
            rdata.append(self.rdata[i])
        expression = expression.set_index(rdata).T
        pdata_columns = list(self.pdata.columns)
        pdata_columns = list(set(pdata_columns) - set(['Sample', 'Samples']))
        technical = copy(pdata_original)
        technical = technical.groupby(pdata_columns).agg(','.join).reset_index()
        technical.columns = pdata_columns + ['technical']
        expression = expression.groupby(pdata_columns).mean(numeric_only=True)
        self.Params['Params']['DataHandling_BeforeStat_1'] = 'Biological replicates abundance estimated as mean of technical measures'
        pdata = expression.index.to_frame().reset_index(drop=True)
        Variables = pdata.loc[:, pdata.columns != 'Biological']
        Variables = Variables.apply(lambda x: '-'.join(x.astype(str)), axis=1)
        self.Params['Params']['DataHandling_BeforeStat_2'] = 'Biological replicates was renamed as BioRep_. More info can be found in pdata'
        Sample_name = 'BioRep_' + \
            pdata['Biological'].astype(str) + \
            '.' + Variables
        expression = expression.T
        expression.columns = Sample_name
        rdata = expression.index.to_frame().reset_index(drop=True)
        expression = expression.set_index(rdata.Accession)

        # Filtering data. Protein must be detected at least in one sample per group
        self.Params['Params']['DataHandling_BeforeStat_3'] = 'Filter proteins: minimum 2 samples/condition'
        nanvalues = expression.replace(0, np.nan)
        nanvalues = nanvalues.notna()
        nanvalues.columns = pdata.Condition
        nanvalues = nanvalues.groupby(nanvalues.columns, axis=1).sum()
        nanvalues[nanvalues <= 1] = np.nan
        nanvalues = nanvalues.dropna()
        expression = expression[expression.index.isin(nanvalues.index)]
        rdata = rdata[rdata['Accession'].isin(nanvalues.index)]
        # New pdata
        pdata['Sample'] = Sample_name
        pdata = pdata.merge(technical, how='left', on=pdata_columns)

        self.pdata = pdata

        # New rdata
        self.rdata = rdata.reset_index(drop=True)
        return expression

    def deps(self):
        """Get a dataframe just with differentially regulated entities.
        """
        from copy import copy
        FoldChange_cutoff = self.FoldChange_cutoff
        PValue_cutoff = self.PValue_cutoff
        deps = copy(self.quant_data)
        deps = deps[deps[self.pvalue] < PValue_cutoff]
        deps = deps.loc[(deps['log2(fc)'] <= -FoldChange_cutoff) |
                        (deps['log2(fc)'] >= FoldChange_cutoff)]
        deps = deps[['gene_name', 'Accession', self.pvalue,
                     f'-log10({self.pvalue})', 'log2(fc)']]
        return (deps)

    def savefile(self, Path: str):
        from copy import copy
        import os
        data = copy(self)
        experimental = data.experimental
        string = '-'.join(data.Conditions)
        with open(Path + '/' + string + '.omics', 'w') as f:
            dfAsString = data.quant_data[['gene_name', 'Accession', self.pvalue, 'log2(fc)', 'TotalMean']].to_csv(
                sep='\t', index=False).replace('\r', "")
            file = str("OmicScope" + "\n" +
                    "This file is the output performed by OmicScope pipeline and can be used as input" +
                    " for group comparisons having the controlling group used as used according to OmicScope." +
                    "Please, cite: Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: a Comprehensive Python" +
                    "package designed for Shotgun Proteomics" +
                    '\nControlGroup:' + '\t' + data.ctrl + '\n' +
                    'Experimental:' + '\t' + '\t'.join(experimental) + '\n' +
                    'Statistics:' + '\t' + self.pvalue + '\n' +
                    'Expression:\n' + '-------\n' +
                    dfAsString)
            f.write(file.replace(os.linesep, '\n'))
    __all__ = [
        'bar_ident',
        'heatmap',
        'correlation',
        'DynamicRange',
        'pca',
        'PPInteractions',
        'bar_protein',
        'boxplot_protein',
        'MAplot',
        'k_trend',
        'volcano',
        'normalization_boxplot'
    ]
