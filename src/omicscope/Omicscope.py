"""  OmicScope Module

OmicScope allows user to import data and perform statistical analysis
for differential gene/protein expression. This is the key module of
all OmicScope workflow, being the input for all other modules.

"""

from .Input import Input


class Omicscope(Input):
    from .GeneralVisualization import DynamicRange
    from .GeneralVisualization import MAplot
    from .GeneralVisualization import bar_ident
    from .GeneralVisualization import bar_protein
    from .GeneralVisualization import bigtrend
    from .GeneralVisualization import boxplot_protein
    from .GeneralVisualization import color_scheme
    from .GeneralVisualization import correlation
    from .GeneralVisualization import find_k
    from .GeneralVisualization import heatmap
    from .GeneralVisualization import pca
    from .GeneralVisualization import volcano
    from .GeneralVisualization import volcano_2cond
    from .GeneralVisualization import volcano_Multicond

    def __init__(self, Table, Method, ControlGroup, ExperimentalDesign='static',
                 pvalue='pAdjusted', PValue_cutoff=0.05,
                 FoldChange_cutoff=0, logTransformed=False, ExcludeKeratins=True,
                 degrees_of_freedom=2, **kwargs):
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
            FoldChange_cutoff (float, optional): Difference cutoff. Defaults to 0.0.
            logTransformed (bool, optional): Abundance values were previously log-transformed. Defaults to False.
            ExcludeKeratins (bool, optional): Keratins proteins is excluded. Defaults to True.
            degrees_of_freedom (int, optional): Degrees of freedom used to run longitudinal analysis. Defaults to 2.
        """
        import pandas as pd
        super().__init__(Table, Method=Method, **kwargs)
        self.PValue_cutoff = PValue_cutoff
        self.FoldChange_cutoff = FoldChange_cutoff
        self.logTransformed = logTransformed
        self.ExcludeKeratins = ExcludeKeratins
        self.pvalue = pvalue
        self.ControlGroup = ControlGroup
        pvalues = ['pvalue', 'pAdjusted', 'pTukey']
        if pvalue not in pvalues:
            raise ValueError("Invalid pvalue specification. Expected one of: %s" % pvalues)
        if 'pdata' in kwargs:
            # If pdata was assigned by user, OmicScope read
            # excel or csv frames.
            try:
                self.pdata = pd.read_excel(kwargs['pdata'])
            except ValueError:
                self.pdata = pd.read_csv(kwargs['pdata'])
        self.define_conditions()

        self.ctrl = self.ControlGroup
        #  Has user already performed statistical analyses?
        statistics = ['Anova (p)', 'pvalue', 'p-value', 'q-value', 'q Value', 'qvalue',
                      'padj', 'p-adjusted', 'p-adj']
        if True in self.rdata.columns.isin(statistics):
            from .Stats.Performed_Stat import imported_stat
            self.quant_data = imported_stat(self, statistics)
            self.pdata['Sample'] = self.pdata['Sample']+'.'+self.pdata['Condition']
            print('User already performed statistical analysis')

        elif ExperimentalDesign == 'static':  # OmicScope perform statistics
            from .Stats.Statistic_Module import perform_static_stat

            # Construct pivot-table considering technical and biological replicates
            self.expression = self.expression()
            # perform stat
            self.quant_data = perform_static_stat(self)
            print('OmicScope performed statistical analysis (Static workflow)')
        elif ExperimentalDesign == 'longitudinal':
            from .Stats.Statistic_Module import perform_longitudinal_stat
            self.expression = self.expression()
            self.degrees_of_freedom = degrees_of_freedom
            self.quant_data = perform_longitudinal_stat(self)
            print('OmicScope performed statistical analysis (Longitudinal workflow)')

        self.deps = self.deps()

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
            Conditions.remove(self.ControlGroup)
            self.experimental = Conditions
        self.Conditions = copy([self.ControlGroup]) + copy(self.experimental)

    def expression(self):
        """Joins the technical replicates and organizes biological
        conditions.
        """
        pdata = []
        for i in self.pdata.columns:
            pdata.append(self.pdata[i])
        expression = self.assay.T
        expression = expression.set_index(pdata).T

        rdata = []
        for i in self.rdata.columns:
            rdata.append(self.rdata[i])
        expression = expression.set_index(rdata).T
        pdata_columns = list(self.pdata.columns)
        pdata_columns = list(set(pdata_columns) - set(['Sample', 'Samples']))
        expression = expression.groupby(pdata_columns).mean(numeric_only=True)
        pdata = expression.index.to_frame().reset_index(drop=True)
        Variables = pdata.loc[:, pdata.columns != 'Biological']
        Variables = Variables.apply(lambda x: '-'.join(x.astype(str)), axis=1)
        Sample_name = 'BioRep_' + \
            pdata['Biological'].astype(str) + \
            '.' + Variables
        expression = expression.T
        expression.columns = Sample_name
        rdata = expression.index.to_frame().reset_index(drop=True)
        expression = expression.set_index(rdata.Accession)

        # New pdata
        pdata['Sample'] = Sample_name
        self.pdata = pdata

        # New rdata
        self.rdata = rdata
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
        data = copy(self)
        experimental = data.experimental
        string = '-'.join(data.Conditions)
        with open(Path + '/' + string + '.omics', 'w') as f:
            dfAsString = data.quant_data[['gene_name', 'Accession', self.pvalue, 'log2(fc)', 'TotalMean']].to_csv(
                sep='\t', index=False)
            f.write("Omics v1.0.0" + "\n" +
                    "This file is the output performed by OmicScope pipeline and can be used as input" +
                    " for group comparisons having the controlling group used as used according to OmicScope." +
                    "Please, cite: Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: a Comprehensive Python" +
                    "package designed for Shotgun Proteomics" +
                    '\nControlGroup:' + '\t' + data.ctrl + '\n' +
                    'Experimental:' + '\t' + '\t'.join(experimental) + '\n' +
                    'Statistics:' + '\t' + self.pvalue + '\n' +
                    'Expression:\n' + '-------\n' +
                    dfAsString)
    __all__ = [
        'bar_ident',
        'heatmap',
        'correlation',
        'DynamicRange',
        'pca',
        'color_scheme',
        'bar_protein',
        'boxplot_protein',
        'MAplot',
        'find_k',
        'bigtrend',
        'volcano',
        'volcano_2cond',
        'volcano_Multicond'
    ]
