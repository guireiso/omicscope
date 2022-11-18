"""  OmicScope Module

OmicScope allows user to import data and perform statistical analysis
for differential gene/protein expression. This is the key module of 
all OmicScope workflow, being the input for all other modules.

"""

from .Input import *


class Omicscope(Input):
    def __init__(self, Table, ControlGroup, Method, FoldChange_cutoff=0,
                 statistics = 'pAdjusted', pdata = '', ExperimentalDesign = 'static',
                 PValue_cutoff=0.05, logTransformed=False, ExcludeKeratins=True, **kwargs):
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
            ControlGroup (str): Control group.
            Method (str): Algorithm/software used to perform quantitative
            FoldChange_cutoff (float, optional): log2(Fold-Change) threshold
            for differentially regulated entities against ControlGroup.
            Defaults to 0.
            PValue_cutoff (float, optional): P-value threshold for differential
            expression analysis . Defaults to 0.05.
            logTransformed (bool, optional): Was entity abundance previously 
            log-normalized?. Defaults False.
            ExcludeKeratins (bool, optional): Drop keratins from dataset.
            Defaults to True.
        """

        super().__init__(Table, ControlGroup, Method, **kwargs)
        self.PValue_cutoff = PValue_cutoff
        self.FoldChange_cutoff = FoldChange_cutoff
        self.PValue_cutoff = PValue_cutoff
        self.logTransformed = logTransformed
        self.ExcludeKeratins = ExcludeKeratins
        if pdata != '':
            import pandas as pd
            try:
                self.pdata = pd.read_excel(pdata)
            except:
                self.pdata = pd.read_csv(pdata)

        # Construct pivot-table considering technical and biological replicates
        self.expression = self.expression()

        #  Has user already performed statistical analyses?
        statistics = [f'Anova (p)', 'pvalue', 'p-value', 'q-value', 'q Value', 'qvalue',
                    'padj', 'p-adjusted', 'p-adj']
        if True in self.rdata.columns.isin(statistics):
            from .Stats.Performed_Stat import imported_stat
            self.quant_data = imported_stat(self, statistics)
            print('User already performed statistical analysis')
        elif ExperimentalDesign == 'static':  # OmicScope perform statistics
            from .Stats.Statistic_Module import perform_static_stat
            self.quant_data = perform_static_stat(self)
            print('OmicScope performed statistical analysis (Static workflow)')
        elif ExperimentalDesign == 'longitudinal':
            from .Stats.Statistic_Module import perform_longitudinal_stat
            self.quant_data = perform_longitudinal_stat(self)
            print('OmicScope performed statistical analysis (Longitudinal workflow)')

        self.deps = self.deps()

    def expression(self):
        """Joins the technical replicates and organizes biological
        conditions.
        """
        from copy import copy
        pdata = []
        for i in self.pdata.columns:
            pdata.append(self.pdata[i])
        expression = self.assay.T
        expression = expression.set_index(pdata).T

        rdata = []
        for i in self.rdata.columns:
            rdata.append(self.rdata[i])
        expression = expression.set_index(rdata)

        # Melt expression data to get all values for each sample and each gene
        expression = expression.melt(ignore_index=False).reset_index()

        # Getting the mean of technical replicates for each biological replicate
        var = list(self.rdata.columns)
        pdata_columns = list(self.pdata.columns)
        pdata_columns = list(set(pdata_columns) - set(['Sample', 'Samples', 'TechRep']))
        for i in pdata_columns:
            var.append(i)

        expression = expression.groupby(var).mean(numeric_only= True).reset_index().rename(
            columns={'value': 'abundance'})

        expression['Sample'] = 'BioRep' + \
            expression['Biological'].astype(str) + \
            '.' + expression['Condition']
        
        # New pdata
        pdata_columns.append('Sample')
        pdata = expression[pdata_columns]
        self.new_pdata = pdata

        # New rdata 
        rdata = expression[copy(self.rdata.columns)]
        self.new_rdata = rdata
        #  data with just biological replicates abundance
        expression = expression[['Sample', 'Accession', 'abundance']]
        expression = expression.pivot(columns='Sample', index='Accession', values='abundance')
        return(expression)

    def deps(self):
        """Get a dataframe just with differentially regulated entities.
        """
        FoldChange_cutoff = self.FoldChange_cutoff
        PValue_cutoff = self.PValue_cutoff
        deps = self.quant_data
        deps = deps[deps['pvalue'] < PValue_cutoff]
        deps = deps.loc[(deps['log2(fc)'] <= -FoldChange_cutoff) |
                        (deps['log2(fc)'] >= FoldChange_cutoff)]
        deps = deps[['gene_name', 'Accession', 'pvalue',
                     '-log10(p)', 'log2(fc)']]
        return(deps)

    def savefile(self, Path: str):
        from copy import copy
        data = copy(self)
        experimental = data.experimental
        string = '-'.join(data.Conditions)
        with open(Path + '/' + string + '.omics', 'w') as f:
            dfAsString = data.quant_data[['gene_name', 'Accession', 'pvalue', 'log2(fc)', 'TotalMean']].to_csv(
                sep='\t', index=False)
            f.write("Omics v1.0.0" + "\n" +
                  "This file is the output performed by OmicScope pipeline and can be used as input" +
                  " for group comparisons having the controling group used as used according to OmicScope." +
                  "Please, cite: Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: a Comprehensive Python package designed for Shotgun Proteomics" +
                  '\nControlGroup:' + '\t' + data.ctrl + '\n' +
                  'Experimental:' + '\t' + '\t'.join(experimental) + '\n' +
                  'Expression:\n' + '-------\n' +
                  dfAsString)
