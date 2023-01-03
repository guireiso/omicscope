"""  OmicScope Module

OmicScope allows user to import data and perform statistical analysis
for differential gene/protein expression. This is the key module of 
all OmicScope workflow, being the input for all other modules.

"""

from ..Input import *


class OmicScope_Peptide(Input):
    def __init__(self, Table, ControlGroup, Method, FoldChange_cutoff=0,
                 PValue_cutoff=0.05, logTransformed=False, ExcludeKeratins=True,
                 **kwargs):
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
        self.FoldChange_cutoff = FoldChange_cutoff
        self.PValue_cutoff = PValue_cutoff
        self.logTransformed = logTransformed
        self.ExcludeKeratins = ExcludeKeratins
        # Construct pivot-table considering technical and biological replicates
        self.expression = self.expression()
        from ..Stats.Statistic_Module import perform_stat
        quant_data = perform_stat(self)
        columns = list(quant_data.columns)
        columns.remove('Modifications')
        quant_data['Modifications'] = quant_data['Modifications'].astype(str)
        quant_data['Accession'] = quant_data['Accession'].str.split('.').str[-1]
        quant_data = quant_data.groupby(columns)['Modifications'].apply('|'.join).reset_index()
        quant_data['gene_name'] = quant_data.Description.str.split('GN=').str[-1]
        quant_data['gene_name'] = quant_data['gene_name'].str.split(' ').str[0]
        import numpy as np
        quant_data['gene_name'] = quant_data['gene_name'].replace(np.nan, '')
        quant_data['gene_name'] = quant_data['gene_name']+'-'+quant_data['Sequence']
        self.quant_data = quant_data

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
        expression = expression.set_index(rdata)

        # Melt expression data to get all values for each sample and each gene
        expression = expression.melt(ignore_index=False).reset_index()

        # Getting the mean of technical replicates for each biological replicate
        expression = expression.groupby(['Sequence', 'Accession', 'Condition', 'Biological']).mean().reset_index().rename(
            columns={'value': 'abundance'})

        expression['Sample'] = 'BioRep' + \
            expression['Biological'].astype(str) + \
            '.' + expression['Condition']

        #  data with just biological replicates abundance
        expression['Accession'] = expression['Sequence'] + '.' + expression['Accession']
        expression = expression[['Sample', 'Accession', 'abundance']]
        # Sum peptides intensities
        expression = expression.groupby(['Sample', 'Accession']).sum().reset_index()
        expression = expression.pivot(columns='Sample', index='Accession', values='abundance')

        # Adjusting rdata Accessions
        rdata = self.rdata
        rdata['Accession'] = rdata['Sequence'] + '.' + rdata['Accession']
        self.rdata = rdata[['Sequence', 'Modifications', 'Accession', 'Description']]
        return (expression)
