import warnings

import os
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


class Omicscope_Snapshot():
    from .GeneralVisualization import PPInteractions
    from .GeneralVisualization import bar_ident
    from .GeneralVisualization import volcano

    def __init__(self, Table,
                 pvalue='pAdjusted', PValue_cutoff=0.05,
                 FoldChange_cutoff=0):
        """ Snapshot
        Snapshot is a module that enables users to import pre-analyzed data into OmicScope quickly.
          To use Snapshot, the user must upload a CSV or Excel file that includes the following
          information in the first and second rows/cells: 'ControlGroup: INSERT_HERE_YOUR_CONTROL'
          and 'Experimental: INSERT_HERE_YOUR_EXPERIMENTAL_GROUPS_SEPARATED_BY_COMMAS', respectively.
          The third row should contain the specific data, including columns for 'Accession', 'gene_name',
          'log2(fc)', and either 'pvalue' or 'pAdjusted'.
        Args:
            Table (str): Path to quantitative data.
            pvalue (str, optional): Statistical parameter to take into account
            to consider entities differentially regulated. Defaults to 'pAdjusted'.
            PValue_cutoff (float, optional): Statistical cutoff. Defaults to 0.05.
            FoldChange_cutoff (float, optional): Difference cutoff. Defaults to 0.0..
        """
        self.Table = Table
        self.PValue_cutoff = PValue_cutoff
        self.FoldChange_cutoff = FoldChange_cutoff
        self.pvalue = pvalue
        self.pdata = None
        self.Params = None
        self.rdata = None
        # Define Conditions and Quant_data
        self.define_conditions_quantdata()

        # Check for correct dataframe build
        if not all(item in self.quant_data.columns for item in ['Accession', 'gene_name', 'log2(fc)',
                                                                pvalue]):
            raise IndexError(f"Maybe some columns are missing in your input (Accession, {pvalue}, gene_name, log2(fc))")
        pvalues = ['pvalue', 'pAdjusted']
        if pvalue not in pvalues:
            raise ValueError("Invalid pvalue specification. Expected one of: %s" % pvalues)

        self.quant_data[f'-log10({pvalue})'] = -np.log10(self.quant_data[pvalue])

        if len(self.deps) == 0:
            print('ATTENTION: There is no differential regulation in your dataset')
        else:
            print('OmicScope identifies: ' + str(len(self.deps)) + ' deregulations')

    def define_conditions_quantdata(self):
        """Extracting conditions and define quant_data
        The study conditions must be specified in the first two rows of the imported file,
         of which the first one must be the ControlGroup and the second the experimental groups.
        """
        try:
            quant_data = pd.read_csv(self.Table, header=2)
            quant_data['TotalMean'] = 1
            self.quant_data = quant_data
            deps = quant_data[quant_data[self.pvalue] <= self.PValue_cutoff]
            self.deps = deps
            Conditions = pd.read_csv(self.Table,
                                     nrows=2, names=['Conditions'], sep='\t')
            Conditions = self.Conditions = Conditions.Conditions.str.split(':').str[1]
            self.ctrl = self.ControlGroup = Conditions[0]
            self.experimental = Conditions[1].split(',')
        except UnicodeDecodeError:
            quant_data = pd.read_excel(self.Table, header=2)
            quant_data['TotalMean'] = 1
            self.quant_data = quant_data
            deps = quant_data[quant_data[self.pvalue] <= self.PValue_cutoff]
            self.deps = deps

            Conditions = pd.read_excel(self.Table,
                                       nrows=2)
            Conditions.iloc[1, :] = Conditions.iloc[0, :]
            Conditions.iloc[1, :] = Conditions.columns
            Conditions = Conditions.dropna(axis=1)
            Conditions.columns = ['Conditions']
            Conditions = Conditions.sort_values(by='Conditions')
            Conditions = self.Conditions = Conditions.Conditions.str.split(':').str[1]
            self.ctrl = self.ControlGroup = Conditions[1]
            self.experimental = Conditions[0].split(',')

    def savefile(self, Path: str):
        from copy import copy
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
        'PPInteractions',
        'volcano'
    ]
