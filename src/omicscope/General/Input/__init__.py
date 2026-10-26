"""   Import data

OmicScope allows user to import data from several Omics Software.
To date, OmicScope is designed to run with output from Progenesis,
PatternLab V, MaxQuant, DIA-NN, Spectronaut, FragPipe, Proteome Discoverer,
RPKM tables (transcriptomics) and a generic data.

Since outputs can largely vary among software, user must specify
from what software the data came from.

"""


class Input():
    def __init__(self, Table, Method, **kwargs):
        """  Input data from several software sources

        Args:
            Table (str): Path to quantitative data.
            Method (str): algorithm/software used to perform quantitative
            proteomics
        """
        Methods = ['Progenesis', 'General', 'PatternLab',
                   'MaxQuant', 'DIA-NN', 'Spectronaut', 'FragPipe',
                   'ProteomeDiscoverer', 'RPKM']
        if not hasattr(self, 'Params'):
            self.Params = {'Params': {}}
        if Method not in Methods:
            raise ValueError("Invalid Method input. Expected one of: %s" % Methods)
        if Method == 'Progenesis':
            from .Progenesis import Input
            self.Params['Params']['ImportMethod'] = 'Progenesis'
            self.Params['Params']['ImportWarning_1'] = 'Quantification using Normalized Abundance'
        if Method == 'ProteomeDiscoverer':
            from .ProteomeDiscoverer import Input
            self.Params['Params']['ImportMethod'] = 'ProteomeDiscoverer'
            self.Params['Params']['ImportWarning_1'] = 'Quantification using Normalized Abundance'
        elif Method == 'General':
            from .General import Input
            self.Params['Params']['ImportMethod'] = 'General'
        elif Method == 'PatternLab':
            from .PatternLab import Input
            self.Params['Params']['ImportMethod'] = 'PatternLab'
            self.Params['Params']['ImportWarning_1'] = 'Drop reverse and contaminant proteins'
            self.Params['Params']['ImportWarning_2'] = 'Quantification using XIC-normalized abundances'
        elif Method == 'MaxQuant':
            from .MaxQuant import Input
            self.Params['Params']['ImportMethod'] = 'MaxQuant'
            self.Params['Params']['ImportWarning_1'] = 'Quantification using LFQ intensities'
            self.Params['Params']['ImportWarning_2'] = 'Drop REV_ and CON_ proteins'
            self.Params['Params']['ImportWarning_3'] = 'Drop proteins with Andromeda Score <= 0'
        elif Method == 'DIA-NN':
            from .DIANN import Input
            self.Params['Params']['ImportMethod'] = 'DIA-NN'
            self.Params['Params']['ImportWarning_1'] = 'Quantification using PG.MaxLFQ (main report) or report.pg_matrix values'
            self.Params['Params']['ImportWarning_2'] = 'Drop decoy and contaminant protein groups; main report filtered at PG.Q.Value <= 0.01'
        elif Method == 'Spectronaut':
            from .Spectronaut import Input
            self.Params['Params']['ImportMethod'] = 'Spectronaut'
            self.Params['Params']['ImportWarning_1'] = 'Quantification using PG.Quantity'
            self.Params['Params']['ImportWarning_2'] = 'Drop decoy and contaminant protein groups'
        elif Method == 'FragPipe':
            from .FragPipe import Input
            self.Params['Params']['ImportMethod'] = 'FragPipe'
            self.Params['Params']['ImportWarning_1'] = 'Quantification using MaxLFQ intensities (TMT-Integrator abundances for isobaric data)'
            self.Params['Params']['ImportWarning_2'] = 'Drop decoy and contaminant proteins'
        elif Method=='RPKM':
            from .RPKM import Input
            self.Params['Params']['ImportMethod'] = 'RPKM - Transcriptomics'
            self.Params['Params']['ImportWarning_1'] = 'LogTransformation using log(RPKM+1)'
        data = Input(Table, **kwargs)
        if Method == 'MaxQuant':
            self.Params['Params']['ImportWarning_1'] = f'Quantification using {data.quant_strategy} columns'
        self.Table = Table
        self.rdata = data.rdata
        self.assay = data.assay
        self.pdata = data.pdata
        self.Conditions = data.Conditions
