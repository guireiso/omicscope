"""   Import data

OmicScope allows user to import data from several Omics Software.
To date, OmicScope is designed to run with output from Progenesis,
PatternLab V, and a generic data.

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
        Methods = ['Progenesis', 'General', 'PatternLab', 'MaxQuant', 'DIA-NN']
        if Method not in Methods:
            raise ValueError("Invalid Method input. Expected one of: %s" % Methods)
        if Method == 'Progenesis':
            from .Progenesis import Input
        elif Method == 'General':
            from .General import Input
        elif Method == 'PatternLab':
            from .PatternLab import Input
        elif Method == 'MaxQuant':
            from .MaxQuant import Input
        elif Method == 'DIA-NN':
            from .DIANN import Input
        data = Input(Table, **kwargs)
        self.Table = Table
        self.rdata = data.rdata
        self.assay = data.assay
        self.pdata = data.pdata
        self.Conditions = data.Conditions
