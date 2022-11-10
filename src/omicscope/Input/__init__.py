"""   Import data

OmicScope allows user to import data from several Omics Software.
To date, OmicScope is designed to run with output from Progenesis,
PatternLab V, and a generic data.

Since outputs can largely vary among software, user must specify
from what software the data came from.

"""


class Input():
    def __init__(self, Table, ControlGroup, Method, **kwargs):
        """  Input data from several software sources

        Args:
            Table (str): Path to quantitative data.
            ControlGroup (str): Control group.
            Method (str): algorithm/software used to perform quantitative
            proteomics
        """
        if Method == 'Progenesis':
            from .Progenesis import Input
        elif Method == 'General':
            from .General import Input
        elif Method == 'PatternLab':
            from .PatternLab import Input
        elif Method == 'MaxQuant':
            from .MaxQuant import Input
        data = Input(Table, ControlGroup, **kwargs)
        self.Table = Table
        self.ctrl = ControlGroup
        self.rdata = data.rdata
        self.assay = data.assay
        self.pdata = data.pdata
        self.Conditions = data.Conditions
        self.experimental = data.experimental
        return(Input)
