__version__ = '1.0.0'
from .Input import *
from .Omicscope import Omicscope
from .EnrichmentAnalysis import *
from .GeneralVisualization import *
from .PeptideLens import *
from .MultipleData import *



if __name__ != "__main__":
    print("OmicScope v " + __version__ + " For help: Insert\n" +
          "If you use  in published research, please cite: 'lallalala'\n" +
          "Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: an Comprehensive Python library for Systems Biology Visualization.\n")
