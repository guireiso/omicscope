__version__ = '1.0.0'
from .Input import *
from .Omicscope import Omicscope
from .EnrichmentAnalysis import *
from .GeneralVisualization import (DynamicRange, MAplot, bar_ident,
                                   bar_protein, boxplot_protein, heatmap, pca,
                                   correlation, volcano, volcano_qvalue )
from .PeptideLens import *
from .MultipleData.circos import circos
from .MultipleData.multiples import multiples

__version__ = get_versions()['version']
del get_versions


if __name__ != "__main__":
    print("OmicScope v " + __version__ + " For help: Insert\n" +
          "If you use  in published research, please cite: 'lallalala'\n" +
          "Reis-de-Oliveira G, Martins-de-Souza D. OmicScope: an Comprehensive Python library for Systems Biology Visualization.\n")
