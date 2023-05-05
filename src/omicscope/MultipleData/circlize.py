from copy import copy

import pandas as pd


def deps(self, pvalue=0.05):
    data = copy(self)
    groups = data.groups
    regulation = data.original
    regulation = [x[x[data.pvalue] < pvalue] for x in regulation]
    regulation = [x[['gene_name', 'log2(fc)']] for x in regulation]
    regulation = [x.rename(columns={'log2(fc)': y}) for x, y in zip(regulation, groups)]
    regulation = pd.concat(regulation)
    regulation = regulation.groupby('gene_name').sum().reset_index()
    return regulation


def enrichment_filtering(self, Term_enriched):
    data = copy(self.enrichment)
    data = [x[x['Term'].str.contains(Term_enriched)] for x in data]
    data = pd.concat(data)
    deps = list(data['Genes'])
    deps = [x.replace("'", '') for x in deps]
    deps = [x.replace("[", '') for x in deps]
    deps = [x.replace("]", '') for x in deps]
    deps = [x.split(', ') for x in deps]
    deps = sum(deps, [])
    deps = list(set(deps))
    return deps


def deps_matrix(df):
    deps = copy(df)
    deps[deps > 0] = 1
    deps[deps < 0] = 0.5
    deps[deps.isna()] = 0
    return deps


def color_matrix(df, colors):
    colmat = copy(df)
    for i, z in zip(range(0, len(colmat.columns)), colors):
        colmat.iloc[:, i] = z
    return colmat


def circlize(matrix, colmat, colors, labels, width=3000, height=3000,
             save=None, vector=True):
    import rpy2.robjects as robjects
    from IPython.display import Image
    from IPython.display import display
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.lib import grdevices
    with localconverter(robjects.default_converter + pandas2ri.converter):
        matrix_2 = robjects.conversion.py2rpy(matrix)
    with localconverter(robjects.default_converter + pandas2ri.converter):
        colmat_2 = robjects.conversion.py2rpy(colmat)
    robjects.globalenv["matrix"] = matrix_2
    robjects.globalenv["colmat"] = colmat_2
    robjects.globalenv["grid.col"] = robjects.StrVector(colors)
    robjects.globalenv["name.grid"] = robjects.StrVector(labels)
    robjects.globalenv["wid"] = width
    robjects.globalenv["hei"] = height
    if save is not None:
        if vector is True:
            string = f'svg("{save}my_plot.svg", width = 10, height = 10)'
        else:
            string = f'png("{save}my_plot.png", width = wid, height = hei, res = 300)'
    else:
        string = '\n'
    with grdevices.render_to_bytesio(grdevices.jpeg, width=3000, height=3000, res=300) as img:
        robjects.r(
            '''
    suppressPackageStartupMessages(library(circlize))
    library(circlize)
    names(grid.col) = name.grid
    circos.plot <- function(mat, colmat,grid.col, y0 = -1.3, y = 0){
    mat1 = as.matrix(mat)
    colmat = as.matrix(colmat)
    mat1 = mat1
    mat1[mat1==0.5] <- 1
      ''' +
            string
            +
            '''
    circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)
  cdm_res = chordDiagram(t(mat1),
                         annotationTrack = c('grid'), grid.col = grid.col,
                         annotationTrackHeight = c(0.05, 0.02),
                         big.gap = 5, small.gap = 1,
                         preAllocateTracks = list(track.height = 0.1),
                         directional = TRUE,
                         link.target.prop = FALSE)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
    sector.name = get.cell.meta.data('sector.index')
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = 'clockwise', niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA)
  col_fun = colorRamp2(c(0.5, 0.75, 1), c('#293375', 'white', '#ba1616'))
  y1 = y0
  y2 = y

  for(i in seq_len(nrow(cdm_res))) {
    if(cdm_res$value1[i] != 0) {
      circos.rect(cdm_res[i, 'x2'], y1, cdm_res[i, 'x2'] - abs(cdm_res[i, 'value1']), y1 + (y2-y1)*0.3,
                  col = col_fun(t(mat)[cdm_res$rn[i], cdm_res$cn[i]]),
                  border = col_fun(t(mat)[cdm_res$rn[i], cdm_res$cn[i]]),
                  sector.index = cdm_res$cn[i], track.index = 1, )}}
    dev.off()
  }

  circos.plot(matrix, colmat, grid.col, y0 = -0.2, y= -2)

  ''')
    display(Image(data=img.getvalue(), format='png', embed=True))
