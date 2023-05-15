list.of.packages <- c("jsonlite", "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(jsonlite)
suppressPackageStartupMessages(library("circlize"))
library(circlize)

pydata = jsonlite::read_json('circlize.json')

#define matrix and colmat
matrix = as.data.frame(fromJSON(pydata$matrix))
colmat = as.data.frame(fromJSON(pydata$colmat))
gene_names = fromJSON(pydata$gene_names)
row.names(matrix) = gene_names
row.names(colmat) = gene_names
# define other parameters
grid.col = fromJSON(pydata$colors)
name.grid = fromJSON(pydata$labels)
wid = pydata$width
hei = pydata$height
save = pydata$save
vector = pydata$vector

names(grid.col) = name.grid

circos.plot <- function(mat, colmat,grid.col, save, vector, y0 = -1.3, y = 0){
    mat1 = as.matrix(mat)
    colmat = as.matrix(colmat)
    mat1 = mat1
    mat1[mat1==0.5] <- 1
    
    if (!(is.null(save))){
      if (vector){
        string = svg(paste0(save, "_my_plot.svg"), width = 10, height = 10)
      }else{
        string = png(paste0(save, "_my_plot.png"), width = wid, height = hei, res = 300)
      }
    }else{
      string = NULL
    }
        
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


circos.plot(matrix, colmat, grid.col, save, vector, y0 = -0.2, y= -2)
