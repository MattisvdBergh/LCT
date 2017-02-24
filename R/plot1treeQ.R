#' @title Plot a growing tree
#'
#' @param resTree A Latent Class Tree object
#' @param edg Edges of the layout of the tree. Can be computed with \code{computeEdges}.
#' @param dirLCTsplits list of directories of an figure for every split. This is returned by \code{plotLCTsplits}.
#' @param name Name of the to be created pdf .
#' @param vsize heighth of the nodes
#' @param vsize2 width of the nodes
#'
#' @description 1 pdf file with a growing tree.
#'
#' @export
plot1treeQ = function(resTree,
                     edg,
                     dirLCTsplits,
                     name = "Tree.pdf",
                     vsize = 8,
                     vsize2 = 0.5){

  exprPlots = list()

  plotsLogical = edg[,2] %in% edg[,1]

  exprPlots[c(TRUE, plotsLogical)] = dirLCTsplits[[2]]

  labelsLogical = c(FALSE, !plotsLogical)
  allLabels = round(c(resTree$splitInfo$Ntot,
                      resTree$splitInfo$CppG * resTree$splitInfo$Ntot))
  labels = labelsLogical
  labels[labelsLogical] = allLabels[labelsLogical]
  labels[!labelsLogical] = ""

  vsizeInf = rep(vsize, length(labelsLogical))
  vsizeInf[labelsLogical] = 3

  pdf(name, height = 15, width = 10)
  for(idxLevels in unique(edg[,3])[-1]){
    idxRows = edg[,3] <= idxLevels

    edgTemp = edg[idxRows,]
    idxNsplits = 1:max(edgTemp[,2])

    qgraph(edgTemp[,-3], layout = layout_as_tree,
           subplots = exprPlots[idxNsplits],
           vsize = vsizeInf[idxNsplits], vsize2 = vsize2,
           labels = labels[idxNsplits], borders = labelsLogical[idxNsplits],
           mar = c(5,4,4,2))
  }
  dev.off()
}
