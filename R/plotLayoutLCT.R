#' Plot the layout of a Latent Class Tree
#'
#'
#' @param resTree Results of a Latent Class Tree analysis in a LCT object
#'
#' @return The layout of a Latent Class Tree
#' @export
plotLayoutLCT <- function(resTree, ...){

  classSizes = round(c(resTree$splitInfo$Ntot,
                       resTree$splitInfo$CppG * resTree$splitInfo$Ntot))

  qgraph::qgraph(cbind(computeEdges(resTree)[,-3], resTree$splitInfo$CppG * resTree$splitInfo$Ntot),
                 layout = igraph::layout_as_tree,
                 labels = classSizes)
}


