#' Compute the layout of a Latent Class Tree
#'
#'
#' @param resTree A Latent Class Tree object
#'
#' @details The \code{computeEdges} function constructs the layout of a LCT, but only the splits.
#' This can be used by \code{igraph} and \code{qgraph} to make a figure without redundant nodes.
#'
#' @return Returns a matrix with in the first column the parent nodes (multiple entries for every split),
#'  the second column the child nodes and in the third column the level of the tree.
#' @export
computeEdgesSplits = function(resTree){

  Splits = resTree$treeSetup$Splits
  cleanNames = resTree$treeSetup$cleanNames

  namesSplitPoints = list()
  E = Enames = matrix(,0,2)

  for(idxLevels in 1:(length(Splits) - 1)){
    namesSplitPoints[[idxLevels]] = cleanNames[[idxLevels]][Splits[[idxLevels]] > 1]
  }

  counterSplitPoints = list(1)
  length(counterSplitPoints) = length(namesSplitPoints)
  counter = 2
  for(idxLevels in 2:length(namesSplitPoints)){
    for(idxClasses in 1:length(namesSplitPoints[[idxLevels]])){
      counterSplitPoints[[idxLevels]][idxClasses] = counter
      counter = counter + 1
    }}


  E = matrix( ,0, 3)
  counter1 = 2
  for(idxLevels in 2:length(namesSplitPoints)){
    for(idxClasses in 1:length(namesSplitPoints[[idxLevels]])){
      cc = namesSplitPoints[[idxLevels]][idxClasses]
      whichCount = substr(cc, 1, nchar(cc) - 1) == namesSplitPoints[[idxLevels - 1]]
      namesSplitPoints
      EtoBe = matrix(c(
        counterSplitPoints[[idxLevels - 1]][whichCount],
        counterSplitPoints[[idxLevels]][idxClasses],
        idxLevels),
        ncol = 3)
      E = rbind(E, EtoBe)
    }}
  return(E)
}
