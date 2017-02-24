#' Compute the layout of a Latent Class Tree
#'
#' @param resTree A Latent Class Tree object
#'
#' @details The \code{computeEdges} function constructs the layout of a LCT to be used by \code{igraph} and \code{qgraph}.
#'
#' @return Returns a matrix with in the first column the parent nodes (multiple entries for every split),
#'  the second column the child nodes and in the third column the level of the tree
#' @export
computeEdges = function(resTree){

  cleanNames = resTree$treeSetup$cleanNames
  Splits = resTree$treeSetup$Splits

  counter1 = counter2 = 1
  E = Enames = matrix(,0,3)
  for(idxLevels in 1:(length(cleanNames) - 1)){
    level1Temp = cleanNames[[idxLevels]]
    pSplits = nchar(level1Temp) == max(nchar(level1Temp))
    level2Temp = cleanNames[[idxLevels + 1]]
    rSplits = Splits[[idxLevels]]>1

    for(idxClasses in 1:length(level1Temp)){
      if(pSplits[idxClasses]){
        if(rSplits[idxClasses]){

          splitClassesLogical = grepl(level1Temp[idxClasses], level2Temp)
          counter2 = counter2 + cumsum(splitClassesLogical)[splitClassesLogical]

          EnamestoBe = cbind(level1Temp[idxClasses], level2Temp[splitClassesLogical], idxLevels)
          EtoBe = cbind(counter1, counter2, idxLevels)

          counter1 = counter1 + 1
          counter2 = max(EtoBe[,2])

          E = rbind(E,EtoBe)
          Enames = rbind(Enames, EnamestoBe)

        } else {counter1 = counter1 + 1}
      }}}
  return(E)
}
