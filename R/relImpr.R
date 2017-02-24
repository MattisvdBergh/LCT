#' Compute relative improvement at the root of the tree
#'
#' Compute relative improvement by adding classes at the root of the tree
#'
#' @param x treeobject
#' @param Criterion relative improvement for which criterion, default is the log likelihood, other options are AIC and BIC
#'
#' @return None
#' @export
relImpr = function(x, criterion = "LL"){

  allImprovement = matrix(, nrow = length(x$splitInfo$IC[,criterion[1]][[1]]) - 2,
                          ncol = length(criterion))

  for(idxCrit in seq_along(criterion)){
    nClass = length(x$splitInfo$IC[,criterion[idxCrit]][[1]])
    improvement = numeric()
    counter = 1
    for(i in 3:nClass){
      improvement[counter] = (x$splitInfo$IC[,criterion[idxCrit]][[1]][i] -
                                x$splitInfo$IC[,criterion[idxCrit]][[1]][i - 1])/
        (x$splitInfo$IC[,criterion[idxCrit]][[1]][2] -
           x$splitInfo$IC[,criterion[idxCrit]][[1]][1])
      counter = counter + 1
    }
    allImprovement[,idxCrit] = improvement
  }
  colnames(allImprovement) = criterion
  return(allImprovement)
}
