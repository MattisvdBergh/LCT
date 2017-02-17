#' Print tree results to the terminal
#'
#' These function display the results of a Latent Class Tree analysis
#'
#' @param dataDir directory of the data or a dataframe with the data
#' @param LG directory of the Latent GOLD executable
#' @param LGS directory of Latent GOLD syntax for a model with 1- and 2-class splits
#' @param decreasing Whether the ordering of classes should be decreasing or not. Defaults to TRUE.
#' @param maxClassSplit1 Maximum size of the first split of the tree. Will be assessed with the criterion given in stopCriterium. Defaults to two.
#' @param maxClassSplit2 Maximum size of each split after the first split of the tree. Defaults to two.
#' @param stopCriterium Criterium to decide on a split. Can be "LL" (logLikelihood), "AIC" or "BIC".
#' @param resultsName Name of a folder which will be created in the working directory and contains all results by Latent GOLD.
#' @param minSampleSize Minimum sample size of a class. If this is below 1, a probability of the total sample size is used.
#' @param itemNames The names of the indicators. If this is not given, the rownames of the datafile will be used.
#' @param nKeepVariables Number of variables to be kept if one wants to explore the results with external variables.
#' @param weight Name of the variable with the weights. When all records are unique observations, this should be one for every observation.
#' @param measurementLevels A character vector being either ordinal or continuous to indicate the measurement level of each variable. It is required when LGS is specified.
#'
#' @return None
#' @export
print.LCT <- function(x, ...){
  cat("tree object","\n\n")
  for(idxLevels in 1:(length(x$Names.clean) - 1)){
    cat("level",idxLevels,": ", x$Names.clean[[idxLevels]], "\n")
  }
  cat("Class Sizes final level:\n")
  print(x$ClassppGlobal[x$finalClasses])
  qgraph::qgraph(computeEdges(x), layout = igraph::layout_as_tree, labels = rownames(x$ClassProportions))
}

makeCleanNames = function(names, finalClasses){

  names.classes.clean = names
  for(k in 2:length(names)){
    row.classes = names[[k]]
    classbranches = matrix(, nrow = k - 1, ncol = length(row.classes))
    for (j in 1:length(row.classes)){
      classbranch = character()
      for(i in 2:k){classbranch[i - 1] = substr(row.classes[j], 1, i)}
      classbranches[,j] = classbranch}
    classbranchesTF = apply(classbranches, 2, function(x){x%in%finalClasses})
    if(k == 2){
      IndexChangedClasses = which(classbranchesTF)
    } else {IndexChangedClasses = which(classbranchesTF, arr.ind = TRUE)[,2]}
    names.classes.clean[[k]][IndexChangedClasses] = classbranches[which(classbranchesTF, arr.ind = TRUE)]
  }
 return(names.classes.clean)
}
