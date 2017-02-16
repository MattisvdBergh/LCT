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
print.LCT <- function(res, ...){
  cat("tree object","\n\n")
  for(idxLevels in 1:(length(res$Names.clean) - 1)){
    cat("level",idxLevels,": ", res$Names.clean[[idxLevels]], "\n")
  }
  cat("Class Sizes final level:\n")
  print(res$ClassppGlobal[res$finalClasses])
  qgraph(computeEdges(res), layout = layout_as_tree, labels = rownames(res$ClassProportions))
}

computeEdges = function(res){
  library(qgraph)
  library(igraph)

  cleanNames = res$Names.clean
  Splits = res$Splits

  counter1 = counter2 = 1
  E = Enames = matrix(,0,2)
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

          EnamestoBe = cbind(level1Temp[idxClasses], level2Temp[splitClassesLogical])
          EtoBe = cbind(counter1, counter2)

          counter1 = counter1 + 1
          counter2 = max(EtoBe[,2])

          E = rbind(E,EtoBe)
          Enames = rbind(Enames, EnamestoBe)

        } else {counter1 = counter1 + 1}
      }}}
  return(E)
}

computeGlobalCpp = function(res){

  cpp = res$ClassProportions
  cpp = t(apply(cpp, 1, as.numeric))

  sizeAllSplits = unlist(res$Splits)[unlist(res$Splits)>1]
  allSplitClasses = unlist(sapply(1:length(sizeAllSplits),
                                  function(i){paste0(res$Splitpoints[i], 1:sizeAllSplits[i])}
  ))

  cppToBe = cpp
  splitsCpp = rownames(cpp)
  names(cppToBe) = paste0(0, 1:length(cppToBe))

  for(row in 2:nrow(cpp)){
    ncharSplit = nchar(splitsCpp[row])

    oldRow = substr(splitsCpp[row], 1, ncharSplit - 1)
    newCol = as.numeric(substr(splitsCpp[row], ncharSplit, ncharSplit))
    cppToBe[row,] = cppToBe[oldRow,newCol] * cppToBe[row,]
  }
  cppGtemp = as.numeric(t(cppToBe))
  tempNames = expand.grid(splitsCpp, 1:ncol(cpp))
  names(cppGtemp) = apply(
    tempNames[order(tempNames[,1]),],
    1, paste, collapse="")

  cppGtoReturn = cppGtemp[allSplitClasses]
  return(cppGtoReturn)
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
