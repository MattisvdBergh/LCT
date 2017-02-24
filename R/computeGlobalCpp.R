computeGlobalCpp = function(ClassProportions, Splits, Splitpoints){

  cpp = ClassProportions
  cpp = t(apply(cpp, 1, as.numeric))

  sizeAllSplits = unlist(Splits)[unlist(Splits)>1]
  allSplitClasses = unlist(sapply(1:length(sizeAllSplits),
                                  function(i){paste0(Splitpoints[i], 1:sizeAllSplits[i])}
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
