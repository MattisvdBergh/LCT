# allSplitPosteriors =
#   res$Posteriors[,
#                  get(
#                    paste0(
#                    "res$Posteriors$W_", names(nClassTotal),""
#                  ]


getPosteriors = function(resTree){
  uniqueLevels = 2:(length(resTree$Names.clean) - 1)
  startPosteriors = apply(resTree$Posteriors, 2, function(x){as.numeric(as.character(x))})
  L.postFinal = list()

  for(cl in uniqueLevels){
    namesLevel = resTree$Names.clean[[cl]]
    branches = substr(namesLevel, 1, nchar(namesLevel) - 1)
    uBranches = unique(branches)
    postTemp = startPosteriors[, paste0("Post_", namesLevel)]

    if(length(uBranches) == 1){
      postFinal = postTemp
    } else {

      maxNcharNamesLevel = max(nchar(namesLevel))
      counter = 1

      for(idxClasses in seq_along(namesLevel)){
        if(nchar(namesLevel[idxClasses]) < maxNcharNamesLevel){
          noSplit = namesLevel[nchar(namesLevel) != max(nchar(namesLevel))][counter]
          counter = counter + 1
          postFinal = cbind(postFinal, postFinal[, paste0("Post_", noSplit)])
          colnames(postFinal)[ncol(postFinal)] = paste0("Post_", noSplit, "_", cl)
        } else {

          postTemp2 = postTemp[, paste0("Post_", namesLevel[idxClasses])]

          postNew = postTemp2 * postFinal[,paste0("Post_",
                                                  substr(namesLevel[idxClasses], 1, nchar(namesLevel[idxClasses]) - 1))]
          postFinal = cbind(postFinal, postNew)
          colnames(postFinal)[ncol(postFinal)] = paste0("Post_", namesLevel[idxClasses])
        }}}}
  return(postFinal)
}

