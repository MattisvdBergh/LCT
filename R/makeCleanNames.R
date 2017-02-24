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
