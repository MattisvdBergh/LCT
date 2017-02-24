expectedValuesGrowthOrd = function(parmsEV,
                                   classes,
                                   vecInd,
                                   nIndependent,
                                   levelsDependent){

  timeSequence = sort(unique(vecInd[,1]))
  basis = timeSequence
  for(idxInd in 2:nIndependent){
    basis = cbind(basis, timeSequence^idxInd)
  }

  toReturn = array(, c(length(classes), levelsDependent, length(timeSequence)),
                   dimnames = list(classes,
                                   (1:levelsDependent),
                                   timeSequence))

  for(colParmsEV in classes){

    etaTemp = matrix(, nrow = nrow(basis), ncol = levelsDependent)
    tempParmsEV = parmsEV[!is.na(parmsEV[,colParmsEV]),colParmsEV]

    for(idxDep in 1:levelsDependent){
      intercept = tempParmsEV[idxDep]
      etas = intercept +
        colSums(tempParmsEV[(levelsDependent+1):length(tempParmsEV)] * t(basis)*(idxDep - 1))
      etaTemp[,idxDep] = etas
    }
    probabilities = apply(etaTemp, 1, function(x){
      probs = numeric(length(x))
      for(i in 1:length(x)){probs[i] = exp(x[i])/sum(exp(x))}
      return(probs)})
    toReturn[colParmsEV,,] = probabilities
  }
  return(toReturn)
}

meanProb = function(prob){
  return(apply(prob, 3, function(x){
    apply(x, 1, function(i){
      sum(i*(0:(length(i) - 1)))}
    )}))
}
