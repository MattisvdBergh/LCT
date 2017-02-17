computeEdges = function(res){

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
