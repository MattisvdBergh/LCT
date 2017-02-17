computeEdgesSplits = function(res){
  namesSplitPoints = list()
  E = Enames = matrix(,0,2)
  for(idxLevels in 1:(length(res$Splits) - 1)){
    namesSplitPoints[[idxLevels]] = res$Names.clean[[idxLevels]][res$Splits[[idxLevels]] > 1]
  }

  counterSplitPoints = list(1)
  length(counterSplitPoints) = length(namesSplitPoints)
  counter = 2
  for(idxLevels in 2:length(namesSplitPoints)){
    for(idxClasses in 1:length(namesSplitPoints[[idxLevels]])){
      counterSplitPoints[[idxLevels]][idxClasses] = counter
      counter = counter + 1
    }}


  E = matrix( ,0, 2)
  counter1 = 2
  for(idxLevels in 2:length(namesSplitPoints)){
    for(idxClasses in 1:length(namesSplitPoints[[idxLevels]])){
      cc = namesSplitPoints[[idxLevels]][idxClasses]
      whichCount = substr(cc, 1, nchar(cc) - 1) == namesSplitPoints[[idxLevels - 1]]
      namesSplitPoints
      EtoBe = matrix(c(
        counterSplitPoints[[idxLevels - 1]][whichCount],
        counterSplitPoints[[idxLevels]][idxClasses]),
        ncol = 2)
      E = rbind(E, EtoBe)
    }}
  return(E)
}
