#' @title Plot the profile plots of a LCT in one image
#'
#' @param resTree A Latent Class Tree object.
#'
#' @description 1 pdf file with a growing tree.
#'
#' @export
plot1TreeExpl = function(resTree,
                         explTree,
                         rYlim = c(0, 1),
                         mar = c(1,1,1,1)){

  nLevels = max(nchar(names(explTree$splitInfo)))
  nLevels = length(resTree$treeSetup$cleanNames) - 1
  Splits = resTree$treeSetup$Splits
  splits = unlist(Splits)
  sizeSplits = splits[splits>1]
  sumSplits = sum(sizeSplits)

  layoutMatrix = matrix(, nrow = nLevels - 1, ncol = sumSplits)
  layoutMatrix[1, sumSplits/2] = 1

  counter = 2
  for(idxLev in 2:(nLevels - 1)){
    nClassRow = sum(Splits[[idxLev - 1]])
    idxParClass = which(!is.na(layoutMatrix[idxLev - 1,]))
    idxColTemp = round(seq(1, ncol(layoutMatrix), length.out = nClassRow+2))[-c(1, nClassRow + 2)]
    idxCol = idxColTemp[Splits[[idxLev]]>1]
    layoutMatrix[idxLev, idxCol] = counter: (counter + sum(Splits[[idxLev]]>1) - 1)
    counter = counter + sum(Splits[[idxLev]]>1)
  }

  idxEmpCol = apply(layoutMatrix, 2, function(x){all(is.na(x))})
  layoutMatrix = layoutMatrix[,!idxEmpCol]

  emptyValue = max(layoutMatrix, na.rm = TRUE) + 1
  lay = layoutMatrix
  lay[is.na(layoutMatrix)] = emptyValue

  layout(lay)
  par(mar = mar)
  plotLCTexpl(resTree = resTree,
              explTree = explTree,
                rYlim = rYlim)
}

