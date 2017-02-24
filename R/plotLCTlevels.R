#' @title Plot the class means or conditinal probabilities of every level of a LCT
#'
#' @param resTree A Latent Class Tree object
#' @param rYlim the ylim of a plot of a continuous variable. By default for categorical variables c(0, 1) is used.
#' @param doPlot Logical, does the object need to be printed to a file as indicated by \code{dev}
#' @param nameTree if \code{doPlot} is TRUE (default), a folder called 'Levels'\code{nameTree} will be created.
#' In this folder of every level of the tree a plot is created, which is called \code{nameTree},
#' followed by the name of the split at hand.
#' @param dev device to create a plot, either "png", or "pdf"
#'
#' @description Produce the class means or conditinal probabilities of every split of a LCT
#'
#' @export
plotLCTLevels = function(resTree,
                         rYlim = c(0,1),
                         doPlot = FALSE,
                         nameTree = "",
                         dev = "png"){

  if(doPlot == TRUE){
  mainDir = getwd();  subDir = paste0("Levels", nameTree)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  }

  evTree = resTree$splitInfo$EV
  namesTree = resTree$treeSetup$cleanNames
  sizeMlevels = resTree$splitInfo$sizeMlevels
  splits = unlist(resTree$treeSetup$Splits)
  spSize = splits[splits>1]
  namesPlots = list()

  if(all(sizeMlevels == 2)){
    evTree = evTree[seq(2, nrow(evTree), 2), ]
  }

  for(idxLevels in 2:length(namesTree)){

    namesLevel = namesTree[[idxLevels]]
    evLevel = evTree[,namesLevel]
    parentClasses = substr(namesLevel, 1,  nchar(namesLevel) - 1)
    pSplits = resTree$treeSetup$Splits[[idxLevels - 1]]

    xaxes = list()
    counter = 1
    for(idxSplits in seq_along(pSplits)){
      xaxes[[idxSplits]] = c(counter:(pSplits[idxSplits] + counter - 1))
      counter = max(unlist(xaxes)) + 1
    }


    if(doPlot == TRUE){
      if(dev == "png"){
    png(paste0(nameTree, "_", idxLevels, ".png"))
      }
      if(dev == "pdf"){
        png(paste0(nameTree, "_", idxLevels, ".png"))
      }
    }


    plot(unlist(xaxes), rep(-10, sum(lengths(xaxes))),
         ylim = rYlim, xaxt = "n", xlab = "", las = 2,
         ylab = "", axes = FALSE, type = "b")

    for(idxSplits in seq_along(xaxes)){
      for(rows in 1:nrow(evLevel)){
        lines(xaxes[[idxSplits]],
              evLevel[rows, xaxes[[idxSplits]]],
              pch = rows, col = rows, type = "b")}
    }

    axis(1, at = unlist(xaxes), labels = namesLevel,
         col = "white", cex.axis = 1.5)
    axis(2, las = 2, cex.axis = 1.5)

    if(doPlot == TRUE){
      dev.off()
    }
  }

  itemNames = rownames(resTree$splitInfo$EV)
  nItems = length(itemNames)

  if(doPlot == TRUE){
  png("Legend.png")
  plot.new()
  legend("center", itemNames, pch = 1:nItems, col = 1:nItems, bty = "n", cex = 2, xpd = TRUE)
  dev.off()
  setwd(mainDir)
  }
}
