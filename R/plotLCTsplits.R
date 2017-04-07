#' @title Plot the class means or conditinal probabilities of every split of a LCT
#'
#' @param resTree A Latent Class Tree object
#' @param rYlim the ylim of a plot of a continuous variable. By default for categorical variables c(0, 1) is used.
#' @param doPlot Logical, does the object need to be printed to a file as indicated by \code{dev}
#' @param nameTree if \code{doPlot} is TRUE (default), a folder called Levels'nameTree' will be created.
#' In this folder of every level of the tree a plot is created, which is called 'nameTree',
#' followed by the name of the split at hand
#' @param dev device to create a plot, either "png", or "pdf"
#'
#' @description Produce a figure for every split of the tree
#'
#'@export
plotLCTSplits = function(resTree,
                         rYlim = c(0, 1),
                         doPlot = FALSE,
                         nameTree = "",
                         dev = "png",
                         doEq = FALSE,
                         idxSP = NULL,
                         ylab = "Mean score",
                         xlab = "Classes",
                         lwd = 1){

  if(doPlot == TRUE){
    mainDir = getwd();  subDir = paste0("Splits", nameTree)
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    setwd(file.path(mainDir, subDir))
  }

  if(any("meanEV" %in% names(resTree$splitInfo))){
    evTreeAll = resTree$splitInfo$meanEV
    evTree = t(do.call(rbind, evTreeAll))
  } else {
  evTreeAll = resTree$splitInfo$EV
  if(all(resTree$splitInfo$sizeMlevels == 2)){
    evTree = evTreeAll[seq(2, nrow(evTreeAll), 2),]
  } else {
    evTree = evTreeAll
  }}

  sp = allSP = resTree$treeSetup$parentClasses
  if(!is.null(idxSP)){
    sp = allSP[idxSP]
  }
  splits = unlist(resTree$treeSetup$Splits)
  spSize = splits[splits>1]
  namesPlots = exprPlots = list()

  for(idxSplits in 1:length(sp)){

    spClasses = apply(expand.grid(sp[idxSplits], 1:spSize[idxSplits]),
                      1,
                      paste0, collapse="")

    evSplit = evTree[,spClasses]

    namesPlots[[idxSplits]] = paste0(getwd(),"/",
                                     paste0(nameTree, "Split", sp[idxSplits], ".png"))

    if(doPlot == TRUE){
      if(dev == "png"){
        png(paste0(nameTree, "_", sp[idxSplits], ".png"))
      }
      if(dev == "pdf"){
        png(paste0(nameTree, "_", sp[idxSplits], ".png"))
      }
    }

    if(any("meanEV" %in% names(resTree$splitInfo))){
      plot(evSplit[,1], ylim = rYlim,
           xlab = xlab, ylab = ylab, lwd = lwd,
           main = sp[idxSplits], type = "l", bty = "n", axes = FALSE)
      for(idxClass in 2:ncol(evSplit)){
        lty = idxClass
        if(lty > 2){lty = lty + 1}
        lines(evSplit[,idxClass], pch = idxClass, col = idxClass, type = "l", lty = lty, lwd = lwd)
      }
      axis(1)
      axis(2, las = 2)
    } else{

    plot(evSplit[1,], ylim = rYlim,
         xlab = "Classes", ylab = ylab,
         type = "b", bty = "n", axes = FALSE)
    for(idxVar in 2:nrow(evSplit)){
      lines(evSplit[idxVar,], pch = idxVar, col = idxVar, type = "b")
    }
    axis(1, at = 1:spSize[idxSplits], substr(spClasses, 2, nchar(spClasses)))
    axis(2, las = 2)

    }

    if(doPlot == TRUE){
      dev.off()
    }

    exprPlots[[idxSplits]] = parse(text = paste("evSplit = matrix(c(", paste(evSplit, collapse = ","), "),
                                  nrow =", as.numeric(paste(nrow(evSplit))), ")
                                   spSize = c(", paste(spSize, collapse = ","), ")
                                   spClasses = c(", paste(spClasses, collapse = ","), ")
                                   rownames(evSplit) = c('", paste(rownames(evSplit),collapse = "','"),"')
                                   plot(evSplit[1,], ylim = c(",
                                                paste(rYlim, collapse = ","),
                                                "),
                                   xlab = '', ylab = '', type = 'b', bty = 'n', axes = FALSE, xlim = c(0, length(spClasses) + 1))
                                   for(idxVar in 2:nrow(evSplit)){
                                   lines(evSplit[idxVar,], pch = idxVar, col = idxVar, type = 'b')
                                   }
                                   axis(1, at = 1:length(spClasses), rep('', length(spClasses)), cex.axis = 0.5)
                                   axis(2, las = 2, cex.axis = 0.5)"))
  }

  itemNames = rownames(evSplit)
  nItems = length(itemNames)

  if(doPlot == TRUE){
    png("Legend.png")
    plot.new()
    legend("center", itemNames, pch = 1:nItems, col = 1:nItems, bty = "n", cex = 2, xpd = TRUE)
    dev.off()
    setwd(mainDir)
  }


  if(doPlot == TRUE & doEq == FALSE){
    return(namesPlots)
  }

  if (doEq == TRUE & doPlot == FALSE){
  return(exprPlots)
  }

  if(all(c(doPlot, doEq) == TRUE)){
    return(list(dirPlots = namesPlots,
                exprPlots = exprPlots))
  }
}
