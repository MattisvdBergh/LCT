#' @title Plot figures with the mean or proportion of the covariates for every split.
#'
#' @param explTree A Latent Class Tree to which the 3-step method is applied (with the function \code{exploreTree}).
#'
#' @description Several figures with the mean or proportion of the covariates for every split.
#'
#'@export
plotLCTexpl = function(explTree,
                       ylab="",
                       ylab2 = "",
                       beside = FALSE, ...){

  splitSize = sapply(1:length(explTree$splitInfo), function(i){length(explTree$splitInfo[[i]]$ClassProb)})
  mLevels = explTree$varInfo$mLevels
  sizeMlevels = explTree$varInfo$sizeMlevels
  parentClasses = names(explTree$splitInfo)

  for(idxSplits in 1:length(explTree$splitInfo)){

    evSplit = explTree$splitInfo[[idxSplits]]$EV
    childClasses = paste0(parentClasses[idxSplits], 1:splitSize[idxSplits])

    for(idxVar in 1:length(evSplit)){
      if(idxVar != 1){par(new = TRUE)}

      ev = evSplit[[idxVar]]

      if(mLevels[idxVar] == "ordinal"){
        bm = barplot(matrix(
          unlist(ev),nrow = sizeMlevels[idxVar]),
          beside = beside, col = rainbow(sizeMlevels[idxVar]), ylim = c(0,1),
          names.arg = substr(childClasses, 2, nchar(childClasses)), las = 1, ylab = ylab)
      }

      if(mLevels[idxVar] == "continuous"){
        allEV = sapply(explTree$splitInfo, function(x) {x$EV[[idxVar]]})
        rEV = range(allEV)
        difRev = 0.1 * diff(rEV)
        ylimcont = c(rEV[1] - difRev, rEV[2] + difRev)
        plot(1:splitSize[idxSplits], ev,
             xlim = c(0, splitSize[idxSplits] + 1), ylim = ylimcont,
             axes= FALSE, type = "b", bty = "n", xlab = "Classes", ylab = "")
        axis(4, las = 2)
        mtext(ylab2, side = 4, line = 3, cex = 2)
      }
    }
  }
}


