#' @title Plot figures with the mean or proportion of the covariates for every split.
#'
#' @param explTree A Latent Class Tree to which the 3-step method is applied (with the function \code{exploreTree}).
#'
#' @description Several figures with the mean or proportion of the covariates for every split.
#'
#'@export
plotExplTreeCov = function(explTree,
                           ylab="Probability of belonging to a class",
                           xlab = "Value of the covariates",
                           ylim = c(0,1),
                           xlim = NULL,
                           lwd = 1,
                           ...){

  splitSize = sapply(1:length(explTree$splitInfo), function(i){length(explTree$splitInfo[[i]]$ClassProb)})
  mLevels = explTree$varInfo$mLevels
  sizeMlevels = explTree$varInfo$sizeMlevels
  parentClasses = names(explTree$splitInfo)

  for(idxSplits in 1:length(explTree$splitInfo)){

    evSplit = explTree$splitInfo[[idxSplits]]$EV
    childClasses = paste0(parentClasses[idxSplits], 1:splitSize[idxSplits])
    if(is.null(xlim)){
      xlim = range(as.numeric(unlist(lapply(evSplit, rownames))))
    }

    for(idxVar in 1:length(evSplit)){
      ev = evSplit[[idxVar]]

      for(idxClass in 1:ncol(ev)){
        if(idxVar != 1 | idxClass != 1){par(new = TRUE)}

        lty = idxClass
        if(lty > 2){lty = lty + 1}

        plot(as.numeric(rownames(ev)), ev[,idxClass], ylim = ylim,
             xlim = xlim, axes= FALSE, type = "l", bty = "n", ylab = "",
             xlab = "", lty = lty, col = idxVar, lwd = lwd)

        # if(idxVar != 2){
        # text(as.numeric(rownames(ev))[1], ev[1,idxClass],
        #      paste0(substr(parentClasses[idxSplits], 2, nchar(parentClasses[idxSplits])),
        #             idxClass),
        #      xpd = TRUE)
        # } else{
        #   text(as.numeric(rownames(ev))[nrow(ev)], ev[nrow(ev),idxClass],
        #        paste0(substr(parentClasses[idxSplits], 2, nchar(parentClasses[idxSplits])),
        #               idxClass),
        #        xpd = TRUE)
        # }


      }
    }
    axis(1)
    axis(2, las = 1)
    mtext(xlab, side = 1, line = 3, cex = par()$cex)
    mtext(ylab, side = 2, line = 3, cex = par()$cex)
  }
}
