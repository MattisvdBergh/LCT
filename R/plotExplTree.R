# plot.explLCT = function(x, ...){
#
#
#
# x = list(splitInfo = explTree.SC2[-length(explTree.SC2)],
#          varInfo = c(explTree.SC2[length(explTree.SC2)],
#                         mLevelsLGdependent = data.frame(explTree.SC2$`0`$MeasurementLevels)))
#
#
#
#
#   for(idxSplits in 1:length(x$splitInfo)){
#
#     evSplit = x$splitInfo[[idxSplits]]$EV
#
#     for(idxVar in length(x$varInfo$mLevelsCovariates):1){
#
#       ev = evSplit[[idxVar]]
#
#       if(idxVar == 1){
#         par(new = TRUE)
#       }
#       nClass = length(x$splitInfo[[idxSplits]]$ClassProb)
#
#       if(x$varInfo$mLevelsLGdependent[idxVar] == "continuous"){
#         allEV = sapply(x$splitInfo, function(x) {x$EV[[idxVar]]})
#         rEV = range(allEV)
#         difRev = 0.1 * diff(rEV)
#         ylimcont = c(rEV[1] - difRev, rEV[2] + difRev)
#         plot(1:nClass, ev, xlim = c(0, nClass + 1), ylim = ylimcont,
#              axes= FALSE, type = "b", bty = "n")
#         axis(4)
#       }
#
#       if(x$varInfo$mLevelsLGdependent[idxVar] == "ordinal"){
#
#         bm = barplot(matrix(unlist(x$splitInfo[[idxSplits]]$EV[[idxVar]]),
#                       nrow = x$varInfo$mLevelsCovariates[idxVar]),
#                              beside = TRUE, col = rainbow(2))
#       }
#
#     }
#
#
# }
#
# explLCTcont = function(){
#
#   plot(x$splitInfo$`0`$EV[[1]], type = "b",
#        xlim = c(0, 2 + 1),
#        ylim = ylimcont,
#        axes = FALSE)
#   axis(4)
#
#   )
# }
#
#
#
#
#
# explLCTbar = function(x, idxS){
#   mp = barplot(1:2,matrix(x$`0`$EV[[idxVar]],                      ,
#                       nrow = x$mLevelsCovariates[idxVar]),
#                xlim = c(0, x$mLevelsCovariates[idxVar] + 1))
#
#
# }
#
#
#
# idxCont  = which(x$`0`$MeasurementLevels == "continuous")
#
#
#
#
#
