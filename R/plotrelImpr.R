#' Title
#'
#' @param resTree A Latent Class Tree object
#'
#' @return null
#' @export
plotRelImpr = function(resTree){
  nClass = length(resTree$splitInfo$IC[,"LL"][[1]])
  plot(1:(nClass - 1),
       c(1, relImpr(resTree)),
       type = "b", ylab = "", xlab = "Number of Classes",
       las = 1, bty ="n", axes = FALSE, ylim = c(0,1))
  axis(1, at = 1:(nClass - 1),
       labels = paste0(1:(nClass - 1),"-", 2:nClass))
  axis(2, las = 2)
}
