plot.relImpr = function(res){
  nClass = length(res$IC[,"LL"][[1]])
  plot(1:(nClass - 1),
       c(1, relImpr(res)),
       type = "b", ylab = "", xlab = "Number of Classes",
       las = 1, bty ="n", axes = FALSE, ylim = c(0,1))
  axis(1, at = 1:(nClass - 1),
       labels = paste0(1:(nClass - 1),"-", 2:nClass))
  axis(2, las = 2)
}
