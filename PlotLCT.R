plot.lct = function(Results, 
                      uptolevel = NULL, 
                      vsize = ""){
  library(igraph)
  library(qgraph)
  
  if(!is.null("uptolevel")){uptolevel = length(Results$Names) - 1}
  
  data = Results$Names[1:uptolevel]
  splits = Results$Splits[1:uptolevel]
  nodesPerLevel = sapply(data,length)
  labels = unlist(data)
  nNodes = length(labels)
  nc = nchar(labels)
  if(vsize == ""){vsize = 15*exp(-nNodes/80)+5}
  nPP = length(unlist(Results$Splits[1:(uptolevel)]))
  list.profile.plots = list()
  level = 1
  
  E = matrix(,0,2)
  for (i in 1:nNodes){
    for (j in which(nc == nc[i]+1)){
      if (grepl(paste0("^",labels[i],"\\d$"),   labels[j])){
        E = rbind(E,c(i,j))
      }}}
  
  g = graph.edgelist(E)
  L = layout.reingold.tilford(g)
  
  SplitSizes = unlist(lapply(Results$Names[1:(uptolevel+1)], function(x){table(substr(x, 1, nchar(x)-1))}))[-1]
  
  Set.Splits = list()
  counter = 1
  for(i in 1:length(SplitSizes)){
    # for(i in 1:10){
    if(SplitSizes[i] > 1){
      Set.Splits[[counter]] = paste0(names(SplitSizes)[i], 1:SplitSizes[i])
      counter = counter + 1
    } else {
      Set.Splits[[counter]] = paste0(names(SplitSizes)[i], 1:SplitSizes[i])
      counter = counter + 1
    }}
  
  EV = Results$EV
  L.EVforPlot = list()
  length(L.EVforPlot) = length(Set.Splits)
  for(i in 1:length(Set.Splits)){
    if(length(Set.Splits[[i]])>1){
      EVforPlot = EV[,Set.Splits[[i]]]
      L.EVforPlot[[i]] = EVforPlot
    }}
  
  L.EVforPlot.temp = L.EVforPlot
  assign("L.EVforPlot", L.EVforPlot, envir = globalenv())
  
  rangeYlim = c(0,1)
  nPP = length(L.EVforPlot)
  list.profile.plots = list()
  
  l.vsize = numeric(length(Set.Splits))
  l.vsize[lengths(Set.Splits)>1] = vsize
  l.vsize[lengths(Set.Splits)==1] = vsize/2
  
  level=1
  list.profile.plots = list()
  for(i in 1:length(L.EVforPlot)){
    if(length(Set.Splits[[i]])>1){
    profile.plots.temp = parse(text = sprintf({
      "plot(%s, EV[1, %s], ylim = c(0,1),
      xaxt=\"n\", xlab=\"\", las=2, ylab=\"\", axes=FALSE, type = \"b\")
      for(rows in 1:nrow(EV)){
        lines(%s, EV[rows, %s], pch = rows, col = rows, type = \"b\")
        }
        axis(1, at=%s, labels = %s,col = \"white\", cex.axis = 0.5)
      %s"},
      paste("c(", paste(1:length(Set.Splits[[i]]), collapse = ","), ")"),
      paste("c(", paste(which(colnames(EV) %in% Set.Splits[[i]]), collapse = ","), ")"),
      paste("c(", paste(1:length(Set.Splits[[i]]), collapse = ","), ")"),
      paste("c(", paste(which(colnames(EV) %in% Set.Splits[[i]]), collapse = ","), ")"),
      paste("c(", paste(1:length(Set.Splits[[i]]), collapse = ","), ")"),
      paste("c(", paste(Set.Splits[[i]], collapse = ","), ")"),
      if(level < max(nchar(Set.Splits[[i]]))){
        level = max(nchar(Set.Splits[[i]]))
        paste("axis(2, cex.axis = 0.5, las = 2)")
      } else {
        paste("")
        }))
    list.profile.plots[[i]] = profile.plots.temp
    } else {
      WhichSplit = Set.Splits[[i]]
      while(!any(colnames(EV) %in% WhichSplit)){
      WhichSplit = substr(WhichSplit, 1, nchar(WhichSplit) - 1)}
      profile.plots.temp = parse(text = sprintf({
        "plot(EV[1, %s], ylim = c(0,1),
      xaxt=\"n\", xlab=\"\", las=2, ylab=\"\", axes=FALSE, type = \"b\")
      for(rows in 1:nrow(EV)){
        points(EV[rows, %s], pch = rows, col = rows, type = \"b\")
        }
        axis(1, at=1, labels = %s,col = \"white\", cex.axis = 0.5)
        %s"},
        paste("c(", paste(which(colnames(EV)==WhichSplit), collapse = ","), ")"),
        paste("c(", paste(which(colnames(EV)==WhichSplit), collapse = ","), ")"),
        paste("c(", paste(Set.Splits[[i]], collapse = ","), ")"),
        if(level < max(nchar(Set.Splits[[i]]))){
          level = max(nchar(Set.Splits[[i]]))
          paste("axis(2, cex.axis = 0.5, las = 2)")
        } else {
          paste("")
      }))
      list.profile.plots[[i]] = profile.plots.temp
      }}

  qgraph(E, layout = L, 
         subplots = list.profile.plots,
         borders = FALSE, labels = FALSE, xpd = TRUE,
         vsize = vsize)
}  

