# rm(list=ls())
# LT = TRUE
# folder = "D:/"
# folder2 = "D:\\"
#
# if(LT == TRUE){
#   folder = "C:/Users/Mattis/"
#   folder2 = "C:\\Users\\Mattis\\"
# }
#
# LG = "D:/Application_Data/LatentGOLD5.1/lg51.exe"
# if(LT == TRUE){LG = "C:/Users/Mattis/Documents/LatentGOLD5.1/lg51.exe"}
#
# MainFolder = "SURFdrive/PhD Mattis/ProjectArtikel2/EVS/"
# MainFolder2 = "SURFdrive\\PhD Mattis\\ProjectArtikel2\\EVS\\"
#
# load(paste0(folder, MainFolder, "Results3split_10_5.Rdata"))
# Plot.PPlct(Results3, uptolevel = 3)

# Results = Results.G2
# # Results = test3
# uptolevel = 4
# filetype = "pdf"
# filename = "test"
# width = 7
# height = 7
# cex.axis = 1
# vsize = ""
# 
# Plot.PPlct(test, timeSequence = 0:22)
# Plot.PPlct(test3, timeSequence = 0:22, uptolevel = 3)

Plot.PPlct = function(Results,
                      uptolevel = 4,
                      timeSequence = 1:4,
                      filetype = "pdf",
                      filename = "test",
                      width = 7,
                      height = 7,
                      cex.axis = 1,
                      vsize = 10){
  # uptolevel = 4; filetype = "pdf"; filename = "test"; width = 7; height = 7

  library(igraph)
  library(qgraph)

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

  nNodes = sum(unlist(Results$Splits)>1)
  nc[!unlist(Results$Splits)>1] = NA

  allClasses = cbind(1:length(labels), labels, unlist(splits))
  importantNodes = allClasses[!nchar(allClasses[,2]) == max(nchar(allClasses[,2])),]

  nc2 = nchar(importantNodes[,2])
  labels2 = importantNodes[,2]

  E = matrix(,0,2)
  for(i in 1:nrow(importantNodes)){
  for (j in which(nc2 == nc2[i]+1)){
    if (grepl(paste0("^",labels2[i],"\\d$"),   labels2[j])){
      E = rbind(E,c(i,j))
    }}}

  # E = matrix(,0,2)
  # counter = 2
  # for(i in 1:nrow(foo)){
  #   if(foo[i,3] > 1){
  #     for(j in counter:(counter + foo[i,3])){
  #       E = rbind(E, c(i, j))
  #       counter = counter + 1
  #     }}
  #   }
  #
  # E = matrix(,0,2)
  # for (i in 1:nNodes){
  #   for (j in which(nc == nc[i]+1)){
  #     if (grepl(paste0("^",labels[i],"\\d$"),   labels[j])){
  #       E = rbind(E,c(i,j))
  #     }}}
  #
  # Etest = matrix(c(1, 2,
  #                  1, 3,
  #                  3, 4),
  #                3, 2, T)
  #
  # gTest = graph.edgelist(Etest)
  # Ltest = layout.reingold.tilford(gTest)
  # Ltest[4,1] = Ltest[4,1] - 0.25

  g = graph.edgelist(E)
  L = layout.reingold.tilford(g)

  SplitSizes = unlist(lapply(Results$Names[1:(uptolevel+1)], function(x){table(substr(x, 1, nchar(x)-1))}))[-1]
  SplitSizes = SplitSizes[1:nrow(L)]

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

  RPT = Results$ParametersTimepoints[,!apply(Results$ParametersTimepoints, 2, anyNA)]

  EV = apply(RPT, 2, function(x){
    exp(x[1] + x[2] * timeSequence + x[3] * timeSequence^2)/
      (1+exp(x[1] + x[2] * timeSequence + x[3] * timeSequence^2))
    })

  L.EVforPlot = list()
  length(L.EVforPlot) = length(Set.Splits)
  for(i in 1:length(Set.Splits)){
    if(length(Set.Splits[[i]])>1){
      EVforPlot = EV[,Set.Splits[[i]]]
      # L.EVforPlot = c(L.EVforPlot, list(EVforPlot))
      L.EVforPlot[[i]] = EVforPlot
    }}

  L.EVforPlot.temp = L.EVforPlot
  assign("L.EVforPlot", L.EVforPlot, envir = globalenv())

  rangeYlim = c(0,1)
  nPP = length(L.EVforPlot)
  list.profile.plots = list()

  EVlist = L.EVforPlot[!sapply(L.EVforPlot, is.null)]

  helpPaste = function(x){
    return(paste("c(",paste(x, collapse = ","), ")"))
  }

  countClass = level = 1
  listPP = list.profile.plots = list()
  for(i in 1:length(EVlist)){
    profile.plots.temp = parse(text = sprintf({
      "plot(timeSequence, EVlist[[%s]][, 1], ylim = c(0,1),
       xaxt=\"n\", xlab=\"\", las=2, ylab=\"\", axes=FALSE, type = \"b\", las=1)
       for(cols in 2:ncol(EVlist[[%s]])){
        lines(timeSequence, EVlist[[%s]][, cols], pch = cols, col = cols, type = \"b\")
       }
       axis(1)
       axis(2)"},
      i,
      i,
      i))
    listPP[[i]] = profile.plots.temp
  }

  finalPP = list()
  finalPP[which(!sapply(L.EVforPlot, is.null))] = listPP
  finalPP[which(sapply(L.EVforPlot, is.null))] = parse(text = "frame()")

  # relative = 15 + 5 / (1 + exp(-2*(L[, 2] - mean(L[, 2]))))
  # base = 1#76 / sum(relative)
  # vsize = relative * base

  # Eweight = rep(0, length(L.EVforPlot))
  # Eweight[which(!sapply(L.EVforPlot, is.null))] = 1
  # E = cbind(E, Eweight[-1])
  timeSequence <<- timeSequence
  EVlist <<- EVlist
  
  noSplits = sapply(finalPP, function(x){ifelse(class(x)=="call", TRUE, FALSE)})
  Lnew = L[!noSplits,]
  Enew = E[!E[,2]%in%which(noSplits),]
  finalPPnew = finalPP[!noSplits]
  
  Enew[3,2]=6
  
  counter = 2
  for(Row in 1:nrow(Enew)){
    rowTemp = Enew[Row,]
    if(diff(rowTemp) == 1){
      counter = Enew[Row,2]
    } else{ 
      counter = counter + 1
      Enew[Row, 2] = counter
  }}

  qgraph(Enew, layout = Lnew,
         subplots = finalPPnew,
         borders = FALSE, labels = FALSE, xpd = TRUE,
         vsize = vsize, mar = c(7,4,5,3))
  
  
  
  # qgraph(E, layout = L,
  #        subplots = finalPP,
  #        borders = FALSE, labels = FALSE, xpd = TRUE,
  #        vsize = vsize, mar = c(7,4,5,3))


}

