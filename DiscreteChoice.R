# Ctrl + f op @#$ voor wijzigingen.


LCT = function(dataDir, LG, LGS, alternativesDir, choicesDir,
               decreasing = TRUE, maxClassSplit1 = 2, maxClassSplit2 = 2,
               stopCriterium = "BIC", resultsName = "", diffStopCriterium = 0,
               minSampleSize = 10){
  
  library(abind)
  # Read the original data
  mydata = read.table(dataDir, header = TRUE)
  ResponseOptions = length(unique(mydata$choice))
  Nitems = length(unique(mydata$set))
  ObservationsPerRespondent = mean(unname(table(mydata$id)))
  weight = mydata$weight

  # The tree is summarized in two lists. 
  # One with the sizes of the splits, and one with the names.  
  namesClasses = splitsClasses = list()
  Hclass = 1        # Level of the hierarchy
  
  mainDir = getwd();  subDir = paste0("Results", resultsName)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  
  
  # Create results folder and put LCT syntax there
  syntax = readLines(LGS)
  makeSyntaxDiCh(dataDir, Hclass, syntax, maxClassSplit1, alternativesDir = alternativesDir, choicesDir = choicesDir)
  
  # Perform LC for 1- and 2-classes ####
  shell(paste(LG, "LCT0.lgs", "/b"))
  
  # Results ####
  
  # Write results to LCT file
  resultsTemp = readLines("LCT0.lst")
  
  helpFun = function(x, greplIdx, idx = NULL) {
    x = x[grepl(greplIdx, x)]
    if (is.null(idx)) idx = seq_along(x)
    return(as.numeric(gsub(",", ".", sapply(strsplit(x, "\t")[idx], `[`, 2))))
  }
  
  Ntot = helpFun(resultsTemp, "Number of cases")[1]
  if(minSampleSize < 1){minSampleSize = Ntot * minSampleSize}
  
  out = getOutputDiCh(resultsTemp, Hclass, maxClassSplit1, Ntot = Ntot, stopCriterium = stopCriterium, decreasing = decreasing, ResponseOptions = ResponseOptions, Nitems = Nitems)
  
  ## Update the tree after the first split
  splitsClasses[[Hclass]] = unname(out$solution)
  namesClasses[[1]] = "0"
  namesClasses[[2]]  = paste0("0", 1:out$solution)
  Splitpoints = 0
  
  # The above procedure repeats, while not all splits at the lowest level are 1
  ## Start while loop ---------------------------------
  
  # If there are classes which can be split up, repeat
  while (!all(splitsClasses[[length(splitsClasses)]] == 1)) {
    # while (Hclass==1 |
    #   !sum(splitsClasses[[length(splitsClasses)]] == 2) == 1){
    
    # Update variables with the position in the tree:
    # sizePreviousSplits keeps track of the split size of every previous split.
    # nPreviousSplits keeps track of the number of previous splits.
    Hclass = Hclass + 1 
    sizePreviousSplits = 1:splitsClasses[[length(splitsClasses)]][1] # which classes of the previous split are to be tested
    nPreviousSplits = 1:length(splitsClasses[[length(splitsClasses)]]) # originating from which class in the previous split
    
    # namesClassesTemp will contain the new names of the split classes and will be later added to namesClasses
    # splitsClassesTemp will contain the split sizes of the new splits and will be later added to splitsClasses
    namesClassesTemp = splitsClassesTemp = numeric() # to be filled up with number of classes at each hierarchical level
    
    # At each hierarchical level, we have to run through the the number of splits performed at the previous level.
    for (wPreS in nPreviousSplits) {
      
      # The size of the previous split must updated
      if (wPreS != 1) {
        sizePreviousSplits = sizePreviousSplits + 
          splitsClasses[[length(splitsClasses)]][wPreS - 1]
      }
      
      # Correction if we switch number of classes
      length(sizePreviousSplits) = splitsClasses[[length(splitsClasses)]][wPreS]
      if (any(is.na(sizePreviousSplits))) { # zie anyNA @#$
        for (k in which(is.na(sizePreviousSplits))) {
          sizePreviousSplits[k] = sizePreviousSplits[k - 1] + 1
        }
      }
      
      # If no previous split is made, no next split has to be checked and we go to the next split
      if (length(namesClasses[[length(namesClasses)]][sizePreviousSplits]) > 1) {
        
        # For every new class of the next previous split
        preC = namesClasses[[length(namesClasses)]][sizePreviousSplits]
        for (iCC in seq_along(namesClasses[[length(namesClasses)]][sizePreviousSplits])) {
        
          ## Make new data ####
          data.temp = read.delim(paste0("H", Hclass-1,
                                        "c", substr(preC[iCC], 1, Hclass - 1), "_sol",
                                        splitsClasses[[length(splitsClasses)]][wPreS],
                                        ".txt"), dec = ",")
          
          idxWeight = which(colnames(data.temp) == "weight")
          
          weight.perClass = data.temp[, idxWeight:(ncol(data.temp) - 1)]
          # weight.perClass[weight.perClass[,1] == ".", ] = NA
          
          if(is.factor(weight.perClass[,1])){
            weight.perClass = apply(weight.perClass, 2, function(x){as.numeric(as.character(gsub(",", ".", x)))})
          }
          weight.perClass[is.na(weight.perClass[,1]), ] = 0
          
          if (is.factor(data.temp$weight)) {
            data.temp$weight = as.numeric(levels(data.temp$weight))[data.temp$weight]
          }
          
          order.Classpp = order(colSums(
            weight.perClass[,1] * weight.perClass[,-1], na.rm = TRUE
          ), decreasing = decreasing)
          NewWeight = weight.perClass[,
                                      order.Classpp[
                                        which(
                                          preC[iCC]==namesClasses[[
                                            length(namesClasses)]][
                                              sizePreviousSplits]
                                        )] + 1] *
            weight.perClass[,1]
          mydata = cbind(data.temp[,1:(idxWeight - 1)], weight = NewWeight)
          write.table(mydata, paste0("mydata", Hclass, preC[iCC], ".dat"),
                      sep = "\t", row.names = FALSE, quote = FALSE)
          
          ## Update LGS ####
          
          # New datafile
          syntax = readLines(LGS)
          makeSyntaxDiCh(dataDir = paste0(
            normalizePath(getwd()), "\\mydata", Hclass, preC[iCC], ".dat"),
            Hclass = Hclass,
            syntax = syntax,
            maxClassSplit1 = maxClassSplit2,
            alternativesDir = alternativesDir, 
            choicesDir = choicesDir,
            CC = preC[iCC])
          
          ## Perform LC for 1- and 2-classes ####
          shell(paste0(LG, " LCT", preC[iCC],".lgs /b"))
          
          ## Results ####
          
          ## Write results to LCT file
          resultsTemp = readLines(paste0("LCT", preC[iCC],".lst"))
          outTemp = getOutputDiCh(resultsTemp, Hclass, maxClassSplit2, CC = preC[iCC], Ntot = Ntot, stopCriterium = stopCriterium, itemNames = itemNames, decreasing = decreasing, ResponseOptions = ResponseOptions, Nitems = Nitems)
          
          # EV = abind(EV, EV.temp, along = 1)
          
          
          idxRbind = 1:2
          idxCbind = 3
          idxAbind = 4
          maxCol = max(ncol(out[[2]]), ncol(outTemp[[2]]))
          maxRow = max(nrow(out[[4]]), nrow(outTemp[[4]]))
          while(ncol(outTemp[[2]]) < maxCol){outTemp[[2]] = cbind(outTemp[[2]], NA)}
          while(ncol(outTemp[[3]]) < maxCol){outTemp[[3]] = cbind(outTemp[[3]], NA)}
          # while(maxRow > nrow(outTemp[[4]])){outTemp[[4]] = rbind(outTemp[[4]], NA)}
          
          temp1 = Map(function(x, y) rbind(x,y), out[idxRbind], outTemp[idxRbind])
          temp2 = Map(function(x, y) cbind(x,y), out[idxCbind], outTemp[idxCbind])
          temp3 = Map(function(x, y, z) abind(x,y, along = z), out[idxAbind], outTemp[idxAbind], 1)
          
          out = c(temp1, temp2, temp3)
          solution = outTemp$solution
          if(solution > 1){Splitpoints = c(Splitpoints, preC[iCC])}
          
          ## Update the names and split sizes if there were splits
          namesClassesTemp = c(namesClassesTemp, paste0(preC[iCC], 1:solution))
          splitsClassesTemp = c(splitsClassesTemp, solution)
        }} else { # Update the names and split sizes if there were no splits
          namesClassesTemp =
            c(namesClassesTemp,
              paste0(namesClasses[[length(namesClasses)]][sizePreviousSplits], 1))
          splitsClassesTemp = c(splitsClassesTemp, 1)
        }}
    # save temporary objects
    namesClasses[[length(namesClasses) + 1]] = namesClassesTemp
    splitsClasses[[length(splitsClasses) + 1]] = unname(splitsClassesTemp)
    print(namesClasses)
  }

  ### Global class proportions can be found in ClassppGlobal
  Classpp = out$Classpp[!apply(out$Classpp, 1, anyNA),]
  Classpp.matrix = t(Classpp)
  Classpp.numeric = as.numeric(Classpp.matrix)
  Classpp.numeric.names = character()
  for(l in 1:ncol(Classpp.matrix)){
    Classpp.numeric.names = c(Classpp.numeric.names,
                              paste0(colnames(Classpp.matrix)[l], 1:nrow(Classpp.matrix)))}  
  Missing = is.na(Classpp.numeric)
  Classpp.numeric = Classpp.numeric[!Missing]
  names(Classpp.numeric) = Classpp.numeric.names[!Missing]
  
  ClassppGlobal = Names.sub.Classes = numeric()
  for(i in 1:length(Classpp.numeric)){
    if(nchar(names(Classpp.numeric)[i])>2){
      Names.sub.Classes = numeric()
      for(j in 1:(nchar(names(Classpp.numeric)[i])-2)){
        Names.sub.Classes [j] = 
          substr(names(Classpp.numeric)[i], 1, nchar(names(Classpp.numeric)[i]) - j)
      }
      ClassppGlobal[i] = 
        prod(c(Classpp.numeric[i], Classpp.numeric[Names.sub.Classes]))
    } else{
      ClassppGlobal[i] = 
        Classpp.numeric[i]
    }
  }
  names(ClassppGlobal) = names(Classpp.numeric)
  
  ### Which classes are split can be found in names.split.classes
  names.split.classes = character()
  for(l in 1:length(Splitpoints)){
    names.split.classes = c(names.split.classes, 
                            paste0(Splitpoints[l],
                                   1:unlist(splitsClasses)[unlist(splitsClasses)!=1][l]))
  }
  
  ### Classes at the lowest level are put in final.classes
  names.split.classes.temp = names.split.classes
  for(i in 1:length(names.split.classes)){
    if(any(nchar(names.split.classes[i]) < nchar(names.split.classes) &
           grepl(names.split.classes[i], names.split.classes))){
      names.split.classes.temp[i] = NA
    }
  }  
  final.classes = names.split.classes.temp[!is.na(names.split.classes.temp)]
  
  
  #############
  
  
  names.classes.clean = namesClasses
  for(k in 2:length(namesClasses)){
    row.classes = namesClasses[[k]]
    classbranches = matrix(, nrow = k - 1, ncol = length(row.classes))
    for (j in 1:length(row.classes)){
      classbranch = character()
      for(i in 2:k){classbranch[i - 1] = substr(row.classes[j], 1, i)}
      classbranches[,j] = classbranch}
    classbranchesTF = apply(classbranches, 2, function(x){x%in%final.classes})
    if(k == 2){
      IndexChangedClasses = which(classbranchesTF)
    } else {IndexChangedClasses = which(classbranchesTF, arr.ind = TRUE)[,2]}
    names.classes.clean[[k]][IndexChangedClasses] = classbranches[which(classbranchesTF, arr.ind = TRUE)]
  }
  
  # End function ####
  results = list(namesClasses, out$IC, splitsClasses, out$EV, out$Classpp, out$Post,
                 Splitpoints, ClassppGlobal, names.split.classes, final.classes)
  names(results) = c("Names", "IC", "Splits", "EV", "Classpp", "Posteriors",
                     "Splitpoints", "ClassppGlobal", "Names.Split.Classes",
                     "Final.Classes")
  setwd(mainDir)
  return(results)
}

# LCT provides a list where:
## $Names contains the class names,
## $IC contains the LL and BIC of every split,
## $Splits contains the size of every split
## $EV contains the estimated conditional probabilities
## $Classpp contains the Class proportions
## $Posteriors is a large matrix with subsequently for every class:
### the weights (proportions) for every response pattern/participant,
### followed by the posterior probabilities of the split of that class
## $Rbic contains the relative BIC improvement
## $Entropy contains the Entropy of each split
## $CE contains the Classification Error of each split







makeSyntaxDiCh = function(dataDir, Hclass, syntax, maxClassSplit1, alternativesDir, choicesDir, CC = 0) {
  
  syntax[grep("infile", syntax)] = capture.output(cat(paste0("infile \'", dataDir, "'")))
  syntax[grep("alternatives", syntax)] = capture.output(cat(paste0("alternatives \'", alternativesDir, "'")))
  syntax[grep("choicesets", syntax)] = capture.output(cat(paste0("choicesets \'", choicesDir, "'")))
  k = 1:maxClassSplit1
  
  # sL is the first line (starting line) of the first model
  # lm is the number of lines of one model.
  # bm is the beginning where the next model should start
  # em is the end where the next model should end.
  sL = which(syntax=="model")[1]
  lm = which(syntax == "end model")[2] - which(syntax == "model")[2] + 1
  bm = which(syntax=="model")[1] + (k - 1) * lm
  em = (which(syntax=="model")[1] - 1) + k * lm
  
  newSyntax = character(em[length(em)])
  newSyntax[1:(sL - 1)] = syntax[1:(sL - 1)]
  newSyntax[sL:em[length(em)]] = rep(syntax[sL:em[1]],
                                    maxClassSplit1)
  tempCsv = strsplit(newSyntax[grep("write", newSyntax)], "write=")
  newSyntax[grep("write", newSyntax)] = paste0(tempCsv[[1]][[1]],
                                               paste0("write='H", Hclass, "c", CC, "_sol", k, ".csv'"))
  
  tempOutfile = strsplit(newSyntax[grep("outfile", newSyntax)], " ")
  tempIdx = grepl(".txt", tempOutfile[[1]])
  tempVal = character(length(tempOutfile[[1]]))
  tempVal[!tempIdx] = tempOutfile[[1]][!tempIdx]
  
  newSyntax[grep("outfile", newSyntax)] =
    paste0(paste0(tempVal[(1:length(tempVal)) < which(tempIdx)],collapse = ''),
           paste0("'H", Hclass, "c", CC, "_sol", k, ".txt'"),
           paste0(tempVal[(1:length(tempVal)) > which(tempIdx)],collapse = ''))
  
  
  newSyntax[grep("Class nominal", newSyntax)] = paste0("      Class nominal ", k, ";")
  write.table(newSyntax, paste0("LCT", CC, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}

getOutputDiCh = function(resultsTemp, Hclass, maxClassSplit1, CC = 0, Ntot = Ntot, stopCriterium = stopCriterium, itemNames = itemNames, decreasing = TRUE, ResponseOptions = ResponseOptions, Nitems = Nitems){
  
  ## What is the lowest Information Criterium?
  LL = helpFun(resultsTemp, "Log-likelihood \\(LL\\)")
  Npar = helpFun(resultsTemp, "Npar")[-1]
  # Ntot = helpFun(resultsTemp, "Number of cases")
  
  names(LL) = 1:maxClassSplit1
  BIC = -2 * LL + log(Ntot) * Npar
  AIC = -2 * LL + 2 * Npar
  AIC3 = -2 * LL + 3 * Npar
  
  IC = matrix(list(LL, BIC, AIC, AIC3), ncol = 4,
              dimnames = list(CC, c("LL", "BIC", "AIC", "AIC3")))
  
  solution = which.min(IC[[1, grep(stopCriterium, colnames(IC))]])
  
  csvTemp = read.csv(paste0("H", Hclass, "c", CC, "_sol", solution, ".csv"),
                     sep = ";", nrows = 1, skip = 8, header = FALSE, dec = ",")
  
  # Class proportions
  Classpp.temp = csvTemp[1, -c(1:5)][1:solution]
  
  order.Classpp = order(Classpp.temp, decreasing = decreasing)
  Classpp = as.matrix(sort(Classpp.temp, decreasing = decreasing))
  colnames(Classpp) = 1:solution
  rownames(Classpp) = CC
  
  # Estimated Values
  evTemp = csvTemp[1,-c(1:(5 + solution))]
  
  idxEVtemp = numeric()
  Start = 1
  End = ResponseOptions * solution
  Length = ResponseOptions * solution
  for(items in 1:Nitems){
    idxEVtemp = c(idxEVtemp, Start:End)
    Start = Start + Length * 2
    End  = End + Length * 2
  }
  
  EV = aperm(array(as.numeric(evTemp[idxEVtemp]), 
                   dim = c(ResponseOptions, solution, Nitems)),
             perm =c(2,1,3))
  
  if(solution > 1){EV = EV[order.Classpp, ,]}
  
  dimnames(EV) = list(paste0(CC, 1:solution),
                      paste0("Choice ", 1:ResponseOptions),
                      paste0("Item ", 1:Nitems))
  
  # Posteriors
  Post.temp = read.delim(paste0("H", Hclass, "c", CC, "_sol", solution, ".txt"),
                         dec = ",")
  Post.unordered = Post.temp[,(ncol(Post.temp) - (solution + 1)): (ncol(Post.temp) - 1)]
  Post = Post.unordered[c(1, order.Classpp + 1)]
  colnames(Post) = c(paste0("W_", CC), paste0("Post_", CC, 1:solution))
  
  ICtemp = Map(helpFun,
               x = list(resultsTemp),
               greplIdx = list("Entr",
                               "Classification log-likelihood",
                               "CLC",
                               "AWE",
                               "ICL-BIC",
                               "L-squared",
                               "X-squared",
                               "Cressie-Read"),
               idx = list(seq(1, 2 * maxClassSplit1, 2),
                          1:maxClassSplit1,
                          1:maxClassSplit1,
                          1:maxClassSplit1,
                          1:maxClassSplit1,
                          1:maxClassSplit1,
                          1:maxClassSplit1,
                          1:maxClassSplit1)
  )
  
  IC = cbind(IC, matrix(ICtemp, nrow = 1,
                        dimnames = list(NULL,
                                        c("Entropy", "CL", "CLC", "AWE",
                                          "ICL-BIC", "Lsquared", "Xsquared", "CR"))))
  
  toReturn = list(IC = IC, Classpp = Classpp,  #rbind
                  Post = Post, #cbind
                  EV = EV, #abind
                  solution = solution)
  return(toReturn)
}

helpFun = function(x, greplIdx, idx = NULL) {
  x = x[grepl(greplIdx, x)]
  if (is.null(idx)) idx = seq_along(x)
  return(as.numeric(gsub(",", ".", sapply(strsplit(x, "\t")[idx], `[`, 2))))
}
