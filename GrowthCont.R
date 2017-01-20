#' Illustration of crayon colors
#'
#' Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param method2order method to order colors (\code{"hsv"} or \code{"cluster"})
#' @param cex character expansion for the text
#' @param mar margin paramaters; vector of length 4 (see \code{\link[graphics]{par}})
#'
#' @return None
#' @export
LCGT = function(dataDir, LG, LGS,
                decreasing = TRUE, maxClassSplit1 = 2,
                maxClassSplit2 = 2, stopCriterium = "BIC",
                resultsName = "", minSampleSize = 5,
                levelsDependent = 2, nIndependent = 2, nCov = 0){

  helpFun = function(x, greplIdx, idx = NULL) {
    x = x[grepl(greplIdx, x)]
    if (is.null(idx)) idx = seq_along(x)
    return(as.numeric(gsub(",", ".", sapply(strsplit(x, "\t")[idx], `[`, 2))))
  }

  # Read the original data
  mydata = read.table(dataDir, header = TRUE)
  weight = mydata$weight

  # The tree is summarized in two lists.
  # One with the sizes of the splits, and one with the names.
  namesClasses = splitsClasses = list()
  Hclass = 1

  # Create results folder and make LCT syntax there
  mainDir = getwd();  subDir = paste0("Results", resultsName)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  
  syntax = readLines(LGS)
  makeSyntax(dataDir, Hclass, syntax, maxClassSplit1)
  
  # Perform LC for 1- and 2-classes ####
  shell(paste(LG, "LCT0.lgs", "/b"))
  
  # foo = function(...){print(1)}
  # obj = 'case1'
  # obj = 'case2'
  # output = switch(obj,
  #        case1 = foo(),
  #        case2 = print(1:10),
  #        case3 = 3)
  
  ## Write results to LCT file
  resultsTemp = readLines("LCT0.lst")
  
  Ntot = helpFun(resultsTemp, "Number of cases")[1]
  if(minSampleSize < 1){minSampleSize = Ntot * minSampleSize}
  out = getOutput(resultsTemp, Hclass, maxClassSplit1, 
                  Ntot = Ntot, stopCriterium = stopCriterium, 
                  decreasing = decreasing, levelsDependent = levelsDependent, 
                  nIndependent = nIndependent, nCov = nCov)
  
  ## Update the tree after the first split
  splitsClasses[[Hclass]] = unname(out$solution)

  classSizes = Ntot * as.numeric(as.character(out$Classpp[1,]))
  classSizes = classSizes[!is.na(classSizes)]
  if(any((classSizes < minSampleSize) == TRUE)){
    solution = 1
  }

  namesClasses[[1]] = "0"
  namesClasses[[2]]  = paste0("0", 1:out$solution)
  Splitpoints = 0
  
  while (!all(splitsClasses[[length(splitsClasses)]] == 1)) {
    
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
    # wPreS is which previous split
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
          
          # temp.fun = function(x, iCC) paste0("H", Hclass-1, "c", substr(iCC, 1, Hclass - 1), "_sol", x, ".txt")
          # files = apply(as.matrix(splitsClasses[[length(splitsClasses)]]), 1, temp.fun, iCC)
          # data.temp = lapply(files, function(x) read.delim(file = x, dec = ","))[[1]]
          
          if (is.factor(data.temp$weight)) {
            data.temp$weight = as.numeric(levels(data.temp$weight))[data.temp$weight]
          }
          colweight = which(colnames(data.temp) == "weight")
          EstimatedWeights = suppressWarnings(
            apply(data.temp[(colweight + 1):(ncol(data.temp)-1)], 2,
                  function(x){as.numeric(as.character(gsub(",", ".", x)))}))
          order.Classpp = order(colSums(data.temp$weight * EstimatedWeights, na.rm = TRUE),
                                decreasing = decreasing)
          NewWeights = EstimatedWeights[, order.Classpp[iCC]] * data.temp$weight
          mydata = cbind(data.temp[, 1:(colweight-1)], weight = NewWeights)
          write.table(mydata, paste0("mydata", Hclass, preC[iCC], ".dat"),
                      sep = "\t", row.names = FALSE, quote = FALSE)
          
          ## Update LGS ####
          makeSyntax(dataDir = paste0(
            normalizePath(getwd()), "\\mydata", Hclass, preC[iCC], ".dat"),
            Hclass = Hclass,
            syntax = syntax,
            maxClassSplit1 = maxClassSplit2,
            CC = preC[iCC])
          
          ## Perform LC for 1- and 2-classes ####
          shell(paste0(LG, " LCT", preC[iCC],".lgs /b"))
          
          ## Write results to LCT file
          resultsTemp = readLines(paste0("LCT", preC[iCC],".lst"))
          outTemp = getOutput(resultsTemp, Hclass, maxClassSplit2, CC = preC[iCC], 
                              Ntot = Ntot, stopCriterium = stopCriterium, decreasing = decreasing,
                              levelsDependent = levelsDependent, nIndependent = nIndependent,
                              nCov = nCov)
          
          idxRbind = 1:4
          idxCbind = 5:9
          maxCol = max(ncol(out[[2]]), ncol(outTemp[[2]]))
          maxRow = max(nrow(out[[4]]), nrow(outTemp[[4]]))
          maxRow2 = max(nrow(out[[7]]), nrow(outTemp[[7]]))
          while(ncol(outTemp[[2]]) < maxCol){outTemp[[2]] = cbind(outTemp[[2]], NA)}
          while(ncol(out[[2]]) < maxCol){out[[2]] = cbind(out[[2]], NA)}
          while(ncol(outTemp[[3]]) < maxCol){outTemp[[3]] = cbind(outTemp[[3]], NA)}
          while(ncol(out[[3]]) < maxCol){out[[3]] = cbind(out[[3]], NA)}
          while(maxCol > ncol(outTemp[[4]])){outTemp[[4]] = cbind(outTemp[[4]], NA)}
          while(maxRow2 > nrow(outTemp[[7]])){outTemp[[7]] = rbind(outTemp[[7]], NA)}
          
          temp1 = Map(function(x, y) rbind(x,y), out[idxRbind], outTemp[idxRbind])
          temp2 = Map(function(x, y) cbind(x,y), out[idxCbind], outTemp[idxCbind])
          
          out = c(temp1, temp2)
          solution = outTemp$solution
          
          sizeParentClass = sum(as.numeric(as.character(outTemp$Post[,1])), na.rm = TRUE)/4
          classSizes =  sizeParentClass * as.numeric(as.character(outTemp$Classpp[1,]))
          classSizes = classSizes[!is.na(classSizes)]
          if(any((classSizes < minSampleSize) == TRUE)){
            solution = 1
          }
          
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
                              paste0(
                                colnames(Classpp.matrix)[l],
                                1:nrow(Classpp.matrix)))
  }
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
  results = list(namesClasses, out$IC, splitsClasses, out$Parms.cpp, out$Classpp, out$Parms.time,
                 out$Post, Splitpoints, ClassppGlobal, names.split.classes, final.classes,
                 Ntot, names.classes.clean)
  names(results) = c("Names", "IC", "Splits", "ParametersClassProportions", "ClassProportions", "ParametersTimepoints",
                     "Posteriors", "Splitpoints", "ClassppGlobal", "Names.Split.Classes", "Final.Classes",
                     "Ntot", "Names.clean")
  setwd(mainDir)
  return(results)
  
}

makeSyntax = function(dataDir, Hclass, syntax, maxClassSplit1, CC = 0) {
  
  syntax[grep("infile", syntax)] = capture.output(cat(paste0("infile \'", dataDir, "'")))
  
  # lm is the number of lines of one model.
  # bm is the beginning where the next model should start
  # em is the end where the next model should end.
  k = 1:maxClassSplit1
  
  lm = which(syntax == "end model")[2] - which(syntax == "model")[2] + 1
  bm = 6 + (k - 1) * lm
  em = 5 + k * lm
  
  newSyntax = character(em[length(em)])
  
  newSyntax[1:5] = syntax[1:5]
  newSyntax[6:em[length(em)]] = rep(syntax[6:em[1]], maxClassSplit1)
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


getOutput = function(resultsTemp, Hclass, maxClassSplit1, 
                     CC = 0, Ntot = Ntot, stopCriterium = stopCriterium, 
                     decreasing = TRUE, levelsDependent = levelsDependent,  nIndependent = nIndependent,
                     nCov = nCov){
  
  helpFun = function(x, greplIdx, idx = NULL) {
    x = x[grepl(greplIdx, x)]
    if (is.null(idx)) idx = seq_along(x)
    return(as.numeric(gsub(",", ".", sapply(strsplit(x, "\t")[idx], `[`, 2))))
  }
  
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
  if(solution>1){
  csvTemp = read.csv(paste0("H", Hclass, "c", CC, "_sol", solution, ".csv"),
                     sep = ",", nrows = 3, skip = 7, header = FALSE, dec = ",", na.strings = "")
  
  ##### covariates
  Parms.temp = as.numeric(as.character(csvTemp[1,-c(1:5)][!is.na(csvTemp[1,-c(1:5)])]))
  SE.temp = as.numeric(as.character(unlist(csvTemp[3, -c(1:5)])))
  
  Classpp.temp = csvTemp[2, -c(1:5)][1:solution]
  SE.Classpp.temp = SE.temp[1:(solution - 1)]
  Parms.cpp.temp = c(0, Parms.temp[1:(solution - 1)])
  
  Parms.tempNoCpp = Parms.temp[-c(1:(solution - 1))]
  
  if(nCov > 0){
    covariates.temp = cbind(0, matrix(Parms.tempNoCpp[1:(nCov* (solution - 1))], ncol = solution - 1, byrow = TRUE))
    SE.cov.temp = cbind(0, matrix(SE.temp[1:(nCov* (solution - 1))], ncol = solution - 1, byrow = TRUE))
    Parms.tempNoCov = Parms.tempNoCpp[-c(1:(nCov*(solution - 1)))]
    SE.tempNoCov = SE.temp[-c(1:(nCov*(solution - 1)+(solution - 1)))]
  }
  
  lp = length(Parms.tempNoCov)
  Parms.time.unordered = matrix(Parms.tempNoCov[-c((lp - solution + 1):lp)], ncol = solution, byrow = TRUE)
  SE.time.unordered = matrix(SE.tempNoCov[-c((lp - solution + 1):lp)], ncol = solution, byrow = TRUE)
  
  order.Classpp = order(Classpp.temp, decreasing = decreasing)
  Classpp = as.matrix(sort(Classpp.temp, decreasing = decreasing))
  Parms.cpp = t(as.matrix(Parms.cpp.temp[order.Classpp]))
  Parms.time = as.matrix(Parms.time.unordered[,order.Classpp])
  SE.cpp = matrix(c(0, as.numeric(as.character(unlist(SE.Classpp.temp))))[order.Classpp], 1)
  SE.time = SE.time.unordered[, order.Classpp]
  covariates  = covariates.temp[,order.Classpp]
  SE.cov = SE.cov.temp[,order.Classpp]

  colnames(covariates) = colnames(SE.cov) = colnames(SE.time) = colnames(Parms.time) = paste0(CC, 1:solution)
  colnames(Classpp) = colnames(Parms.cpp) = 1:solution
  rownames(Classpp) = rownames(Parms.cpp) = CC

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
                               "ICL-BIC"),
               idx = list(seq(1, 2 * maxClassSplit1, 2),
                          1:maxClassSplit1,
                          1:maxClassSplit1,
                          1:maxClassSplit1,
                          1:maxClassSplit1)
  )

  IC = cbind(IC, matrix(ICtemp, nrow = 1,
                        dimnames = list(NULL,
                                        c("Entropy", "CL", "CLC", "AWE",
                                          "ICL-BIC"))))

  toReturn = list(IC = IC, Parms.cpp = Parms.cpp, Classpp = Classpp,  SE.cpp = SE.cpp, #rbind
                  Parms.time = Parms.time, Post = Post, SE.time = SE.time, covariates = covariates, SE.cov = SE.cov, #cbind
                  solution = solution)} else { toReturn = list(IC = matrix(NA, ncol = 9, nrow = 1), Parms.cpp = matrix(NA, 1, 1), Classpp = matrix(NA, 1, 1),  SE.cpp= matrix(NA, 1, 1), #rbind
                                    Parms.time = matrix(NA, 4, 1), Post = matrix(NA, nrow = length(weight)), SE.time = matrix(NA, 4, 1), covariates = matrix(NA, nCov, 1), SE.cov = matrix(NA, nCov, 1), #cbind
                                    solution = solution)}
  return(toReturn)
}

