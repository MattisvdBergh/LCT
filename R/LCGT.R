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
LCGT = function(dataDir,
                LG,
                LGS = NULL,
                dependent = NULL,
                independent = NULL,
                caseid = NULL,
                decreasing = TRUE,
                maxClassSplit1 = 2,
                maxClassSplit2 = 2,
                stopCriterium = "BIC",
                resultsName = "",
                minSampleSize = 5,
                levelsDependent = 2,
                weight = "weight",
                nKeepVariables = 0,
                Covariates3step = NULL,
                sets = 16,
                iterations = 50){

  #########################################################################
  ######## Check aspects of the data and write to a results folder ########
  #########################################################################

  # Read the original data from a directory
  if(!is.data.frame(dataDir)){
    mydata = read.table(dataDir, header = TRUE)
  } else {
    # or use it as a dataframe
    mydata = dataDir
  }
  stopifnot(is.data.frame(mydata))

  # check the weight variable
  if(any(colnames(mydata) == weight)){
    colWeight = which(colnames(mydata) == weight)
  } else { # else make a weight vector of ones and add this to the data
    mydata = cbind(mydata, weight = 1)
    colWeight = ncol(mydata)
    colnames(mydata)[ncol(mydata)] = weight
  }
  weights = mydata[,colWeight]

  # Create results folder and make LCT syntax there
  mainDir = getwd();  subDir = paste0("Results", resultsName)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))

  # Measurement Levels (in)dependent variables
  nIndependent = length(independent)
  mLevelsRdependent = class(mydata[,dependent])
  mLevelsLGdependent = ifelse(mLevelsRdependent ==  "numeric", "continuous", "")
  mLevelsRindependent = sapply(mydata[,independent], class)
  mLevelsLGindependent = ifelse(mLevelsRindependent ==  "numeric", "", "numeric")

  # The tree is summarized in two lists.
  # One with the sizes of the splits, and one with the names.
  namesClasses = splitsClasses = list()
  Hclass = 1

  # Write the first complete dataset
  write.table(mydata, file = "mydataTotal.dat",
              sep = "\t", row.names = FALSE, quote = FALSE)
  dataDir = paste0(getwd(), "/mydataTotal.dat")

  ###########################################################
  ######## Make Latent GOLD syntax and run the model ########
  ###########################################################

  # check if a syntax is predefined or make a basic syntax
  if(is.null(LGS)){
    syntax = makeNewSyntaxLCGT(dataDir = dataDir,
                               weight = weight,
                               independent = independent,
                               dependent = dependent,
                               caseid = caseid,
                               mLevelsLGindependent = mLevelsLGindependent,
                               mLevelsLGdependent = mLevelsLGdependent,
                               sets = sets, iterations = iterations)
  } else {
    syntax = readLines(LGS)
  }

  makeSyntax(dataDir = dataDir,
             Hclass = Hclass,
             syntax = syntax,
             maxClassSplit1 = maxClassSplit1,
             nKeepVariables = nKeepVariables,
             Covariates3step = Covariates3step)

  # Perform LC for 1- and 2-classes
  shell(paste(LG, "LCT0.lgs", "/b"))

  ################################
  ######## Return Results ########
  ################################
  ## Read results from LCT file
  resultsTemp = readLines("LCT0.lst")

  Ntot = helpFun(resultsTemp, "Number of cases")[1]
  if(minSampleSize < 1){minSampleSize = Ntot * minSampleSize}
  out = getOutput(resultsTemp, Hclass, maxClassSplit1,
                  Ntot = Ntot, stopCriterium = stopCriterium,
                  decreasing = decreasing, levelsDependent = levelsDependent,
                  nIndependent = nIndependent)

  ## Update the tree after the first split
  splitsClasses[[Hclass]] = unname(out$solution)
  namesClasses[[1]] = "0"
  namesClasses[[2]]  = paste0("0", 1:out$solution)
  Splitpoints = 0

  #########################################################
  ######## Start a while loop to repeat the splits ########
  #########################################################

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

          ###############################
          ######## Make new Data ########
          ###############################
          data.temp = read.delim(paste0("H", Hclass-1,
                                        "c", substr(preC[iCC], 1, Hclass - 1), "_sol",
                                        splitsClasses[[length(splitsClasses)]][wPreS],
                                        ".txt"), dec = ",")
          nClassesSplit = splitsClasses[[length(splitsClasses)]][wPreS]

          if (is.factor(data.temp$weight)) {
            data.temp$weight = as.numeric(levels(data.temp$weight))[data.temp$weight]
          }
          colweight = which(colnames(data.temp) == "weight")
          EstimatedWeights = suppressWarnings(
            apply(data.temp[,paste0("Cluster.", 1:nClassesSplit)],
                  2,
                  function(x){as.numeric(as.character(gsub(",", ".", x)))}))
          order.Classpp = order(colSums(data.temp$weight * EstimatedWeights, na.rm = TRUE),
                                decreasing = decreasing)
          NewWeights = EstimatedWeights[, order.Classpp[iCC]] * data.temp$weight
          mydata = data.temp[, -c((ncol(data.temp) - nClassesSplit):(ncol(data.temp)))]
          mydata[, colweight] = NewWeights
          write.table(mydata, paste0("mydata", Hclass, preC[iCC], ".dat"),
                      sep = "\t", row.names = FALSE, quote = FALSE)

          ###################################################
          ######## Make new syntax and estimate models ######
          ###################################################
          makeSyntax(dataDir = paste0(
            normalizePath(getwd()), "\\mydata", Hclass, preC[iCC], ".dat"),
            Hclass = Hclass,
            syntax = syntax,
            maxClassSplit1 = maxClassSplit2,
            CC = preC[iCC],
            nKeepVariables = nKeepVariables,
            Covariates3step = Covariates3step)

          # Perform LC for 1- and 2-classes
          shell(paste0(LG, " LCT", preC[iCC],".lgs /b"))

          ################################
          ######## Retrieve Results ######
          ################################
          resultsTemp = readLines(paste0("LCT", preC[iCC],".lst"))
          outTemp = getOutput(resultsTemp, Hclass, maxClassSplit2, CC = preC[iCC],
                              Ntot = Ntot, stopCriterium = stopCriterium, decreasing = decreasing,
                              levelsDependent = levelsDependent, nIndependent = nIndependent)

          idxRbind = 1:3
          idxCbind = 4:5
          maxCol = max(ncol(out[[2]]), ncol(outTemp[[2]]))
          maxRow = max(nrow(out[[4]]), nrow(outTemp[[4]]))
          while(ncol(outTemp[[2]]) < maxCol){outTemp[[2]] = cbind(outTemp[[2]], NA)}
          while(ncol(out[[2]]) < maxCol){out[[2]] = cbind(out[[2]], NA)}
          while(ncol(outTemp[[3]]) < maxCol){outTemp[[3]] = cbind(outTemp[[3]], NA)}
          while(ncol(out[[3]]) < maxCol){out[[3]] = cbind(out[[3]], NA)}
          while(maxRow > nrow(outTemp[[4]])){outTemp[[4]] = rbind(outTemp[[4]], NA)}

          temp1 = Map(function(x, y) rbind(x,y), out[idxRbind], outTemp[idxRbind])
          temp2 = Map(function(x, y) cbind(x,y), out[idxCbind], outTemp[idxCbind])

          out = c(temp1, temp2)
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

  ######################################################
  ######## After while loop prepare some results #######
  ######################################################

  # ### Global class proportions can be found in ClassppGlobal
  # Classpp = out$Classpp[!apply(out$Classpp, 1, anyNA),]
  # Classpp.matrix = t(Classpp)
  # Classpp.numeric = as.numeric(Classpp.matrix)
  # Classpp.numeric.names = character()
  # for(l in 1:ncol(Classpp.matrix)){
  #   Classpp.numeric.names = c(Classpp.numeric.names,
  #                             paste0(
  #                               colnames(Classpp.matrix)[l],
  #                               1:nrow(Classpp.matrix)))
  # }
  # Missing = is.na(Classpp.numeric)
  # Classpp.numeric = Classpp.numeric[!Missing]
  # names(Classpp.numeric) = Classpp.numeric.names[!Missing]
  #
  # ClassppGlobal = Names.sub.Classes = numeric()
  # for(i in 1:length(Classpp.numeric)){
  #   if(nchar(names(Classpp.numeric)[i])>2){
  #     Names.sub.Classes = numeric()
  #     for(j in 1:(nchar(names(Classpp.numeric)[i])-2)){
  #       Names.sub.Classes [j] =
  #         substr(names(Classpp.numeric)[i], 1, nchar(names(Classpp.numeric)[i]) - j)
  #     }
  #     ClassppGlobal[i] =
  #       prod(c(Classpp.numeric[i], Classpp.numeric[Names.sub.Classes]))
  #   } else{
  #     ClassppGlobal[i] =
  #       Classpp.numeric[i]
  #   }
  # }
  # names(ClassppGlobal) = names(Classpp.numeric)

  ClassppGlobal = computeGlobalCpp(ClassProportions = out$Classpp,
                                   Splits = splitsClasses,
                                   Splitpoints = Splitpoints)

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

makeNewSyntaxLCGT = function(dataDir,
                             weight,
                             independent,
                             dependent,
                             caseid,
                             mLevelsLGindependent,
                             mLevelsLGdependent,
                             sets = 16,
                             iterations = 50){

  newSyntaxToBe = capture.output(cat(paste("
//LG5.1//
version = 5.1
infile '", dataDir,"'

model
options
   maxthreads=all;
   algorithm
      tolerance=1e-008 emtolerance=0,01 emiterations=250 nriterations=50;
   startvalues
      seed=0 sets=", sets," tolerance=1e-005 iterations=", iterations,";
   bayes
      categorical=1 variances=1 latent=1 poisson=1;
   montecarlo
      seed=0 sets=0 replicates=500 tolerance=1e-008;
   quadrature  nodes=10;
   missing  includeall;
   output
      parameters=first estimatedvalues=model
      write
variables
   caseid ", caseid,";
   caseweight ", weight,";
   dependent;
   independent;
   latent
      Cluster nominal 1;
equations
   Cluster <- 1;
end model
")))


  newSyntaxToBe[grep("\\bdependent\\b", newSyntaxToBe)] =
    capture.output(cat(paste0("   dependent ", paste(dependent, mLevelsLGdependent, collapse = ", "), ";", sep = "")))
  newSyntaxToBe[grep("\\bindependent\\b", newSyntaxToBe)] =
    capture.output(cat(paste0("   independent ", paste(independent, mLevelsLGindependent, collapse = ", "), ";", sep = "")))

  indSyntax = character()
  for(idxInd in 1:length(independent)){
    indSyntax[idxInd] = paste0("+ ", independent[idxInd], " | Cluster")
  }

  newSyntaxToBe[length(newSyntaxToBe)] = paste0(dependent, "<- 1 | Cluster",
                                                paste0(indSyntax, collapse = ""),";")

  newSyntaxToBe[length(newSyntaxToBe) + 1] = "end model"
  newSyntaxToBe[length(newSyntaxToBe) + 1] = ""
  return(newSyntaxToBe)
}

makeSyntax = function(dataDir,
                      Hclass,
                      syntax,
                      maxClassSplit1,
                      nKeepVariables,
                      Covariates3step,
                      CC = 0) {

  syntax[grep("infile", syntax)] = capture.output(cat(paste0("infile \'", dataDir, "'")))
  syntax[grep("write", syntax)] =  paste0(
    "write = 'H", Hclass, "c", CC, "_sol", 1, ".csv'
    writeestimatedvalues='ev", CC, "_sol", 1, ".txt'
    outfile 'H", Hclass, "c", CC, "_sol", 1, ".txt' classification ",
    ifelse(nKeepVariables>0,
           paste0("keep ", paste0(Covariates3step, collapse = ", ")),
           ""),";"
  )

  # lm is the number of lines of one model.
  lm = which(syntax == "end model")[1] - which(syntax == "model")[1]
  syntaxModel = syntax[which(syntax == "model") :which(syntax == "end model")]
  newSyntax = syntax

  for(nClassSplit in 2:maxClassSplit1){
    # ls = length(newSyntax)
    newSyntax = c(newSyntax, syntaxModel)
    newSyntax[grep("write", newSyntax)][nClassSplit] =  paste0(
      "write = 'H", Hclass, "c", CC, "_sol", nClassSplit, ".csv'
      writeestimatedvalues='ev", CC, "_sol", nClassSplit, ".txt'
      outfile 'H", Hclass, "c", CC, "_sol", nClassSplit, ".txt' classification ",
      ifelse(nKeepVariables>0,
             paste0("keep ", paste0(Covariates3step, collapse = ", ")),
             ""),";"
    )
    newSyntax[grep("Cluster nominal", newSyntax)][nClassSplit] = paste0("      Cluster nominal ", nClassSplit, ";")
  }

  write.table(newSyntax, paste0("LCT", CC, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}

getOutput = function(resultsTemp, Hclass, maxClassSplit1,
                     CC = 0, Ntot = Ntot, stopCriterium = stopCriterium,
                     decreasing = TRUE, levelsDependent = levelsDependent,
                     nIndependent = nIndependent){

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

  ncolCSV = max(count.fields(paste0("H", Hclass, "c", CC, "_sol", solution, ".csv"), sep = ","))

  csvTemp = read.table(paste0("H", Hclass, "c", CC, "_sol", solution, ".csv"),
                       header = FALSE, col.names = paste0("V",seq_len(ncolCSV)), sep =",", fill = TRUE)

  ### DummyFirst coding parameters
  rowParms = grep("Parameters", csvTemp[,5])
  Parms.temp = as.numeric(as.character(csvTemp[
    rowParms,-c(1:5)][
      !is.na(csvTemp[rowParms,-c(1:5)])]))

  Parms.cpp.temp = c(0, Parms.temp[1:(solution - 1)])

  r1 = matrix(1:((levelsDependent - 1)*solution), levelsDependent - 1, solution)
  r2 = matrix((1 + (levelsDependent - 1)*solution):
                ((levelsDependent - 1)*solution + nIndependent*solution),
              nIndependent, solution, byrow = TRUE)
  r12 = rbind(r1, r2)
  Parms.time.unordered = Parms.temp[-c(1:(solution - 1))][r12]
  dim(Parms.time.unordered) = dim(r12) # maak er een matrix van met de juiste dimensies
  Parms.time.unordered = rbind(0, Parms.time.unordered)

  rowEV = grep("EstimatedValues", csvTemp[,5])
  Classpp.temp = csvTemp[rowEV, -c(1:5)][1:solution]

  order.Classpp = order(Classpp.temp, decreasing = decreasing)
  Classpp = as.matrix(sort(Classpp.temp, decreasing = decreasing))
  Parms.cpp = t(as.matrix(Parms.cpp.temp[order.Classpp]))
  Parms.time = as.matrix(Parms.time.unordered[,order.Classpp])

  colnames(Parms.time) = paste0(CC, 1:solution)
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

  toReturn = list(IC = IC, Parms.cpp = Parms.cpp, Classpp = Classpp,  #rbind
                  Parms.time = Parms.time, Post = Post, #cbind
                  solution = solution)
  return(toReturn)
}

computeGlobalCpp = function(ClassProportions, Splits, Splitpoints){

  # cpp = res$ClassProportions
  cpp = ClassProportions
  cpp = t(apply(cpp, 1, as.numeric))
  # Splits = res$Splits
  # Splitpoints = res$Splitpoints

  sizeAllSplits = unlist(Splits)[unlist(Splits)>1]
  allSplitClasses = unlist(sapply(1:length(sizeAllSplits),
                                  function(i){paste0(Splitpoints[i], 1:sizeAllSplits[i])}
  ))

  cppToBe = cpp
  splitsCpp = rownames(cpp)
  names(cppToBe) = paste0(0, 1:length(cppToBe))

  for(row in 2:nrow(cpp)){
    ncharSplit = nchar(splitsCpp[row])

    oldRow = substr(splitsCpp[row], 1, ncharSplit - 1)
    newCol = as.numeric(substr(splitsCpp[row], ncharSplit, ncharSplit))
    cppToBe[row,] = cppToBe[oldRow,newCol] * cppToBe[row,]
  }
  cppGtemp = as.numeric(t(cppToBe))
  tempNames = expand.grid(splitsCpp, 1:ncol(cpp))
  names(cppGtemp) = apply(
    tempNames[order(tempNames[,1]),],
    1, paste, collapse="")

  cppGtoReturn = cppGtemp[allSplitClasses]
  return(cppGtoReturn)
}


helpFun = function(x, greplIdx, idx = NULL) {
  x = x[grepl(greplIdx, x)]
  if (is.null(idx)) idx = seq_along(x)
  return(as.numeric(gsub(",", ".", sapply(strsplit(x, "\t")[idx], `[`, 2))))
}



