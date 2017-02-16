#' Estimate a Latent Class Tree model
#'
#' Estimate a Latent Class Tree model with Latent GOLD 5.1
#'
#' @param dataDir directory of the data or a dataframe with the data
#' @param LG directory of the Latent GOLD executable
#' @param LGS directory of Latent GOLD syntax for a model with 1- and 2-class splits
#' @param decreasing Whether the ordering of classes should be decreasing or not. Defaults to TRUE.
#' @param maxClassSplit1 Maximum size of the first split of the tree. Will be assessed with the criterion given in stopCriterium. Defaults to two.
#' @param maxClassSplit2 Maximum size of each split after the first split of the tree. Defaults to two.
#' @param stopCriterium Criterium to decide on a split. Can be "LL" (logLikelihood), "AIC" or "BIC".
#' @param resultsName Name of a folder which will be created in the working directory and contains all results by Latent GOLD.
#' @param minSampleSize Minimum sample size of a class. If this is below 1, a probability of the total sample size is used.
#' @param itemNames The names of the indicators. If this is not given, the rownames of the datafile will be used.
#' @param nKeepVariables Number of variables to be kept if one wants to explore the results with external variables.
#' @param weight Name of the variable with the weights. When all records are unique observations, this should be one for every observation.
#' @param measurementLevels A character vector being either ordinal or continuous to indicate the measurement level of each variable. It is required when LGS is specified.
#'
#' @return None
#' @export
LCT = function(dataDir, LG,
               LGS = NULL,
               decreasing = TRUE,
               maxClassSplit1 = 2, maxClassSplit2 = 2,
               stopCriterium = "BIC", minSampleSize = 5,
               resultsName = "", itemNames = NULL, measurementLevels = NULL,
               nKeepVariables = 0, Covariates3step = NULL, weight = "weight",
               sets = 16, iterations = 50
){

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
  } else # else make a weight vector of ones and add this to the data
  {
    mydata = cbind(mydata, weight = 1)
    colWeight = ncol(mydata)
    colnames(mydata)[ncol(mydata)] = weight
  }
  weights = mydata[,colWeight]

  # check the names of the indicators
  if(is.null(itemNames)){itemNames = colnames(mydata)[1:(colWeight - 1)]}
  mydata = mydata[,c(itemNames, Covariates3step, weight)]

  # Create results folder and make this the working directory
  mainDir = getwd();  subDir = paste0("Results", resultsName)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))

  # Determine the measurement level of each variable (if no measurement level is given, the two options will be ordinal and continuous)
  if(is.null(measurementLevels)){
    mLevelsRclasses = sapply(mydata[,itemNames], class)
    mLevelsLGdependent = ifelse(mLevelsRclasses ==  "numeric", "continuous", "ordinal")
  } else {
    mLevelsLGdependent = measurementLevels
  }
  measurementLevels = sapply(mydata[,itemNames], function(x){length(unique(x))})
  mLevels = measurementLevels

  # Write the first complete dataset
  write.table(mydata, file = "mydataTotal.dat",
              sep = "\t", row.names = FALSE, quote = FALSE)
  dataDir = paste0(getwd(), "/mydataTotal.dat")

  # One list with the sizes of the splits, and one list with the names
  # Hclass tracks the level of the tree
  namesClasses = splitsClasses = list()
  Hclass = 1

  ###########################################################
  ######## Make Latent GOLD syntax and run the model ########
  ###########################################################

  # check if a syntax is predefined or make a basic syntax
  if(is.null(LGS)){
    syntax = makeNewSyntax(dataDir, itemNames, weight, mLevelsLGdependent,
                           sets = sets, iterations = iterations)
  } else {
    syntax = readLines(LGS)
  }

  # Make a syntax for the root of the tree
  makeSyntax(dataDir,
             Hclass,
             syntax,
             maxClassSplit1,
             nKeepVariables,
             Covariates3step)

  # Perform LC for 1- and 2-classes
  shell(paste(LG, "LCT0.lgs", "/b"))


  ################################
  ######## Return Results ########
  ################################

  # Read results from LCT file
  resultsTemp = readLines("LCT0.lst")

  Ntot = helpFun(resultsTemp, "Number of cases")[1]
  if(minSampleSize < 1){minSampleSize = Ntot * minSampleSize}

  out = getOutput(resultsTemp, Hclass, maxClassSplit1, Ntot = Ntot,
                  stopCriterium = stopCriterium, itemNames = itemNames,
                  decreasing = decreasing, mLevelsLGdependent = mLevelsLGdependent,
                  mLevels = mLevels, measurementLevels = measurementLevels)

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
    sizePreviousSplits = 1:splitsClasses[[length(splitsClasses)]][1]
    nPreviousSplits = 1:length(splitsClasses[[length(splitsClasses)]])

    # namesClassesTemp will contain the new names of the split classes and will be later added to namesClasses
    # splitsClassesTemp will contain the split sizes of the new splits and will be later added to splitsClasses
    namesClassesTemp = splitsClassesTemp = numeric()

    # At each hierarchical level, we have to run through the number of splits performed at the previous level.
    # wPreS is which previous split
    for (wPreS in nPreviousSplits) {

      # The size of the previous split must updated
      if (wPreS != 1) {
        sizePreviousSplits = sizePreviousSplits +
          splitsClasses[[length(splitsClasses)]][wPreS - 1]
      }
      # Correction if we switch number of classes
      length(sizePreviousSplits) = splitsClasses[[length(splitsClasses)]][wPreS]
      if (anyNA(sizePreviousSplits)) {
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

          if (is.factor(data.temp$weight)) {
            data.temp$weight = as.numeric(levels(data.temp$weight))[data.temp$weight]
          }
          colweight = which(colnames(data.temp) == "weight")
          EstimatedWeights = suppressWarnings(
            apply(data.temp[(colweight + 1 + nKeepVariables):(ncol(data.temp)-1)], 2,
                  function(x){as.numeric(as.character(gsub(",", ".", x)))}))
          order.Classpp = order(colSums(data.temp$weight * EstimatedWeights, na.rm = TRUE),
                                decreasing = decreasing)
          NewWeights = EstimatedWeights[, order.Classpp[iCC]] * data.temp$weight

          dataPattern = data.temp[, 1:(colweight-1)]
          if(nKeepVariables > 0){
            keepVariables = data.temp[,(colweight+1):(colweight + nKeepVariables)]
            mydata = cbind(dataPattern, keepVariables, weight = NewWeights)
          } else {
            mydata = cbind(data.temp[, 1:(colweight-1)], weight = NewWeights)
          }
          write.table(mydata, paste0("mydata", Hclass, preC[iCC], ".dat"),
                      sep = "\t", row.names = FALSE, quote = FALSE)

          ###################################################
          ######## Make new syntax and estimate models ######
          ###################################################

          ## Update LGS
          makeSyntax(dataDir = paste0(
            normalizePath(getwd()), "\\mydata", Hclass, preC[iCC], ".dat"),
            Hclass,
            syntax,
            maxClassSplit2,
            CC = preC[iCC],
            nKeepVariables,
            Covariates3step)

          ## Perform LC for 1- and 2-classes
          shell(paste0(LG, " LCT", preC[iCC],".lgs /b"))


          ################################
          ######## Retrieve Results ######
          ################################

          ## Write results to LCT file
          resultsTemp = readLines(paste0("LCT", preC[iCC],".lst"))
          outTemp = getOutput(resultsTemp, Hclass, maxClassSplit2, CC = preC[iCC],
                              Ntot = Ntot, stopCriterium = stopCriterium,
                              itemNames = itemNames, decreasing = decreasing,
                              mLevelsLGdependent = mLevelsLGdependent, mLevels = mLevels,
                              measurementLevels = measurementLevels)

          idxRbind = 1:2
          idxCbind = 3:4
          maxCol = max(ncol(out[[2]]), ncol(outTemp[[2]]))
          maxRow = max(nrow(out[[4]]), nrow(outTemp[[4]]))
          while(ncol(outTemp[[2]]) < maxCol){outTemp[[2]] = cbind(outTemp[[2]], NA)}
          while(ncol(outTemp[[3]]) < maxCol){outTemp[[3]] = cbind(outTemp[[3]], NA)}
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

  # # Global class proportions can be found in ClassppGlobal
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
  results = list(namesClasses, out$IC, splitsClasses, out$EV, out$Classpp, out$Post,
                 Splitpoints, ClassppGlobal, names.split.classes, final.classes,
                 Ntot, names.classes.clean)
  names(results) = c("Names", "IC", "Splits", "EV", "ClassProportions",
                     "Posteriors", "Splitpoints", "ClassppGlobal", "namesSplitClasses", "finalClasses",
                     "Ntot", "Names.clean")
  class(results) = "LCT"
  setwd(mainDir)
  return(results)

}



makeNewSyntax = function(dataDir, itemNames, weight, mLevelsLGdependent,
                         sets = 16, iterations = 50){

  newSyntaxToBe = capture.output(cat(paste("
//LG5.1//
version = 5.1
infile '", dataDir,"'

model
options
   maxthreads=2;
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
   caseweight ", weight,";
   dependent;
   latent
      Cluster nominal 1;
equations
   Cluster <- 1;
end model
")))


  newSyntaxToBe[grep("dependent", newSyntaxToBe)] =
    capture.output(cat(paste0("   dependent ", paste(itemNames, mLevelsLGdependent, collapse = ", "), ";", sep = "")))

  newSyntaxToBe[length(newSyntaxToBe)] = paste0("   ", itemNames[1]," <- 1 + Cluster;")
  for(idxVar in 2:length(itemNames)){
    newSyntaxToBe[length(newSyntaxToBe) + 1] = paste0("   ", itemNames[idxVar]," <- 1 + Cluster;")
  }

  for(idxVarCon in which(mLevelsLGdependent == "continuous")){
    newSyntaxToBe[length(newSyntaxToBe) + 1] = paste0("   ", itemNames[idxVarCon]," | Cluster;")
  }

  newSyntaxToBe[length(newSyntaxToBe) + 1] = "end model"
  newSyntaxToBe[length(newSyntaxToBe) + 1] = ""
  return(newSyntaxToBe)
}

makeSyntax = function(dataDir, Hclass, syntax, maxClassSplit1,
                      nKeepVariables, Covariates3step, CC = 0) {

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

getOutput = function(resultsTemp, Hclass, maxClassSplit1, CC = 0,
                     Ntot = Ntot, stopCriterium = stopCriterium,
                     itemNames = itemNames, decreasing = TRUE,
                     mLevelsLGdependent = mLevelsLGdependent, mLevels = mLevels,
                     measurementLevels = measurementLevels){

  ## What is the lowest Information Criterium?
  LL = helpFun(resultsTemp, "Log-likelihood \\(LL\\)")
  Npar = helpFun(resultsTemp, "Npar")[-1]

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

  rowEV = grep("EstimatedValues", csvTemp[,5])

  # Class proportions
  Classpp.temp = csvTemp[rowEV, -c(1:5)][1:solution]
  order.Classpp = order(Classpp.temp, decreasing = decreasing)
  Classpp = as.matrix(Classpp.temp[order.Classpp])
  colnames(Classpp) = 1:solution
  rownames(Classpp) = CC

  # Estimated Values
  evTemp = csvTemp[rowEV,-c(1:(5 + solution))]
  evTemp = evTemp[!is.na(evTemp)]

  ordVar = which(mLevelsLGdependent == "ordinal")
  conVar = which(mLevelsLGdependent == "continuous")

  count = 1
  EVtemp2 = matrix(, nrow = 0, ncol = solution)
  for(idxVar in seq_along(itemNames)){
    if(idxVar %in% ordVar){
      EVtemp = matrix(evTemp[count : (count + (mLevels[idxVar] * solution) - 1)],
                      ncol = solution)
      rownames(EVtemp) = paste0(itemNames[idxVar], ".", 1:measurementLevels[[idxVar]])
      EVtemp2 = rbind(EVtemp2, EVtemp)

      count = count + (mLevels[idxVar] * solution)
    }
    if(idxVar %in% conVar){
      EVtemp = matrix(evTemp[count : (count + solution - 1)],
                      ncol = solution)
      rownames(EVtemp) = paste0(itemNames[idxVar])
      EVtemp2 = rbind(EVtemp2, EVtemp)
      count = count + solution
    }
  }

  # Sort if needed
  if(solution > 1){
    EV = EVtemp2[,order.Classpp]
  } else {
    EV = EVtemp2
  }

  colnames(EV) = paste0(CC, 1:solution)

  # Posteriors
  Post.temp = read.delim(paste0("H", Hclass, "c", CC, "_sol", solution, ".txt"),
                         dec = ",")
  Post.unordered = Post.temp[,(ncol(Post.temp) - (solution + 1)): (ncol(Post.temp) - 1)]
  Post = Post.unordered[c(1, order.Classpp + 1)]
  colnames(Post) = c(paste0("W_", CC), paste0("Post_", CC, 1:solution))

  if(any(mLevelsLGdependent == "continuous")){
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
                                          c("Entropy",
                                            "CL",
                                            "CLC",
                                            "AWE",
                                            "ICL-BIC"))))

  } else{ ICtemp = Map(helpFun,
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
                                        c("Entropy", "CL",
                                          "CLC",
                                          "AWE",
                                          "ICL-BIC",
                                          "Lsquared",
                                          "Xsquared",
                                          "CR"))))

  }


  toReturn = list(IC = IC, Classpp = Classpp,  #rbind
                  Post = Post, EV = EV, #cbind
                  solution = solution)
  return(toReturn)
}

helpFun = function(x, greplIdx, idx = NULL) {
  x = x[grepl(greplIdx, x)]
  if (is.null(idx)) idx = seq_along(x)
  return(as.numeric(gsub(",", ".", sapply(strsplit(x, "\t")[idx], `[`, 2))))
}

relImpr = function(res, criterion = "LL"){

  allImprovement = matrix(, nrow = length(res$IC[,criterion[1]][[1]]) - 2,
                          ncol = length(criterion))

  for(idxCrit in seq_along(criterion)){
    nClass = length(res$IC[,criterion[idxCrit]][[1]])
    improvement = numeric()
    counter = 1
    for(i in 3:nClass){
      improvement[counter] = (res$IC[,criterion[idxCrit]][[1]][i] -
                                res$IC[,criterion[idxCrit]][[1]][i - 1])/
        (res$IC[,criterion[idxCrit]][[1]][2] -
           res$IC[,criterion[idxCrit]][[1]][1])
      counter = counter + 1
    }
    allImprovement[,idxCrit] = improvement
  }
  colnames(allImprovement) = criterion
  return(allImprovement)
}

makePPperLevel = function(res, PPname = ""){
  EV = apply(res$EV, 2, as.numeric)
  rownames(EV) = rownames(res$EV)
  for(idxLevels in 2:length(res$Names.clean)){
    classesTemp = res$Names.clean[[idxLevels]]
    splits = res$Splits[[idxLevels - 1]]
    xaxes = list()
    counter = 1
    for(idxSplits in seq_along(splits)){
      xaxes[[idxSplits]] = c(counter:(splits[idxSplits] + counter - 1))
      counter = max(unlist(xaxes)) + 1
    }
    EVtemp = EV[,classesTemp]

    pdf(paste0(PPname, "plot", idxLevels, ".pdf"))
    plot(unlist(xaxes), rep(-1, sum(lengths(xaxes))),
         ylim = range(EVtemp), xaxt = "n", xlab = "", las = 2,
         ylab = "", axes = FALSE, type = "b")

    for(idxSplits in xaxes){
      for(rows in 1:nrow(EV)){
        lines(idxSplits, EV[rows, classesTemp[idxSplits]],
              pch = rows, col = rows, type = "b")}
    }
    axis(1, at = unlist(xaxes), labels = substr(classesTemp, 2, nchar(classesTemp)),
         col = "white", cex.axis = 1.5)
    axis(2, las = 2, cex.axis = 1.5)
    dev.off()
  }
  pdf(paste0(PPname, "legend.pdf"))
  frame()
  legend("top",
         legend = rownames(EV),
         pch = 1:nrow(EV),
         col = 1:nrow(EV),
         bty = "n",
         cex = 1,
         xpd = TRUE)
  dev.off()
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


