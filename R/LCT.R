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
LCT = function(dataDir,
               LG,
               LGS = NULL,
               itemNames = NULL,
               mLevels = NULL,
               weight = "weight",
               resultsName = "",
               maxClassSplit1 = 2,
               maxClassSplit2 = 2,
               decreasing = TRUE,
               stopCriterium = "BIC",
               minSampleSize = 5,
               nKeepVariables = 0,
               namesKeepVariables = NULL,
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
  } else # else make a weight vector of ones and add this to the data
  {
    mydata = cbind(mydata, weight = 1)
    colWeight = ncol(mydata)
    colnames(mydata)[ncol(mydata)] = weight
  }
  weights = mydata[,colWeight]

  # check the names of the indicators
  if(is.null(itemNames)){itemNames = colnames(mydata)[1:(colWeight - 1)]}
  mydata = mydata[,c(itemNames, namesKeepVariables, weight)]

  # Create results folder and make this the working directory
  mainDir = getwd();  subDir = paste0("Results", resultsName)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))

  # Determine the measurement level of each variable (if no measurement level is given, the two options will be ordinal and continuous)
  if(is.null(mLevels)){
    mLevelsRclasses = sapply(mydata[,itemNames], class)
    mLevelsLGdependent = ifelse(mLevelsRclasses ==  "numeric", "continuous", "ordinal")
  } else {
    mLevelsLGdependent = mLevels

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
             namesKeepVariables)

  # Perform LC for 1- and 2-classes
  shell(paste(LG, "LCT0.lgs", "/b"))


  ################################
  ######## Return Results ########
  ################################

  # Read results from LCT file
  resultsTemp = readLines("LCT0.lst")

  Ntot = helpFun(resultsTemp, "Number of cases")[1]
  if(minSampleSize < 1){minSampleSize = Ntot * minSampleSize}

  out = getOutputLCT(resultsTemp, Hclass, maxClassSplit1, Ntot = Ntot,
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
            namesKeepVariables)

          ## Perform LC for 1- and 2-classes
          shell(paste0(LG, " LCT", preC[iCC],".lgs /b"))


          ################################
          ######## Retrieve Results ######
          ################################

          ## Write results to LCT file
          resultsTemp = readLines(paste0("LCT", preC[iCC],".lst"))
          outTemp = getOutputLCT(resultsTemp, Hclass, maxClassSplit2, CC = preC[iCC],
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



helpFun = function(x, greplIdx, idx = NULL) {
  x = x[grepl(greplIdx, x)]
  if (is.null(idx)) idx = seq_along(x)
  return(as.numeric(gsub(",", ".", sapply(strsplit(x, "\t")[idx], `[`, 2))))
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


