#' Estimate a Latent Class Tree model
#'
#' Estimate a Latent Class Tree model with Latent GOLD 5.1
#'
#' @param Dataset A dataframe with the data or
#' the filepath of the data.
#' @param LG Filepath of the Latent GOLD executable
#' @param LGS Filepath of Latent GOLD syntax for a model with 1- and 2-class splits
#' @param itemNames The names of the indicators. If this is not given, all column names of the datafile will be used.
#' @param mLevels A character vector being either "ordinal" or "continuous" to indicate the measurement level of each variable. It is required when LGS is specified.
#' @param weight Name of the variable with the weights. When all records are unique observations, this should be one for every observation.
#' @param resultsName Name of a folder which will be created in the working directory and contains all results by Latent GOLD.
#'
#' @param maxClassSplit1 Maximum size of the first split of the tree. Will be assessed with the criterion given in stopCriterium. Defaults to two.
#' @param maxClassSplit2 Maximum size of each split after the first split of the tree. Defaults to two.
#' @param decreasing Whether the ordering of classes should be decreasing or not. Defaults to TRUE.
#' @param stopCriterium Criterium to decide on a split. Can be "LL" (logLikelihood), "AIC" or "BIC".
#' @param minSampleSize Minimum sample size of a class. If this is below 1, a probability of the total sample size is used.
#'
#' @param nKeepVariables Number of variables to be kept if one wants to explore the results with external variables.
#' @param namesKeepVariables Number of variables to be kept if one wants to explore the results with external variables.
#'
#' @param sets Name of the variable with the weights. When all records are unique observations, this should be one for every observation.
#' @param iterations A character vector being either ordinal or continuous to indicate the measurement level of each variable. It is required when LGS is specified.
#'
#' @details The \code{LCT} function constructs a LCT model by sequentially estimating 2-class models with Latent GOLD 5.1.
#' This can be done automatically for standard models, but for more complex models a customized Latent GOLD syntax can be provided.
#' The model size of the root can be increased with \code{maxClassSplit1} and the remaining splits with \code{maxClassSplit2}.
#'
#' @return Results of a Latent Class Tree analysis in an object of class \code{'LCT'}, which is a named list with two named lists.
#' The first list contains information on the setup of the tree and the second list contains information on every split.
#'
#'
#' @export
LCT = function(Dataset,
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
  if(!is.data.frame(Dataset)){
    mydata = utils::read.table(Dataset, header = TRUE)
  } else {
    # or use it as a dataframe
    mydata = Dataset
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
    mLevels = ifelse(mLevelsRclasses ==  "numeric", "continuous", "ordinal")
  }
  sizeMlevels = sapply(mydata[,itemNames], function(x){length(unique(x))})

  # One list with the sizes of the splits, and one list with the names
  # Hclass tracks the level of the tree
  namesClasses = splitsClasses = list()
  Hclass = 1

  # Write the first complete dataset
  utils::write.table(mydata, file = "mydataTotal.dat",
              sep = "\t", row.names = FALSE, quote = FALSE)
  dataDir = paste0(getwd(), "/mydataTotal.dat")


  ###########################################################
  ######## Make Latent GOLD syntax and run the model ########
  ###########################################################

  # check if a syntax is predefined or make a basic syntax
  if(is.null(LGS)){
    syntax = makeNewSyntax(dataDir = dataDir,
                           itemNames = itemNames,
                           weight = weight,
                           mLevels = mLevels,
                           sets = sets,
                           iterations = iterations)
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

  out = getOutputLCT(resultsTemp,
                     Hclass,
                     maxClassSplit1,
                     Ntot = Ntot,
                     stopCriterium = stopCriterium,
                     itemNames = itemNames,
                     decreasing = decreasing,
                     mLevels = mLevels,
                     sizeMlevels = sizeMlevels)

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
          data.temp = utils::read.delim(paste0("H", Hclass-1,
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
          utils::write.table(mydata, paste0("mydata", Hclass, preC[iCC], ".dat"),
                      sep = "\t", row.names = FALSE, quote = FALSE)

          ###################################################
          ######## Make new syntax and estimate models ######
          ###################################################

          ## Update LGS
          makeSyntax(dataDir = paste0(normalizePath(getwd()),
                                      "\\mydata",
                                      Hclass, preC[iCC], ".dat"),
                     Hclass = Hclass,
                     syntax = syntax,
                     maxClassSplit1 = maxClassSplit2,
                     CC = preC[iCC],
                     nKeepVariables = nKeepVariables,
                     namesKeepVariables = namesKeepVariables)

          ## Perform LC for 1- and 2-classes
          shell(paste0(LG, " LCT", preC[iCC],".lgs /b"))

          ################################
          ######## Retrieve Results ######
          ################################

          ## Write results to LCT file
          resultsTemp = readLines(paste0("LCT", preC[iCC],".lst"))
          outTemp = getOutputLCT(resultsTemp = resultsTemp,
                                 Hclass = Hclass,
                                 maxClassSplit1 = maxClassSplit2,
                                 CC = preC[iCC],
                                 Ntot = Ntot,
                                 stopCriterium = stopCriterium,
                                 itemNames = itemNames,
                                 decreasing = decreasing,
                                 mLevels = mLevels,
                                 sizeMlevels = sizeMlevels)

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

  names.classes.clean = makeCleanNames(names = namesClasses,
                                       finalClasses = final.classes)

    # End function ####
  results = list(treeSetup = list(Names = namesClasses,
                                  cleanNames = names.classes.clean,
                                  Splits = splitsClasses,
                                  childClasses = names.split.classes,
                                  parentClasses = Splitpoints,
                                  finalClasses = final.classes),
                 splitInfo = list(IC = out$IC,
                                  EV = out$EV,
                                  Cpp = out$Classpp,
                                  CppG = ClassppGlobal,
                                  Posteriors = out$Post,
                                  Ntot = Ntot,
                                  sizeMlevels = sizeMlevels))
  class(results) = "LCT"
  setwd(mainDir)
  return(results)
}




helpFun = function(x, greplIdx, idx = NULL) {
  x = x[grepl(greplIdx, x)]
  if (is.null(idx)) idx = seq_along(x)
  return(as.numeric(gsub(",", ".", sapply(strsplit(x, "\t")[idx], `[`, 2))))
}



