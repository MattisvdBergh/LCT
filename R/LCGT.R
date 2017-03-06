#' Estimate a Latent Class Growth Tree model
#'
#' Estimate a Latent Class Growth Tree model with Latent GOLD 5.1
#'
#'
#' @param Dataset A dataframe with the data or
#' the filepath of the data.
#' @param LG Filepath of the Latent GOLD executable
#' @param LGS Filepath of Latent GOLD syntax for a model with 1- and 2-class splits
#' @param dependent The name of the dependent variable.
#' @param levelsDependent The number of options of response options of the dependent variable.
#' Should be one for a continuous dependent variable.
#' @param independent Column names of the independent variables in the dataset
#' @param caseid Name of the variable that indicates different repsondents
#' @param weight Name of the variable with the weights. When all records are unique observations,
#' this should be 1 for every observation.
#' @param resultsName Name of a folder which will be created in the working directory and
#' contains all results by Latent GOLD.
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
LCGT = function(Dataset,
                LG,
                LGS = NULL,
                dependent = NULL,
                levelsDependent = 2,
                independent = NULL,
                caseid = NULL,
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
  vecInd = mydata[,independent]

  # The tree is summarized in two lists.
  # One with the sizes of the splits, and one with the names.
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
             namesKeepVariables = namesKeepVariables)

  # Perform LC for 1- and 2-classes
  shell(paste(LG, "LCT0.lgs", "/b"))

  ################################
  ######## Return Results ########
  ################################

  # Read results from LCT file
  resultsTemp = readLines("LCT0.lst")

  Ntot = helpFun(resultsTemp, "Number of cases")[1]
  if(minSampleSize < 1){minSampleSize = Ntot * minSampleSize}
  out = getOutputLCGT(resultsTemp, Hclass, maxClassSplit1,
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
          data.temp = utils::read.delim(paste0("H", Hclass-1,
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
          utils::write.table(mydata, paste0("mydata", Hclass, preC[iCC], ".dat"),
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
            namesKeepVariables = namesKeepVariables)

          # Perform LC for 1- and 2-classes
          shell(paste0(LG, " LCT", preC[iCC],".lgs /b"))

          ################################
          ######## Retrieve Results ######
          ################################
          resultsTemp = readLines(paste0("LCT", preC[iCC],".lst"))
          outTemp = getOutputLCGT(resultsTemp, Hclass, maxClassSplit2, CC = preC[iCC],
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

  LEV = list()
  for(idxSplit in 1:length(Splitpoints)){
    classes = apply(expand.grid(Splitpoints[idxSplit],
                                1:(unlist(splitsClasses)[unlist(splitsClasses)>1])[idxSplit]),
                    1,
                    function(x){paste0(x, collapse = "")})
    LEV[[idxSplit]] = expectedValuesGrowthOrd(parmsEV = out$Parms.time,
                                 classes = classes,
                                 vecInd = vecInd,
                                 nIndependent = nIndependent,
                                 levelsDependent = levelsDependent)
    }

  meanEV = lapply(LEV, meanProb)


  # End function ####
  results = list(treeSetup = list(Names = namesClasses,
                                  cleanNames = names.classes.clean,
                                  Splits = splitsClasses,
                                  childClasses = names.split.classes,
                                  parentClasses = Splitpoints,
                                  finalClasses = final.classes),
                 splitInfo = list(IC = out$IC,
                                  Cpp = out$Classpp,
                                  CppG = ClassppGlobal,
                                  parmsCpp = out$Parms.cpp,
                                  parmsInd = out$Parms.time,
                                  vecInd = vecInd,
                                  EV = LEV,
                                  meanEV = meanEV,
                                  Posteriors = out$Post,
                                  Ntot = Ntot))
  class(results) = "LCT"
  setwd(mainDir)
  return(results)

}




helpFun = function(x, greplIdx, idx = NULL) {
  x = x[grepl(greplIdx, x)]
  if (is.null(idx)) idx = seq_along(x)
  return(as.numeric(gsub(",", ".", sapply(strsplit(x, "\t")[idx], `[`, 2))))
}



