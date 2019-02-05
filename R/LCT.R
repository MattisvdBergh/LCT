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
               iterations = 50,
               dec = ",",
               sep = ";"){

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
  sizeMlevels = sapply(mydata[,itemNames], function(x){length(unique(x[!is.na(x)]))})

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
                     sizeMlevels = sizeMlevels,
                     dec = dec,
                     sep = sep)

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
                                        ".txt"), dec = dec)

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
                                 sizeMlevels = sizeMlevels,
                                 dec = dec,
                                 sep = sep)

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



computeGlobalCpp = function(ClassProportions, Splits, Splitpoints){

  cpp = ClassProportions
  cpp = t(apply(cpp, 1, as.numeric))

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

makeNewSyntax = function(dataDir,
                         itemNames,
                         weight,
                         mLevels,
                         sets = 16,
                         iterations = 50){

  newSyntaxToBe = utils::capture.output(cat(paste("

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
parameters=first estimatedvalues=model reorderclasses
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
    utils::capture.output(cat(paste0("   dependent ", paste(itemNames, mLevels, collapse = ", "), ";", sep = "")))

  newSyntaxToBe[length(newSyntaxToBe)] = paste0("   ", itemNames[1]," <- 1 + Cluster;")
  for(idxVar in 2:length(itemNames)){
    newSyntaxToBe[length(newSyntaxToBe) + 1] = paste0("   ", itemNames[idxVar]," <- 1 + Cluster;")
  }

  for(idxVarCon in which(mLevels == "continuous")){
    newSyntaxToBe[length(newSyntaxToBe) + 1] = paste0("   ", itemNames[idxVarCon]," | Cluster;")
  }

  newSyntaxToBe[length(newSyntaxToBe) + 1] = "end model"
  newSyntaxToBe[length(newSyntaxToBe) + 1] = ""
  return(newSyntaxToBe)
}

makeSyntax = function(dataDir,
                      Hclass,
                      syntax,
                      maxClassSplit1,
                      nKeepVariables,
                      namesKeepVariables,
                      CC = 0) {

  syntax[grep("infile", syntax)] = utils::capture.output(cat(paste0("infile \'", dataDir, "'")))
  syntax[grep("write", syntax)] =  paste0(
    "write = 'H", Hclass, "c", CC, "_sol", 1, ".csv'
    writeestimatedvalues='ev", CC, "_sol", 1, ".txt'
    outfile 'H", Hclass, "c", CC, "_sol", 1, ".txt' classification ",
    ifelse(nKeepVariables>0,
           paste0("keep ", paste0(namesKeepVariables, collapse = ", ")),
           ""),";"
  )

  # lm is the number of lines of one model.
  lm = which(syntax == "end model")[1] - which(syntax == "model")[1]
  syntaxModel = syntax[which(syntax == "model") :which(syntax == "end model")]
  newSyntax = syntax

  for(nClassSplit in 2:maxClassSplit1){
    newSyntax = c(newSyntax, syntaxModel)
    newSyntax[grep("write", newSyntax)][nClassSplit] =  paste0(
      "write = 'H", Hclass, "c", CC, "_sol", nClassSplit, ".csv'
      writeestimatedvalues='ev", CC, "_sol", nClassSplit, ".txt'
      outfile 'H", Hclass, "c", CC, "_sol", nClassSplit, ".txt' classification ",
      ifelse(nKeepVariables>0,
             paste0("keep ", paste0(namesKeepVariables, collapse = ", ")),
             ""),";"
    )
    newSyntax[grep("Cluster nominal", newSyntax)][nClassSplit] = paste0("      Cluster nominal ", nClassSplit, ";")
  }

  utils::write.table(newSyntax, paste0("LCT", CC, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}

getOutputLCT = function(resultsTemp,
                        Hclass,
                        maxClassSplit1,
                        CC = 0,
                        Ntot = Ntot,
                        stopCriterium = stopCriterium,
                        itemNames = itemNames,
                        decreasing = TRUE,
                        mLevels = mLevels,
                        sizeMlevels = sizeMlevels,
                        dec = dec,
                        sep = sep){

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

  ncolCSV = max(utils::count.fields(paste0("H", Hclass, "c", CC, "_sol", solution, ".csv"), sep = sep))

  csvTemp = utils::read.table(paste0("H", Hclass, "c", CC, "_sol", solution, ".csv"),
                              header = FALSE, col.names = paste0("V",seq_len(ncolCSV)),
                              sep = sep, fill = TRUE, dec = dec)

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

  ordVar = which(mLevels == "ordinal")
  conVar = which(mLevels == "continuous")

  count = 1
  EVtemp2 = matrix(, nrow = 0, ncol = solution)
  for(idxVar in seq_along(itemNames)){
    if(idxVar %in% ordVar){
      EVtemp = matrix(evTemp[count : (count + (sizeMlevels[idxVar] * solution) - 1)],
                      ncol = solution)
      rownames(EVtemp) = paste0(itemNames[idxVar], ".", 1:sizeMlevels[[idxVar]])
      EVtemp2 = rbind(EVtemp2, EVtemp)

      count = count + (sizeMlevels[idxVar] * solution)
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
  Post.temp = utils::read.delim(paste0("H", Hclass, "c", CC, "_sol", solution, ".txt"),
                                dec = dec)
  colWeight = which(colnames(Post.temp) == "weight")
  Post.unordered = Post.temp[,c(colWeight, (ncol(Post.temp) - (solution)): (ncol(Post.temp) - 1))]
  Post = Post.unordered[c(1, order.Classpp + 1)]
  colnames(Post) = c(paste0("W_", CC), paste0("Post_", CC, 1:solution))

  if(any(mLevels == "continuous")){
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

makeCleanNames = function(names, finalClasses){

  names.classes.clean = names
  for(k in 2:length(names)){
    row.classes = names[[k]]
    classbranches = matrix(, nrow = k - 1, ncol = length(row.classes))
    for (j in 1:length(row.classes)){
      classbranch = character()
      for(i in 2:k){classbranch[i - 1] = substr(row.classes[j], 1, i)}
      classbranches[,j] = classbranch}
    classbranchesTF = apply(classbranches, 2, function(x){x%in%finalClasses})
    if(k == 2){
      IndexChangedClasses = which(classbranchesTF)
    } else {IndexChangedClasses = which(classbranchesTF, arr.ind = TRUE)[,2]}
    names.classes.clean[[k]][IndexChangedClasses] = classbranches[which(classbranchesTF, arr.ind = TRUE)]
  }
  return(names.classes.clean)
}

