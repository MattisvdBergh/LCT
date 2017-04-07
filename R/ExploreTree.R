#' 3-step method Latent Class Trees
#'
#' @param resTree A Latent Class Tree object
#' @param LG Path to the Latent GOLD executable
#' @param syntaxExpl Directory of Latent GOLD syntax for a model for the 3-step analysis. This can be used for more complex models, but
#' @param measurementLevels Whether the ordering of classes should be decreasing or not. Defaults to TRUE.
#' @param dirTreeResults Maximum size of the first split of the tree. Will be assessed with the criterion given in stopCriterium. Defaults to two.
#' @param ResultsFolder folder with results
#'
#' @description Compute the relation of some covariates with each class by use of the 3-step method.
#'
#' @return None
#' @export
exploreTree = function(resTree,
                       LG = LG,
                       syntaxExpl = NULL,
                       sizeMlevels = c(2, 1),
                       method = "bch",
                       dirTreeResults = getwd(),
                       ResultsFolder = "exploreTree",
                       Covariates = c("age", "sex"),
                       posCovVal = NULL,
                       mLevels = c("ordinal", "continuous"),
                       weight = "weight",
                       analysis = "dependent"){
  mainDir = getwd()

  pSplits = resTree$treeSetup$parentClasses
  Splits = resTree$treeSetup$Splits
  namesClasses = resTree$treeSetup$Names

  allSplits = unlist(Splits)
  splitsLogical = allSplits>1
  splitSizes = allSplits[splitsLogical]
  hClasses = nchar(pSplits)
  names(splitSizes) = pSplits

  # All files with posteriors
  dirPost = paste0("H", hClasses,
                   "c", pSplits,
                   "_sol", splitSizes, ".txt")

  dir.create(file.path(getwd(), ResultsFolder), showWarnings = FALSE)
  setwd(file.path(mainDir, ResultsFolder))

  if(is.null(syntaxExpl)){
    syntaxExpl = makeNewSyntaxExplore(dataDir = dirPost[1],
                                      method = method,
                                      mLevels = mLevels,
                                      Covariates = Covariates,
                                      weight = weight,
                                      analysis = analysis)
  }

  Results3step = run3step(dirPost = dirPost,
                          syntaxExpl = syntaxExpl,
                          splitSizes = splitSizes,
                          pSplits = pSplits,
                          dirTreeResults = dirTreeResults,
                          Covariates = Covariates,
                          analysis = analysis,
                          posCovVal = posCovVal)

  toReturn = reformResults(splitSizes = splitSizes,
                           ParmsTotal = Results3step$ParmsTotal,
                           ProfileTotal = Results3step$ProfileTotal,
                           EV = Results3step$EV,
                           WaldTotal = Results3step$WaldTotal,
                           sizeMlevels = sizeMlevels,
                           mLevels = mLevels,
                           Covariates = Covariates,
                           pSplits = pSplits,
                           analysis = analysis)
  setwd(mainDir)
  return(toReturn)
}

makeNewSyntaxExplore = function(dataDir = dataDir,
                                method = "bch",
                                mLevels = mLevels,
                                Covariates = Covariates,
                                weight = "weight",
                                analysis = "dependent"){

  syntaxExploreToBe = utils::capture.output(cat(paste("
//LG5.1//
version = 5.1
infile '", dataDir,"'

model
 options
   missing  excludeall;
   step3 proportional ", method, ";
 output
      parameters=first betaopts=wl standarderrors=robust profile probmeans=posterior estimatedvalues=model
      write;
 variables
   caseweight ", weight,";
   dependent;
   latent Cluster nominal posterior = ( Cluster#1 Cluster#2 ) ;
 equations
end model"
  )))

  mLevelsToPaste = ifelse(mLevels == "continuous", "", mLevels)

  if(analysis == "dependent"){
    syntaxExploreToBe[grep("dependent", syntaxExploreToBe)] =
      utils::capture.output(cat(paste0("   dependent ", paste(Covariates, mLevels, collapse = ", "), ";", sep = "")))
    syntaxExploreToBe[length(syntaxExploreToBe)] = paste0("   ", Covariates[1]," <- 1 + Cluster;")
    if(length(Covariates)>1){
      for(idxCov in 2:length(Covariates)){
        syntaxExploreToBe[length(syntaxExploreToBe) + 1] = paste0("   ", Covariates[idxCov]," <- 1 + Cluster;")
      }}
    for(idxCovCont in which(mLevels == "continuous")){
      syntaxExploreToBe[length(syntaxExploreToBe) + 1] = paste0("   ", Covariates[idxCovCont],";")
    }}
  if(analysis == "covariates"){
    syntaxExploreToBe[grep("dependent", syntaxExploreToBe)] =
      utils::capture.output(cat(paste0("   independent ", paste(Covariates, mLevelsToPaste, collapse = ", "), ";", sep = "")))
    syntaxExploreToBe[length(syntaxExploreToBe)] = paste0("   Cluster <- 1 + ", paste(Covariates, collapse = "+"),";")
  }
  syntaxExploreToBe[length(syntaxExploreToBe) + 1] = "end model"
  return(syntaxExploreToBe)
}

# function to remove missings
rna = function(x){
  x[!is.na(x)]
}

run3step = function(dirPost,
                    syntaxExpl,
                    splitSizes,
                    pSplits,
                    dirTreeResults,
                    Covariates,
                    posCovVal,
                    analysis){

  EV = Wald = Profile = Parms = list()
  nCov = length(Covariates)
  # to run the three step
  for(reps in 1:length(dirPost)){

    syntaxExpl[grep("infile", syntaxExpl)] =
      utils::capture.output(cat(paste0("infile \'", dirTreeResults, "/", dirPost[reps], "'")))

    syntaxExpl[grep("write", syntaxExpl)] =  paste0(
      "write = 'Results",
      pSplits[reps],
      ".csv' writeestimatedvalues='ev",
      pSplits[reps],
      ".txt';")

    ncharDirPost = nchar(strsplit(dirPost[reps], ".", fixed = TRUE)[[1]][1])
    syntaxExplTemp = sub("Cluster#1 Cluster#2",
                         paste0("Cluster#",
                                1:splitSizes[reps],
                                collapse = " "),
                         syntaxExpl)

    syntaxDir = paste0("sThree", reps, ".lgs")
    write.table(syntaxExplTemp,
                syntaxDir,
                row.names = FALSE,
                quote = FALSE,
                col.names = FALSE)

    shell(paste(LG, syntaxDir, "/b"))

    ncolCSV = max(count.fields(paste0("results",
                                      pSplits[reps],
                                      ".csv"), sep = ","))

    Results = read.table(paste0("results",
                                pSplits[reps],
                                ".csv"),
                         col.names = paste0("V",seq_len(ncolCSV)),
                         header = FALSE, sep =",", fill = TRUE)



    rowParms = which(Results[,5]=="Parameters")
    rowWald = which(Results[,5]=="WaldStatistics")
    rowProfile = which(Results[,5]=="Profile")
    Wald[[reps]] = rna(Results[rowWald,-c(1:5)])
    ParmsTemp = rna(Results[rowParms,-c(1:5)])
    Profile[[reps]] = rna(Results[rowProfile,-c(1:5)])

    matParms = rbind(0, matrix(ParmsTemp, nrow = splitSizes[reps] - 1))

    colnames(matParms) = c("Intercept", Covariates)
    rownames(matParms) = paste0(names(splitSizes)[reps], 1:splitSizes[reps])
    Parms[[reps]] = matParms

    if(analysis == "dependent"){
    EV[[reps]] =  unique(read.table(paste0("ev",
                                           pSplits[reps],
                                           ".txt"),
                                    sep = ",", header = TRUE))
    }

    if(analysis == "covariates"){
      dataSplit = read.delim(paste0(dirTreeResults, "/", dirPost[reps]))
      meanCov = apply(dataSplit[,Covariates], 2, function(x){mean(x, na.rm = TRUE)})
      betaAllCov = betaAllCovMean = list()
      length(betaAllCov) = length(betaAllCovMean) = nCov
      names(betaAllCov) = names(betaAllCovMean) = Covariates

      for(idxCov in 1:nCov){
        if(is.null(posCovVal)){
          posCovVal = sort(unique(rna(dataSplit[,Covariates[idxCov]])))
        }
        betaCovVal = sapply(matParms[,Covariates[idxCov]], function(x){x * posCovVal})
        betaCovMean = matParms[,Covariates[idxCov]] * meanCov[idxCov]
        rownames(betaCovVal) = posCovVal
        betaAllCov[[idxCov]] = betaCovVal
        betaAllCovMean[[idxCov]] = betaCovMean
      }

      EVcov = list()
      length(EVcov) = nCov
      names(EVcov) = Covariates
      for(idxCov in 1:nCov){
        expCov = list()
        for(idxClass in 1:splitSizes[reps]){
          expCov[[idxClass]] =  exp(matParms[idxClass,1] +
                                      betaAllCov[[idxCov]][,idxClass] +
                                      sum(sapply(betaAllCovMean[-idxCov],
                                                 function(x){x[idxClass]})))
        }
        allExp = Reduce("+", expCov)
        EVcov[[idxCov]] = sapply(expCov, function(x){x/allExp})
      }
      EV[[reps]] = EVcov
    }
  }

  WaldTotal = Wald[!sapply(Wald, is.null)]
  ParmsTotal = Parms
  ProfileTotal = Profile[!sapply(Profile, is.null)]

  names(WaldTotal) =
    names(ParmsTotal) =
    names(ProfileTotal) = pSplits

  toReturn3Step = list(WaldTotal = WaldTotal,
                       EV = EV,
                       ParmsTotal = ParmsTotal,
                       ProfileTotal = ProfileTotal)
}

reformResults = function(splitSizes,
                         ParmsTotal,
                         ProfileTotal,
                         EV, WaldTotal,
                         sizeMlevels,
                         mLevels,
                         Covariates,
                         pSplits,
                         analysis){
  splitInfo = list()

  for(idxSplits in 1:length(splitSizes)){
    nclassTemp = splitSizes[idxSplits]
    ParmsTemp = ParmsTotal[[idxSplits]]
    ProfileTemp = ProfileTotal[[idxSplits]]
    EVTemp = EV[[idxSplits]]

    EV1split = Parms1split = list()
    classProbTemp = ProfileTemp[1:nclassTemp]
    classProbTempOrder = order(ProfileTemp[1:nclassTemp], decreasing = TRUE)
    classProb = classProbTemp[classProbTempOrder]

    if(analysis == "dependent"){
    for(j in 1:length(Covariates)){
      EV1split[[j]] = EVTemp[,grepl(Covariates[j], colnames(EVTemp))]

      if(mLevels[j] == "ordinal"){
        EV1split

        Parms1split[[j]] = ParmsTemp[1:((sizeMlevels[j] - 1) + nclassTemp  - 1)]
        ParmsTemp = ParmsTemp[-c(1:((sizeMlevels[j] - 1) * nclassTemp))]
      }

      if(mLevels[j] == "continuous"){
        Parms1split[[j]] = ParmsTemp[1:nclassTemp]
        ParmsTemp = ParmsTemp[-c(1:((sizeMlevels[j] - 1) * nclassTemp))]
      }
    }
    }

    if(analysis == "covariates"){
      EV1split = EV[[idxSplits]]
      Parms1split = ParmsTotal[[idxSplits]]
    }

    splitInfo[[idxSplits]] = list(ClassProb = classProb,
                                  EV = EV1split,
                                  Parameters = Parms1split,
                                  Wald = WaldTotal)
  }
  names(splitInfo) = pSplits

  toReturn = list(splitInfo = splitInfo,
                  varInfo = list(mLevels = mLevels,
                                 sizeMlevels = sizeMlevels))
  return(toReturn)
}
