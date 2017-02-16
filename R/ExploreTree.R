#' Explore a Latent Class Tree model
#'
#' Explore a Latent Class Tree model with Latent GOLD 5.1
#'
#' @param res Treeobject from
#' @param LG Path to the Latent GOLD executable
#' @param levelsCovariates Vector with the levels of each of the covariates (1 = continuous, 1< = ordinal)
#' @param measurementLevels Whether the ordering of classes should be decreasing or not. Defaults to TRUE.
#' @param dirTreeResults Maximum size of the first split of the tree. Will be assessed with the criterion given in stopCriterium. Defaults to two.
#' @param ResultsFolder folder with results
#'
#' @return None
#' @export
exploreTree = function(res,
                       syntaxThree = NULL,
                       mLevelsCovariates = c(1, 2),
                       method = "bch",
                       dirTreeResults = getwd(),
                       ResultsFolder = "exploreTree",
                       Covariates = c("age", "sex"),
                       mLevelsLGdependent = c("continuous", "ordinal"),
                       weight = "weight"){
  mainDir = getwd()

  # unique splits
  uSplits = res$Splitpoints

  # All files with posteriors
  LdirPost = list()
  allSplits = numeric()
  for(i in 1:length(res$Splits)){
    allSplits = c(allSplits, res$Splits[[i]])
    LdirPost[[i]] = paste0("H", i,
                           "c", res$Names[[i]],
                           "_sol", res$Splits[[i]], ".txt")
  }

  # To determine which files of posteriors really belong to a split
  LrealSplits = sapply(LdirPost, function(x){
    nCharLevel = nchar(x)
    substr(x, nCharLevel - 4, nCharLevel - 4) != 1
  })

  # Make a correction for all real splits
  realSplitsLogical = unlist(LrealSplits)
  dirPost = unlist(LdirPost)[realSplitsLogical]
  realSplits = allSplits[realSplitsLogical]

  # function to remove missings
  rna = function(x){
    x[!is.na(x)]
  }

  dir.create(file.path(getwd(), ResultsFolder), showWarnings = FALSE)
  setwd(file.path(mainDir, ResultsFolder))

  if(is.null(syntaxThree)){
    syntaxThree = makeNewSyntaxExplore(dataDir = dirPost[1],
                                       mLevelsLGdependent = mLevelsLGdependent,
                                       Covariates = Covariates)
  }

  nClass = numeric()
  EV = Wald = Profile = Parms = list()

  # to run the three step
  for(reps in 1:length(dirPost)){

    syntaxThree[grep("infile", syntaxThree)] =
      capture.output(cat(paste0("infile \'", dirTreeResults, "/", dirPost[reps], "'")))

    syntaxThree[grep("write", syntaxThree)] =  paste0(
      "write = 'Results", uSplits[reps], ".csv' writeestimatedvalues='ev", uSplits[reps], ".txt';")


    ncharDirPost = nchar(strsplit(dirPost[reps], ".", fixed = TRUE)[[1]][1])
    nClass[reps] = substr(dirPost[reps], ncharDirPost, ncharDirPost)
    syntaxThreeTemp = sub("Cluster#1 Cluster#2", paste0("Cluster#", 1:nClass[reps], collapse = " "), syntaxThree)
    syntaxDir = paste0("sThree", reps, ".lgs")
    write.table(syntaxThreeTemp, syntaxDir, row.names = FALSE, quote = FALSE, col.names = FALSE)
    shell(paste(LG, syntaxDir, "/b"))

    ncolCSV = max(count.fields(paste0(
      "results",
      res$Splitpoints[reps],
      ".csv"), sep = ","))

    Results = read.table(paste0(
      "results",
      res$Splitpoints[reps],
      ".csv"),
      header = FALSE, col.names = paste0("V",seq_len(ncolCSV)), sep =",", fill = TRUE)

    EV[[reps]] =  unique(read.table(paste0(
      "ev",
      res$Splitpoints[reps],
      ".txt"), sep = ",", header = TRUE))

    rowParms = which(Results[,5]=="Parameters")
    rowWald = which(Results[,5]=="WaldStatistics")
    rowProfile = which(Results[,5]=="Profile")
    Wald[[reps]] = Results[rowWald,-c(1:5)]
    Wald[[reps]] = rna(Wald[[reps]])
    Parms[[reps]] = Results[rowParms,-c(1:5)]
    Parms[[reps]] = rna(Parms[[reps]])
    Profile[[reps]] = Results[rowProfile,-c(1:5)]
    Profile[[reps]] = rna(Profile[[reps]])
  }

  nClassTotal = as.numeric(nClass[!is.na(nClass)])
  WaldTotal = Wald[!sapply(Wald, is.null)]
  ParmsTotal = Parms[!sapply(Parms, is.null)]
  ProfileTotal = Profile[!sapply(Profile, is.null)]

  names(nClassTotal) =
    names(WaldTotal) =
    names(ParmsTotal) =
    names(ProfileTotal) = unlist(res$Names[-length(res$Names)])[realSplitsLogical]

  toReturn = reformResults(nClassTotal = nClassTotal,
                           ParmsTotal = ParmsTotal,
                           ProfileTotal = ProfileTotal,
                           EV = EV,
                           WaldTotal = WaldTotal,
                           mLevelsCovariates = mLevelsCovariates,
                           mLevelsLGdependent = mLevelsLGdependent,
                           Covariates = Covariates)
  setwd(mainDir)
  return(toReturn)
}

makeNewSyntaxExplore = function(dataDir = dataDir,
                                method = "bch",
                                mLevelsLGdependent = mLevelsLGdependent,
                                Covariates = Covariates,
                                weight = "weight"){

  syntaxExploreToBe = capture.output(cat(paste("
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

  syntaxExploreToBe[grep("dependent", syntaxExploreToBe)] =
    capture.output(cat(paste0("   dependent ", paste(Covariates, mLevelsLGdependent, collapse = ", "), ";", sep = "")))
  syntaxExploreToBe[length(syntaxExploreToBe)] = paste0("   ", Covariates[1]," <- 1 + Cluster;")
  if(length(Covariates)>1){
  for(idxCov in 2:length(Covariates)){
    syntaxExploreToBe[length(syntaxExploreToBe) + 1] = paste0("   ", Covariates[idxCov]," <- 1 + Cluster;")
  }}
  for(idxCovCont in which(mLevelsLGdependent == "continuous")){
    syntaxExploreToBe[length(syntaxExploreToBe) + 1] = paste0("   ", Covariates[idxCovCont],";")
  }
  syntaxExploreToBe[length(syntaxExploreToBe) + 1] = "end model"
  return(syntaxExploreToBe)
}

reformResults = function(nClassTotal, ParmsTotal, ProfileTotal,
                         EV, WaldTotal, mLevelsCovariates,
                         mLevelsLGdependent, Covariates){
  toReturn = list()

  for(idxSplits in 1:length(nClassTotal)){
    nclassTemp = nClassTotal[idxSplits]
    ParmsTemp = ParmsTotal[[idxSplits]]
    ProfileTemp = ProfileTotal[[idxSplits]]
    EVTemp = EV[[idxSplits]]

    EV1split = Parms1split = list()
    classProbTemp = ProfileTemp[1:nclassTemp]
    classProbTempOrder = order(ProfileTemp[1:nclassTemp], decreasing = TRUE)
    classProb = classProbTemp[classProbTempOrder]


    for(j in 1:length(Covariates)){
      EV1split[[j]] = EVTemp[,grepl(Covariates[j], colnames(EVTemp))]


      if(mLevelsLGdependent[j] == "nominal"){
        # EV1split[[j]] = matrix(
          # EVTemp[1:(mLevelsCovariates[j] * nclassTemp)],
          # ncol = nclassTemp)
        Parms1split[[j]] = ParmsTemp[1:((mLevelsCovariates[j] - 1) * nclassTemp)]
        ParmsTemp = ParmsTemp[-c(1:((mLevelsCovariates[j] - 1) * nclassTemp))]
        # EVTemp = EVTemp[-c(1:(mLevelsCovariates[j] * nclassTemp))]
      }

      if(mLevelsLGdependent[j] == "ordinal"){
        Parms1split[[j]] = ParmsTemp[1:((mLevelsCovariates[j] - 1) + nclassTemp  - 1)]
        ParmsTemp = ParmsTemp[-c(1:((mLevelsCovariates[j] - 1) * nclassTemp))]
      }

      if(mLevelsLGdependent[j] == "continuous"){
        # EV1split[[j]] = EVTemp[1:nclassTemp]
        Parms1split[[j]] = ParmsTemp[1:nclassTemp]
        ParmsTemp = ParmsTemp[-c(1:((mLevelsCovariates[j] - 1) * nclassTemp))]
        # EVTemp = EVTemp[1:nclassTemp]
      }
    }
    toReturn[[idxSplits]] = list(ClassProb = classProb,
                                 EV = EV1split,
                                 Parameters = Parms1split,
                                 Wald = WaldTotal[[idxSplits]],
                                 MeasurementLevels = mLevelsLGdependent)
  }
  names(toReturn) = names(nClassTotal)
  return(toReturn)
}
