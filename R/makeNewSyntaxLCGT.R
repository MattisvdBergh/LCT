makeNewSyntaxLCGT = function(dataDir,
                             weight,
                             independent,
                             dependent,
                             caseid,
                             mLevelsLGindependent,
                             mLevelsLGdependent,
                             sets = 16,
                             iterations = 50){

  newSyntaxToBe = utils::capture.output(cat(paste("
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
    utils::capture.output(cat(paste0("   dependent ", paste(dependent, mLevelsLGdependent, collapse = ", "), ";", sep = "")))
  newSyntaxToBe[grep("\\bindependent\\b", newSyntaxToBe)] =
    utils::capture.output(cat(paste0("   independent ", paste(independent, mLevelsLGindependent, collapse = ", "), ";", sep = "")))

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
