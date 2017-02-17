makeNewSyntax = function(dataDir, itemNames, weight, mLevels,
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
    capture.output(cat(paste0("   dependent ", paste(itemNames, mLevels, collapse = ", "), ";", sep = "")))

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
