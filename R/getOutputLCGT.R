getOutputLCGT = function(resultsTemp, Hclass, maxClassSplit1,
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
