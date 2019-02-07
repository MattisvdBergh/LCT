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
                        sep = sep,
                        ...){
  
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