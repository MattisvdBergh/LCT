dataDir = "C:/Users/Huub/Desktop/TryOut/data.dat"
LG = "C:/Users/Huub/Latent_Gold_5.1/LatentGOLD5.1/lg51.exe"

library(poLCA)

N = 10000
cp1 = matrix(c(seq(0.1, 1, 0.1),
              1 - seq(0.1, 1, 0.1)),
            ncol = 2)

cp2 = matrix(c(seq(1, 0.1, -0.1),
              1 - seq(1, 0.1, -0.1)),
            ncol = 2)

cp3 = matrix(c(rep(0.25, 5), 1 - rep(0.25, 5),
               1 - rep(0.25, 5), rep(0.25, 5)),
             ncol = 2)

probs1 = do.call(list, replicate(10, cp1, simplify = FALSE))
probs2 = do.call(list, replicate(10, cp2, simplify = FALSE))
probs3 = do.call(list, replicate(10, cp3, simplify = FALSE))

simdat1 = poLCA.simdata(N=N,c(probs1))
simdat2 = poLCA.simdata(N=N,c(probs1, probs2))
simdat3 = poLCA.simdata(N=N,c(probs1, probs3))

probsEx <- list(matrix(c(0.6,0.1,0.3,     0.6,0.3,0.1,     0.3,0.1,0.6    ),ncol=3,byrow=TRUE), # Y1
              matrix(c(0.2,0.8,         0.7,0.3,         0.3,0.7        ),ncol=2,byrow=TRUE), # Y2
              matrix(c(0.3,0.6,0.1,     0.1,0.3,0.6,     0.3,0.6,0.1    ),ncol=3,byrow=TRUE), # Y3
              matrix(c(0.1,0.1,0.5,0.3, 0.5,0.3,0.1,0.1, 0.3,0.1,0.1,0.5),ncol=4,byrow=TRUE), # Y4
              matrix(c(0.1,0.1,0.8,     0.1,0.8,0.1,     0.8,0.1,0.1    ),ncol=3,byrow=TRUE)) # Y5
simdatEx <- poLCA.simdata(N=1000, probsEx, P=c(0.2,0.3,0.5))

simOrd = as.data.frame(apply(simdatEx$dat, 2, as.factor))
To2 = LCT(dataDir, LG = LG, resultsName = "SimOrd2", maxClassSplit1 = 2)

### cont
simCont = as.data.frame(matrix(rnorm(10000, 1:20, 30:10), 1000, 10))
colnames(simCont) = paste0("y", 1:10)
Tc2 = LCT(simCont, LG = LG, resultsName = "SimCont2", maxClassSplit1 = 2)

### ord + cont
simOrdCont = cbind(as.data.frame(apply(simdat1$dat, 2, as.character)), cont = simCont[,1])
Toc2 = LCT(simCont, LG = LG, resultsName = "SimOrdCont2", maxClassSplit1 = 2)





