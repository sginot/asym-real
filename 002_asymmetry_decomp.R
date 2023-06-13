# Application of asymmetry decomposition script

# Required packages:
library(abind)
source("../repo/001_functions.R")
source("../repo/002_asym_components.R")
source("Rfunctions1.txt")
#-------------------------------------------------------------------------------

reord_LM <- locate.reorder(shape = pA$mshape,
                           along = 2)

decomp_asym <- mv.asym(A = pA$rotated,
                       reorder.LM = reord_LM,
                       Nrep = 1)

plot(decomp_asym$PCA.sym$x[,1:2], 
     asp = 1)

plot(decomp_asym$PCA.asym$x[,1:2], 
     xlim = c(-max(abs(decomp_asym$PCA.asym$x[, 1])), 
              0),
     ylim = c(-max(abs(decomp_asym$PCA.asym$x[, 2])), 
              max(abs(decomp_asym$PCA.asym$x[, 2]))))
abline(h = 0,
       v = 0,
       col = "gray")


DA <- decomp_asym$M[, 3]
FA <- decomp_asym$M[, 1]

cor.test(csiz.av, FA)
cor.test(bf, FA)

cor.test(csiz.av, DA)
cor.test(bf, DA)

linear_mod <- lm(bf ~ DA)

DA_sq <- DA^2

quadratic_mod <- lm(bf ~ DA + DA_sq)

summary(quadratic_mod)

newvals <- seq(-0.2, 
               0,
               by = 0.001)

pred_quadra <- predict(quadratic_mod, newdata = data.frame(DA = newvals,
                                                           DA_sq = newvals^2))

plot(DA, bf)
lines(newvals, pred_quadra)
