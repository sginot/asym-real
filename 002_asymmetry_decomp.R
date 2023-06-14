# Application of asymmetry decomposition script

# Required packages:
library(abind)
library(scales)
source("../repo/001_functions.R")
source("../repo/002_asym_components.R")
source("Rfunctions1.txt")

#-------------------------------------------------------------------------------
# Global asymmetry decomposition

reord_LM <- locate.reorder(shape = pA$mshape,
                           along = 2)

decomp_asym <- mv.asym(A = pA$rotated,
                       reorder.LM = reord_LM,
                       Nrep = 1,
                       along = 2,
                       indiv.fac = ID)

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

#-------------------------------------------------------------------------------
# Check for mandibles only and check for effect of TA

shpmirshp <- decomp_asym$Procrustes$rotated

plot(shpmirshp[,1,], 
     shpmirshp[,2,],
     asp = 1)

plot(shpmirshp[,2,], 
     shpmirshp[,3,], 
     asp = 1)

distall <- decomp_asym$Procrustes$intereucl.dist

# Total asymmetry (TA) = distance between a shape and its mirror, averaged 
# across replicates.

mat_dist <- as.matrix(distall)

TA <- rep(NA,
          dim(mat_dist)[1]/4) # 2 shapes * 2 mirrored shapes

for (i in 1:length(TA)) {
        
        TA[i] <- mean(c(mat_dist[i*2-1, 
                                 i*2-1 + dim(mat_dist)[1]/2],
                        mat_dist[i*2, 
                                 i*2 + dim(mat_dist)[1]/2]))
}


linear_mod2 <- lm(bf ~ TA)

summary(linear_mod2)

TA_sq <- TA^2

quadratic_mod2 <- lm(bf ~ TA + TA_sq)

summary(quadratic_mod2)


# Mandibular LMs are 16:33 (because landmarks 8:10 were removed previously)

asym_mandi <- mv.asym(A = pA$rotated[16:33,,],
                       reorder.LM = reord_LM[16:33]-15,
                       Nrep = 1,
                       along = 2,
                       indiv.fac = ID)

plot(asym_mandi$Procrustes$rotated[,2,], 
     asym_mandi$Procrustes$rotated[,3,], 
     asp = 1)

plot(asym_mandi$PCA.sym$x[, 1:2])
plot(asym_mandi$PCA.asym$x[, 1:2])

DA_mandi <- asym_mandi$M[, 3]
FA_mandi <- asym_mandi$M[, 1]

linear_mod3 <- lm(bf ~ DA_mandi)

summary(linear_mod3)

DA_mandi_sq <- DA_mandi^2

quadratic_mod3 <- lm(bf ~ DA_mandi + DA_mandi_sq)

summary(quadratic_mod3)

# Non-mandibular LMs are -c(16:33)

asym_head <- mv.asym(A = pA$rotated[-c(16:33),,],
                      reorder.LM = c(reord_LM[1:15], reord_LM[34:35]-18),
                      Nrep = 1,
                      along = 2,
                      indiv.fac = ID)

plot(asym_head$Procrustes$rotated[,2,], 
     asym_head$Procrustes$rotated[,3,], 
     asp = 1)

plot(asym_head$PCA.sym$x[, 1:2])
plot(asym_head$PCA.asym$x[, 1:2])

DA_head <- asym_head$M[, 3]
FA_head <- asym_head$M[, 1]


linear_mod4 <- lm(bf ~ DA_head)

summary(linear_mod4)

DA_head_sq <- DA_head^2

quadratic_mod4 <- lm(bf ~ DA_head + DA_head_sq)

summary(quadratic_mod4)

#-------------------------------------------------------------------------------
# (M)ANOVAs for significance of global / partial / individual coordinate asym

summary(manova(lm(prcomp(asym_head$matshp)$x[,1:50] ~ 
                          asym_head$indiv.fac.mirror * asym_head$mirror.fac)))
summary(manova(lm(prcomp(asym_mandi$matshp)$x[,1:50] ~ 
                          asym_mandi$indiv.fac.mirror * asym_mandi$mirror.fac)))
summary(manova(lm(prcomp(decomp_asym$matshp)$x[,1:50] ~ 
                          decomp_asym$indiv.fac.mirror * decomp_asym$mirror.fac)))

summary(aov(lm(decomp_asym$matshp ~ 
                       asym_mandi$indiv.fac.mirror * asym_head$mirror.fac)))

# Make visual representation

mshp <- decomp_asym$Procrustes$mshape


# For INTER INDIVIDUAL VARIATION

MSx <- MSy <- MSz <- rep(NA, 35)


for (i in 1:35) {
   
smx <- summary(aov(lm(decomp_asym$matshp[, i] ~ 
                       asym_mandi$indiv.fac.mirror * asym_head$mirror.fac)))

smy <- summary(aov(lm(decomp_asym$matshp[, i+35] ~ 
                              asym_mandi$indiv.fac.mirror * asym_head$mirror.fac)))

smz <- summary(aov(lm(decomp_asym$matshp[, i+70] ~ 
                              asym_mandi$indiv.fac.mirror * asym_head$mirror.fac)))

if (smx[[1]]$`Pr(>F)`[1] < 0.05) {
        MSx[i] <- smx[[1]]$`Mean Sq`[1]
} else {MSx[i] <- 0}
        
if (smy[[1]]$`Pr(>F)`[1] < 0.05) {
        MSy[i] <- smy[[1]]$`Mean Sq`[1]
} else {MSy[i] <- 0}

if (smz[[1]]$`Pr(>F)`[1] < 0.05) {
        MSz[i] <- smz[[1]]$`Mean Sq`[1]
} else {MSz[i] <- 0}

}


plot(decomp_asym$Procrustes$mshape[,1],
     decomp_asym$Procrustes$mshape[,2],
     asp = 1,
     pch = 21,
     cex = 2,
     bg = alpha("firebrick",
                alpha =  MSx / max(MSx)),
     lwd = 2,
     xlab = "",
     ylab = "Right <=> Left")

mtext("Posterior <=> Anterior", 
      side = 1, 
      col = "firebrick", 
      line = 3,
      font = 2)

plot(decomp_asym$Procrustes$mshape[,2],
     decomp_asym$Procrustes$mshape[,3],
     asp = 1,
     pch = 21,
     cex = 2,
     bg = alpha("forestgreen",
                alpha =  MSy / max(MSy)),
     lwd = 2,
     ylab = "Ventral <=> Dorsal",
     xlab = "")

mtext("Right <=> Left", 
      side = 1, 
      col = "forestgreen", 
      line = 3,
      font = 2)

plot(decomp_asym$Procrustes$mshape[,3],
     decomp_asym$Procrustes$mshape[,1],
     asp = 1,
     pch = 21,
     cex = 2,
     bg = alpha("blue",
                alpha =  MSz / max(MSz)),
     lwd = 2,
     xlab = "",
     ylab = "Posterior <=> Anterior")

mtext("Ventral <=> Dorsal", 
      side = 1, 
      col = "blue", 
      line = 3,
      font = 2)

# For DIRECTIONAL ASYMMETRY

MSx <- MSy <- MSz <- rep(NA, 35)


for (i in 1:35) {
        
        smx <- summary(aov(lm(decomp_asym$matshp[, i] ~ 
                                      asym_mandi$indiv.fac.mirror * asym_head$mirror.fac)))
        
        smy <- summary(aov(lm(decomp_asym$matshp[, i+35] ~ 
                                      asym_mandi$indiv.fac.mirror * asym_head$mirror.fac)))
        
        smz <- summary(aov(lm(decomp_asym$matshp[, i+70] ~ 
                                      asym_mandi$indiv.fac.mirror * asym_head$mirror.fac)))
        
        if (smx[[1]]$`Pr(>F)`[2] < 0.05) {
                MSx[i] <- smx[[1]]$`Mean Sq`[2]
        } else {MSx[i] <- 0}
        
        if (smy[[1]]$`Pr(>F)`[2] < 0.05) {
                MSy[i] <- smy[[1]]$`Mean Sq`[2]
        } else {MSy[i] <- 0}
        
        if (smz[[1]]$`Pr(>F)`[2] < 0.05) {
                MSz[i] <- smz[[1]]$`Mean Sq`[2]
        } else {MSz[i] <- 0}
        
}


plot(decomp_asym$Procrustes$mshape[,1],
     decomp_asym$Procrustes$mshape[,2],
     asp = 1,
     pch = 21,
     cex = 2,
     bg = alpha("firebrick",
                alpha =  MSx / max(MSx)),
     lwd = 2,
     xlab = "",
     ylab = "Right <=> Left")

mtext("Posterior <=> Anterior", 
      side = 1, 
      col = "firebrick", 
      line = 3,
      font = 2)

plot(decomp_asym$Procrustes$mshape[,2],
     decomp_asym$Procrustes$mshape[,3],
     asp = 1,
     pch = 21,
     cex = 2,
     bg = alpha("forestgreen",
                alpha =  MSy / max(MSy)),
     lwd = 2,
     ylab = "Ventral <=> Dorsal",
     xlab = "")

mtext("Right <=> Left", 
      side = 1, 
      col = "forestgreen", 
      line = 3,
      font = 2)

plot(decomp_asym$Procrustes$mshape[,3],
     decomp_asym$Procrustes$mshape[,1],
     asp = 1,
     pch = 21,
     cex = 2,
     bg = alpha("blue",
                alpha =  MSz / max(MSz)),
     lwd = 2,
     xlab = "",
     ylab = "Posterior <=> Anterior")

mtext("Ventral <=> Dorsal", 
      side = 1, 
      col = "blue", 
      line = 3,
      font = 2)

# For FLUCTUATING ASYMMETRY

MSx <- MSy <- MSz <- rep(NA, 35)


for (i in 1:35) {
        
        smx <- summary(aov(lm(decomp_asym$matshp[, i] ~ 
                                      asym_mandi$indiv.fac.mirror * asym_head$mirror.fac)))
        
        smy <- summary(aov(lm(decomp_asym$matshp[, i+35] ~ 
                                      asym_mandi$indiv.fac.mirror * asym_head$mirror.fac)))
        
        smz <- summary(aov(lm(decomp_asym$matshp[, i+70] ~ 
                                      asym_mandi$indiv.fac.mirror * asym_head$mirror.fac)))
        
        if (smx[[1]]$`Pr(>F)`[3] < 0.05) {
                MSx[i] <- smx[[1]]$`Mean Sq`[3]
        } else {MSx[i] <- 0}
        
        if (smy[[1]]$`Pr(>F)`[3] < 0.05) {
                MSy[i] <- smy[[1]]$`Mean Sq`[3]
        } else {MSy[i] <- 0}
        
        if (smz[[1]]$`Pr(>F)`[3] < 0.05) {
                MSz[i] <- smz[[1]]$`Mean Sq`[3]
        } else {MSz[i] <- 0}
        
}


plot(decomp_asym$Procrustes$mshape[,1],
     decomp_asym$Procrustes$mshape[,2],
     asp = 1,
     pch = 21,
     cex = 2,
     bg = alpha("firebrick",
                alpha =  MSx / max(MSx)),
     lwd = 2,
     xlab = "",
     ylab = "Right <=> Left")

mtext("Posterior <=> Anterior", 
      side = 1, 
      col = "firebrick", 
      line = 3,
      font = 2)

plot(decomp_asym$Procrustes$mshape[,2],
     decomp_asym$Procrustes$mshape[,3],
     asp = 1,
     pch = 21,
     cex = 2,
     bg = alpha("forestgreen",
                alpha =  MSy / max(MSy)),
     lwd = 2,
     ylab = "Ventral <=> Dorsal",
     xlab = "")

mtext("Right <=> Left", 
      side = 1, 
      col = "forestgreen", 
      line = 3,
      font = 2)

plot(decomp_asym$Procrustes$mshape[,3],
     decomp_asym$Procrustes$mshape[,1],
     asp = 1,
     pch = 21,
     cex = 2,
     bg = alpha("blue",
                alpha =  MSz / max(MSz)),
     lwd = 2,
     xlab = "",
     ylab = "Posterior <=> Anterior")

mtext("Ventral <=> Dorsal", 
      side = 1, 
      col = "blue", 
      line = 3,
      font = 2)

