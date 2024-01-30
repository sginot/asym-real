# Re-analysis of the data following reviews
# Main change will be the separate analysis of mandibles and head
# For the asymmetry part this means switching from object symmetry to matching
# symmetry. For the modularity and integration part this means mostly focusing 
# covariation patterns.

# Run scripts 001 to 003 before this one.

#-------------------------------------------------------------------------------
# Required packages
source(file = "Rfunctions1.txt") # Functions from Claude 2008
library(scales)
library(car)
library(rgl)
library(geomorph)
library(EMMLi)
library(paleomorph)

#-------------------------------------------------------------------------------
# Raw data organisation and superimposition

LM_template[-c(8:10), 1] # Shows the list of landmarks, to identify the ones 
                        # belonging to each mandible
dim(arr) # Array containing raw coordinates, with 3 noisy LMs removed.

LM_LM <- c(17, 21:23, 25, 30:33) # Index of landmarks for left mandible
LM_RM <- c(16, 18:20, 24, 26:29) # Index of LMs for right mandible

A_head <- arr[-c(LM_LM, LM_RM),,] # Array containing only head LMs

A_LM <- arr[c(LM_LM),,] # Array containing only left mandible landmarks
A_LM[, 2, ] <- -1 * A_LM[, 2, ] # Mirror configurations along medio-lateral axis

A_RM <- arr[c(LM_RM),,] # Array containing only right mandible landmarks

A_mand <- array(data = NA, #Empty array to gather mandible shape data
                dim = c(dim(A_LM)[1:2], 
                        dim(A_LM)[3] * 2))

for (i in 1:dim(A_LM)[3]) { 
  A_mand[,, i * 2] <- A_LM[,, i]
  A_mand[,, i * 2 - 1] <- A_RM[,, i]
} 
# Filling empty array, alternating left and right shapes. Consecutive elements
# of the array belong to the same individual.

pA_mand <- pgpa(A_mand) # Partial Procrustes sumperimposition

shp_mand <- orp(pA_mand$rotated) # Orthogonal projection into Euclidean space

side <- as.factor(rep(c("R", "L"), 
                      length(ID)))

mshp_RM <- mshape(shp_mand[,, which(side == "R")])
# Mean shape of shapes with odd index -> Right mandible

mshp_LM <- mshape(shp_mand[,, which(side == "L")])
# Mean shape of shapes with even index -> Left mandible

pA_head <- pgpa(A_head) # Partial Procrustes sumperimposition
shp_head <- orp(pA_head$rotated) # Orthogonal projection into Euclidean space
mshp_head <- pA_head$mshape # Mean shape of head LMs

ID_head <- ID # Recycle vector from script 001
ID_mand <- rep(ID, each = 2) # Same vector, with duplicated elements because
                            # each individual has two mandibles.

replic_head <- replic
replic_mand <- rep(replic, each = 2)

#-------------------------------------------------------------------------------
# Asymmetry analyses

# For mandibles, use matching symmetry
bilatsym_mand <- bilat.symmetry(A = shp_mand,
                                ind = ID_mand,
                                side = side,
                                replicate = replic_mand,
                                object.sym = F)

# For head, must define pairs of landmarks to use object symmetry analysis
LM_pairs <- list()

# Plot and text to assess visually the pairs on LMs
plot(shp_head[,2,], 
     shp_head[,3,], 
     asp = 1)

text(shp_head[,2,1], 
     shp_head[,3,1], 
     col = 2, 
     pch = 19, 
     cex = 2)

LM_pairs[[1]] <- c(3, 12)
LM_pairs[[2]] <- c(4, 11)
LM_pairs[[3]] <- c(5, 9)
LM_pairs[[4]] <- c(6, 8)
LM_pairs[[5]] <- c(7, 10)
LM_pairs[[6]] <- c(14, 15)
LM_pairs[[7]] <- c(16, 17)

LM_pairs <- matrix(unlist(LM_pairs), 
                     ncol = 2, 
                     byrow = T)

bilatsym_head <- bilat.symmetry(A = shp_head,
                                ind = ID_head,
                                replicate = replic_head,
                                object.sym = T,
                                land.pairs = LM_pairs)

bilatsym_head2D <- bilat.symmetry(A = shp_head[,2:3,],
                                ind = ID_head,
                                replicate = replic_head,
                                object.sym = T,
                                land.pairs = LM_pairs)

plot(bilatsym_head2D)

#-------------------------------------------------------------------------------
# Show covariance patterns via plots, separated between overall variation,
# symmetric component and asymmetric component.

# Make matrices averaging replicates

# individual:side for the mandible
mat_shp_mand <- matrix(NA, 
                      nrow = dim(shp_mand)[3], 
                      ncol = dim(shp_mand)[1] * dim(shp_mand)[2])

for (i in 1:dim(mat_shp_mand)[1]) {
  mat_shp_mand[i,] <- c(t(shp_mand[,,i]))
}

mat_shp_mand_av <- apply(X = mat_shp_mand, 
                         MARGIN = 2, 
                         FUN = tapply, 
                         INDEX = ID_mand:side, 
                         mean)

# individual average for the head
mat_shp_head <- matrix(NA, 
                       nrow = dim(shp_head)[3], 
                       ncol = dim(shp_head)[1] * dim(shp_head)[2])

for (i in 1:dim(mat_shp_head)[1]) {
  mat_shp_head[i,] <- c(t(shp_head[,,i]))
}

mat_shp_head_av <- apply(X = mat_shp_head, 
                         MARGIN = 2, 
                         FUN = tapply, 
                         INDEX = ID_head, 
                         mean)

# Make global coordinate matrix, combining left and right mand with head for
# each individual

mat_overall <- matrix(NA, 
                      nrow = dim(mat_shp_head_av)[1], 
                      ncol = dim(mat_shp_head_av)[2] + 
                        dim(mat_shp_mand_av)[2] * 2)
  
for (i in 1:dim(mat_overall)[1]) {
  mat_overall[i, 1:51] <- mat_shp_head_av[i, ]
  mat_overall[i, 52:78] <- mat_shp_mand_av[grep(pattern = ":L", 
                                                x = rownames(mat_shp_mand_av))[i],]
  mat_overall[i, 79:105] <- mat_shp_mand_av[grep(pattern = ":R", 
                                                 x = rownames(mat_shp_mand_av))[i],]
}

rownames(mat_overall) <- rownames(mat_shp_head_av)

# Compute covariance and correlation matrices

color <- hcl.colors(30, 
           palette = "viridis",
           alpha = NULL, 
           rev = FALSE, 
           fixup = TRUE)

cov_overall <- abs(cov(mat_overall))
cor_overall <- abs(cor(mat_overall))

image(x = 1:dim(cov_overall)[1], 
      y = 1:dim(cov_overall)[1], 
      z = cov_overall, 
      col = color,
      main = "Landmarks covariance matrix",
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n")

abline(h = c(51.5, 78.5),
       v = c(51.5, 78.5),
       col = "grey",
       lwd = 3)

image(x = 1:dim(cor_overall)[1], 
      y = 1:dim(cor_overall)[1], 
      z = cor_overall, 
      col = color,
      main = "Landmarks correlation matrix",
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n")

abline(h = c(51.5, 78.5),
       v = c(51.5, 78.5),
       col = "grey",
       lwd = 3)

# Compute congruence index correlation matrix
A_overall <- array(NA,
                   dim = c(dim(mat_overall)[2] / 3, 
                           3, 
                           dim(mat_overall)[1]))

for (i in 1:dim(A_overall)[3]) {
  A_overall[,,i] <- matrix(mat_overall[i,], 
                           ncol = 3, 
                           byrow = T)
}

congru_overall <- abs(dotcorr(A_overall))

image(x = 1:dim(congru_overall)[1], 
      y = 1:dim(congru_overall)[1], 
      z = congru_overall, 
      col = color,
      main = "Landmarks congruence matrix",
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n")

abline(v = c(17.5, 26.5),
       h = c(17.5, 26.5),
       lwd = 3,
       col = "grey")

label_head <- DF[-c(LM_LM, LM_RM), 10]
label_LM <- DF[LM_LM, 10]
label_RM <- DF[LM_RM, 10]

lab_congru <- c(label_head, label_LM, label_RM)

axis(side = 1, 
     at = 1:dim(congru_overall)[2], 
     labels = lab_congru, 
     las = 2)

axis(side = 2, 
     at = 1:dim(congru_overall)[2], 
     labels = lab_congru, 
     las = 2)

#-------------------------------------------------------------------------------
# Check if asymmetric components from head and mandi are correlated

A_asym <- array(NA,
                dim = c(dim(bilatsym_head$asymm.shape)[1] + 
                          dim(bilatsym_mand$asymm.shape)[1],
                        dim(bilatsym_head$asymm.shape)[2:3]))

for (i in 1:dim(A_asym)[3]) {
  A_asym[1:17,,i] <- bilatsym_head$asymm.shape[,,i]
  A_asym[18:26,,i] <- bilatsym_mand$asymm.shape[,,i]
}

congru_asym <- abs(dotcorr(A_asym))

image(x = 1:dim(congru_asym)[1], 
      y = 1:dim(congru_asym)[1], 
      z = congru_asym, 
      col = color,
      main = "Landmarks congruence matrix",
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n")

abline(v = 17.5,
       h = 17.5,
       lwd = 3,
       col = "grey")

#-------------------------------------------------------------------------------
# Same for symmetric components from head and mandi

A_sym <- array(NA,
                dim = c(dim(bilatsym_head$symm.shape)[1] + 
                          dim(bilatsym_mand$symm.shape)[1],
                        dim(bilatsym_head$symm.shape)[2:3]))

for (i in 1:dim(A_sym)[3]) {
  A_sym[1:17,,i] <- bilatsym_head$symm.shape[,,i]
  A_sym[18:26,,i] <- bilatsym_mand$symm.shape[,,i]
}

congru_sym <- abs(dotcorr(A_sym))

image(x = 1:dim(congru_sym)[1], 
      y = 1:dim(congru_sym)[1], 
      z = congru_sym, 
      col = color,
      main = "Landmarks congruence matrix",
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n")

abline(v = 17.5,
       h = 17.5,
       lwd = 3,
       col = "grey")

#-------------------------------------------------------------------------------
# Test for 2BPLS relation between shape components and bf

bf_sym_head_pls <- two.b.pls(A1 = bilatsym_head$symm.shape[,,-which(is.na(bf2))], 
                        A2 = bf2[-which(is.na(bf2))])
bf_asym_head_pls <- two.b.pls(A1 = bilatsym_head$asymm.shape[,,-which(is.na(bf2))], 
                        A2 = bf2[-which(is.na(bf2))])

bf_sym_mand_pls <- two.b.pls(A1 = bilatsym_mand$symm.shape[,,-which(is.na(bf2))], 
                        A2 = bf2[-which(is.na(bf2))])
bf_asym_mand_pls <- two.b.pls(A1 = bilatsym_mand$asymm.shape[,,-which(is.na(bf2))], 
                             A2 = bf2[-which(is.na(bf2))])

#-------------------------------------------------------------------------------
# Test for modularity within the head only
# Then test for integration between each mandible and with head modules if signi

modu_head <- DF[which(DF[,6] == 1 | DF[,6] ==  4),
                c(6, 10)]

mod_test_head <- modularity.test(A = shp_head,
                                 partition.gp = modu_head[,1],
                                 iter = 999, 
                                 CI = T)

integ_test_head <- integration.test(A = shp_head,
                                    partition.gp = modu_head[,1],
                                    iter = 999)

plot(A_overall[,2:3,1], asp= 1)

part_overall <- as.factor(c(rep("head", 17), 
                            rep("LM", 9), 
                            rep("RM", 9)))

part_overall2 <- as.factor(c(modu_head[,1], 
                            rep("LM", 9), 
                            rep("RM", 9)))

levels(part_overall2)[1:2] <- c("head", "sensory")

integ_test_3mod <- integration.test(A = A_overall,
                                    partition.gp = part_overall,
                                    iter = 999)

integ_test_4mod <- integration.test(A = A_overall,
                                    partition.gp = part_overall2,
                                    iter = 999)
#-------------------------------------------------------------------------------
# New matrix plot with reordered LMs, and integration test results

modo_head <- modu_head[order(modu_head[,1]),]

lab_congro <- lab_congru
lab_congro[1:17] <- lab_congru[order(modu_head[,1])]

congro_overall <- congru_overall
congro_overall[1:17, 1:17] <- congro_overall[order(modu_head[,1]),
                                             order(modu_head[,1])]

m <- congro_overall

m[upper.tri(congro_overall, diag = T)] <- NA


rpls_val <- paste("r-PLS =", 
                  round(integ_test_4mod$r.pls.mat, 3))

p_val <- paste("P =", integ_test_4mod$pairwise.P.values)

parts <- levels(part_overall2)

image(x = 1:dim(m)[1], 
      y = 1:dim(m)[1], 
      z = m, 
      col = color,
      main = "Landmarks congruence matrix",
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n")

abline(v = c(8.5, 17.5, 26.5),
       h = c(8.5, 17.5, 26.5),
       lwd = 3,
       col = "grey")

axis(side = 1, 
     at = 1:dim(congro_overall)[2], 
     labels = lab_congro, 
     las = 2)

axis(side = 2, 
     at = 1:dim(congro_overall)[2], 
     labels = lab_congro, 
     las = 2)

text(x = c(4.5, 13, 22, 31), 
     y = 34.5,
     labels = parts,
     font = 2)

text(y = c(4.5, 13, 22, 31), 
     x = 1.5,
     labels = parts,
     font = 2, 
     srt = 90)

text(x = c(rep(4.5, 3), rep(13, 2), 22), 
     y = c(13, 22, 31, 22, 31, 31),
     labels = rpls_val, 
     pos = 3)

text(x = c(rep(4.5, 3), rep(13, 2), 22), 
     y = c(13, 22, 31, 22, 31, 31),
     labels = p_val,
     pos = 1)

#-------------------------------------------------------------------------------
# 3D deformations

import.tps <- function(x) {
  sc <- scan(x, what = "character")
  
  
}
sc <- scan("data/LM_templates/LM_left_mandi_only.tps", what = "character")
NLM <- sc