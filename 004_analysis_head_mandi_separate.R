# Re-analysis of the data following reviews
# Main change will be the separate analysis of mandibles and head
# For the asymmetry part this means switching from object symmetry to matching
# symmetry. For the modularity and integration part this means mostly focusing 
# covariation patterns.

# Run scripts 001 to 003 before this one.

#-------------------------------------------------------------------------------
# Required packages
#-------------------------------------------------------------------------------

source(file = "Rfunctions1.txt") # Functions from Claude 2008
library(scales)
library(car)
library(rgl)
library(geomorph)
library(EMMLi)
library(paleomorph)

#-------------------------------------------------------------------------------
# Raw data organisation and superimposition
#-------------------------------------------------------------------------------

LM_template[-c(8:10), 1] # Shows the list of landmarks, to identify the ones 
                        # belonging to each mandible
dim(arr) # Array containing raw coordinates, with 3 noisy LMs removed 
        #(see script 001).

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
# of the array therefore belong to the same individual.

pA_mand <- pgpa(A_mand) # Partial Procrustes sumperimposition

shp_mand <- orp(pA_mand$rotated) # Orthogonal projection into Euclidean space
mshp_mand <- pA_mand$mshape

side <- as.factor(rep(c("R", "L"), 
                      length(ID)))

mshp_RM <- mshape(shp_mand[,, which(side == "R")])
# Mean shape of shapes with odd index -> Right mandible

mshp_LM <- mshape(shp_mand[,, which(side == "L")])
# Mean shape of shapes with even index -> Left mandible (mirrored)

pA_head <- pgpa(A_head) # Partial Procrustes sumperimposition
shp_head <- orp(pA_head$rotated) # Orthogonal projection into Euclidean space
mshp_head <- pA_head$mshape # Mean shape of head LMs

ID_head <- ID # Recycle vector from script 001
ID_mand <- rep(ID, each = 2) # Same vector, with duplicated elements because
                            # each individual has two mandibles.

replic_head <- replic
replic_mand <- rep(replic, each = 2)

#-------------------------------------------------------------------------------
# Define output folder for plots
#-------------------------------------------------------------------------------

output_folder <- "./Figures"

#-------------------------------------------------------------------------------
# Define modules for the split up coordinates matrices
#-------------------------------------------------------------------------------

modu_head <- DF[which(DF[,6] == 1 | DF[,6] ==  4), 
                c(6, 10)] # 2-modules split for the head

modo_head <- modu_head[order(modu_head[,1]),] # Reordered landmarks

part_overall <- as.factor(c(rep("head", 17), 
                            rep("LM", 9), 
                            rep("RM", 9)))

part_overall2 <- as.factor(c(modu_head[,1], 
                             rep("LM", 9), 
                             rep("RM", 9)))

levels(part_overall2)[1:2] <- c("head", "sensory")

#-------------------------------------------------------------------------------
# For analyses not requiring replicates, shapes are averaged across replicates
# for each individual
# Make matrices averaging replicates:
#-------------------------------------------------------------------------------

# individual:side for the mandible, left and right mandibles of one individual
# should not be averaged together!

mat_shp_mand <- matrix(NA, 
                       nrow = dim(shp_mand)[3], 
                       ncol = dim(shp_mand)[1] * dim(shp_mand)[2])

for (i in 1:dim(mat_shp_mand)[1]) {mat_shp_mand[i,] <- c(t(shp_mand[,,i]))}

mat_shp_mand_av <- apply(X = mat_shp_mand, 
                         MARGIN = 2, 
                         FUN = tapply, 
                         INDEX = ID_mand:side, 
                         mean)

# individual average for the head

mat_shp_head <- matrix(NA, 
                       nrow = dim(shp_head)[3], 
                       ncol = dim(shp_head)[1] * dim(shp_head)[2])

for (i in 1:dim(mat_shp_head)[1]) {mat_shp_head[i,] <- c(t(shp_head[,,i]))}

mat_shp_head_av <- apply(X = mat_shp_head, 
                         MARGIN = 2, 
                         FUN = tapply, 
                         INDEX = ID_head, 
                         mean)

# Make global coordinate matrix, combining left and right mandible with head for
# each individual (each row represents one individual)

mat_overall <- matrix(NA, 
                      nrow = dim(mat_shp_head_av)[1], 
                      ncol = dim(mat_shp_head_av)[2] + 
                        dim(mat_shp_mand_av)[2] * 2)

for (i in 1:dim(mat_overall)[1]) {
  
  mat_overall[i, 1:51] <- mat_shp_head_av[i, ]
  
  mat_overall[i, 52:78] <- 
    mat_shp_mand_av[grep(pattern = ":L", 
                         x = rownames(mat_shp_mand_av))[i],]
  
  mat_overall[i, 79:105] <- 
    mat_shp_mand_av[grep(pattern = ":R", 
                         x = rownames(mat_shp_mand_av))[i],]
}

rownames(mat_overall) <- rownames(mat_shp_head_av)

#-------------------------------------------------------------------------------
# Asymmetry analyses
#-------------------------------------------------------------------------------

# For mandibles, use matching symmetry

bilatsym_mand <- bilat.symmetry(A = shp_mand,
                                ind = ID_mand,
                                side = side,
                                replicate = replic_mand,
                                object.sym = F)

# For head, must define pairs of landmarks to use object symmetry analysis

# Plot and text to assess visually the pairs on LMs
plot(shp_head[,2,], 
     shp_head[,3,], 
     asp = 1)

text(shp_head[,2,1], 
     shp_head[,3,1], 
     col = 2, 
     pch = 19, 
     cex = 2)

# Matching pairs of landmarks are the following:

LM_pairs <- list()

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

# Run symmetry analyses from geomorph

bilatsym_head <- bilat.symmetry(A = shp_head,
                                ind = ID_head,
                                replicate = replic_head,
                                object.sym = T,
                                land.pairs = LM_pairs)

plot(bilatsym_head)

# Alternative if 3D plotting via RGL does not work:

#bilatsym_head2D <- bilat.symmetry(A = shp_head[,2:3,],
#                                ind = ID_head,
#                                replicate = replic_head,
#                                object.sym = T,
#                                land.pairs = LM_pairs)

#plot(bilatsym_head2D)


#Test asymmetry in head dorsal/sensory landmarks only

shp_head_sens <- shp_head[which(modu_head == 4),,]

bilatsym_head_sens <- bilat.symmetry(A = shp_head_sens,
                                     ind = ID_head,
                                     replicate = replic_head,
                                     object.sym = T,
                                     land.pairs = matrix(c(1,5, 2,4, 3,6, 8,9), 
                                                         ncol = 2, 
                                                         byrow = T),
                                     RRPP = T)

#Test asymmetry in head ventral landmarks only

shp_head_ventral <- shp_head[which(modu_head == 1),,]

bilatsym_head_ventral <- bilat.symmetry(A = shp_head_ventral,
                                     ind = ID_head,
                                     replicate = replic_head,
                                     object.sym = T,
                                     land.pairs = matrix(c(3,6, 4,5, 7,8), 
                                                         ncol = 2, 
                                                         byrow = T),
                                     RRPP = T)

#-------------------------------------------------------------------------------
# Test for 2BPLS relation between shape components and bf
#-------------------------------------------------------------------------------

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
# Then test for integration between each mandible and with head modules if 
# significant head modularity is found
#-------------------------------------------------------------------------------

mod_test_head <- modularity.test(A = shp_head,
                                 partition.gp = modu_head[,1],
                                 iter = 999, 
                                 CI = T)

integ_test_head <- integration.test(A = shp_head,
                                    partition.gp = modu_head[,1],
                                    iter = 999)

integ_test_3mod <- integration.test(A = A_overall,
                                    partition.gp = part_overall,
                                    iter = 999)

integ_test_4mod <- integration.test(A = A_overall,
                                    partition.gp = part_overall2,
                                    iter = 999)
#-------------------------------------------------------------------------------
# Compute covariance / correlation patterns
#-------------------------------------------------------------------------------

# Compute covariance and correlation matrices for each variable, ie each 
# landmark coordinate is considered independently

cov_overall <- abs(cov(mat_overall)) # Absolute values bc sign does not matter
cor_overall <- abs(cor(mat_overall))

# Compute congruence coefficient correlation matrix, ie correlation matrix 
# between landmarks (combining their 3 dimensions)

A_overall <- array(NA,
                   dim = c(dim(mat_overall)[2] / 3, 
                           3, 
                           dim(mat_overall)[1])) # Create empty array

for (i in 1:dim(A_overall)[3]) {
  
  A_overall[,,i] <- matrix(mat_overall[i,], 
                           ncol = 3, 
                           byrow = T)
  } # Fill array with values, ie reconstitute a configuration array, but with
    # modified order of landmarks, as found in mat_overall

congru_overall <- abs(dotcorr(A_overall))

# Reorder LMs in the congruence matrix so that they are grouped according to
# the two head modules.
congro_overall <- congru_overall
congro_overall[1:17, 1:17] <- congro_overall[order(modu_head[,1]),
                                             order(modu_head[,1])]

m <- congro_overall
m[upper.tri(congro_overall, diag = T)] <- NA # Make half empty matrix

#-------------------------------------------------------------------------------
# Covariance / correlation matrix heatmap plots
#-------------------------------------------------------------------------------

# Define color scale 
color <- hcl.colors(30, 
                    palette = "viridis",
                    alpha = NULL, 
                    rev = FALSE, 
                    fixup = TRUE)


# Covariance matrix, with individual coordinates as variables

pdf(file = paste(output_folder, 
                 "covariance_heatmap.pdf",
                 sep = "/"),
    height = 8,
    width = 8)

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

dev.off()

# Correlation matrix, with individual coordinates as variables

pdf(file = paste(output_folder, 
                 "correlation_heatmap.pdf",
                 sep = "/"),
    height = 8,
    width = 8)

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

dev.off()

# Congruence matrix, with landmarks congruence coefficients as variables

pdf(file = paste(output_folder, 
                 "congruence_heatmap.pdf",
                 sep = "/"),
    height = 8,
    width = 8)

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

dev.off()

# Congruence matrix plot with reordered LMs, and integration test results

lab_congro <- lab_congru
lab_congro[1:17] <- lab_congru[order(modu_head[,1])] # Order head LMs names by 
                                                    # the module they belong to

rpls_val <- paste("r-PLS =", 
                  round(integ_test_4mod$r.pls.mat, 3))

p_val <- paste("P =", integ_test_4mod$pairwise.P.values)

parts <- levels(part_overall2)

pdf(file = paste(output_folder, 
                 "LM_correlation_module_congruence.pdf",
                 sep = "/"),
    height = 8,
    width = 11)

par(mar = c(4, 8, 6, 4))

image(x = 1:dim(m)[1], 
      y = 1:dim(m)[1], 
      z = m, 
      col = color,
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n")

title(main = "Landmarks congruence matrix",
      line = 5)

abline(v = c(8.5, 17.5, 26.5),
       h = c(8.5, 17.5, 26.5),
       lwd = 3,
       col = "grey")

axis(side = 1, 
     at =  c(4.5, 13, 22, 31), 
     labels = parts, 
     font = 2,
     las = 1)

axis(side = 2, 
     at =  c(4.5, 13, 22, 31), 
     labels = parts, 
     font = 2, las = 3)

axis(side = 4, 
     at = 1:dim(congro_overall)[2], 
     labels = lab_congro, 
     las = 2)

axis(side = 3, 
     at = 1:dim(congro_overall)[2], 
     labels = lab_congro, 
     las = 2)

legend(x = -5,
       y = 35.5,
       xpd = T,
       legend = round(seq(min(na.omit(c(m))), 
                          max(na.omit(c(m))), 
                          length = length(color)), 2),
       col = color,
       pch = 15)

text(x = c(rep(4.5, 3), rep(13, 2), 22), 
     y = c(13, 22, 31, 22, 31, 31),
     labels = rpls_val, 
     pos = 3)

text(x = c(rep(4.5, 3), rep(13, 2), 22), 
     y = c(13, 22, 31, 22, 31, 31),
     labels = p_val,
     pos = 1)

dev.off()
# Possible artifact from alignment of left and right mandibles. May artificially
# increase correlation between homologous landmarks


#-------------------------------------------------------------------------------
# Check if asymmetric components from head and mandi are correlated
#-------------------------------------------------------------------------------

tbpls_asym_head_mand <- two.b.pls(bilatsym_head$asymm.shape,
                                  bilatsym_mand$asymm.shape)

# Gather both mandi and head asym component shapes in an array

A_asym <- array(NA,
                dim = c(dim(bilatsym_head$asymm.shape)[1] + 
                          dim(bilatsym_mand$asymm.shape)[1],
                        dim(bilatsym_head$asymm.shape)[2:3]))

for (i in 1:dim(A_asym)[3]) {
  A_asym[1:17,,i] <- bilatsym_head$asymm.shape[,,i]
  A_asym[18:26,,i] <- bilatsym_mand$asymm.shape[,,i]
}

# Compute congruence matrix
congru_asym <- abs(dotcorr(A_asym))

# Plot resulting heatmap
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
#-------------------------------------------------------------------------------

tbpls_sym_head_mand <- two.b.pls(bilatsym_head$symm.shape,
                                 bilatsym_mand$symm.shape)
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
# 3D deformations: compute landmark templates and meshes
#-------------------------------------------------------------------------------

# Make function to import single tps files

import.single.tps <- function(x) {
  
  sc <- scan(x, what = "character")
  
  raw_coo <- as.numeric(sc[-grep(pattern = "=",
                                 x = sc)])
  
  p <- as.numeric(sub(pattern = "LM=",
                      replacement = "", 
                      x = sc[grep(pattern = "LM=", sc)]))
  
  k <- length(raw_coo) / p
  
  m <- matrix(data = raw_coo,
              nrow = p,
              ncol = k,
              byrow = T)
}

# Import 3D templates (Note they are not in the same coordinate system as 
# aligned shapes used up until here)

LM_templ <- import.single.tps("data/LM_templates/LM_left_mandi_only.tps")
RM_templ <- import.single.tps("data/LM_templates/LM_right_mandi_only.tps")
head_templ <- import.single.tps("data/LM_templates/LM_head_only.tps")

# 3D plots to check import went fine

plot3d(LM_templ, aspect = T)
plot3d(RM_templ, aspect = T)
plot3d(head_templ, aspect = T)

# Import 3D PLY meshes (ASCII)

mesh_head <- read.ply("head_mesh.ply")
mesh_LM <- read.ply("left_mandible_mesh.ply")
mesh_RM <- read.ply("right_mandible_mesh.ply")

# Warp meshes by using landmarks digitized on them and warping them to target
# configurations

refMesh_head <- warpRefMesh(mesh = mesh_head, 
                            mesh.coord = head_templ, 
                            ref = mshp_head, 
                            color = 2,
                            centered = F) # No problem for the head

# For mandibles, since there is strong directional asymmetry, the right mandible
# mesh was chosen arbitrarily, because the left mandible configurations were 
# mirrored in the beginning to match right mandible configurations
refMesh_RM <- warpRefMesh(mesh = mesh_RM, 
                          mesh.coord = RM_templ, 
                          ref = mshape(bilatsym_mand$symm.shape), 
                          color = 2,
                          centered = F)


#-------------------------------------------------------------------------------
# Compare mean shapes of right and left mandi, then compute and add or substract 
# the L-R difference matrix N times, to show amplified differences
#-------------------------------------------------------------------------------

# Mandible shape difference matrix. Reminder: Left mandible configurations were 
# mirrored to match right configurations
diff_mand <- mshp_LM - mshp_RM

# Amplify differences by adding/subtracting difference matrix
x2_LM <- mshp_LM + 2 * diff_mand

x2_RM <- mshp_RM - 2 * diff_mand

# Warp reference mesh to the deformed configurations

mesh_mshp_LM <- warpRefMesh(mesh = mesh_LM, 
                            mesh.coord = LM_templ, 
                            ref = mshp_LM, 
                            color = "red1",
                            centered = F)

mesh_x2_LM <- warpRefMesh(mesh = mesh_LM, 
                            mesh.coord = LM_templ, 
                            ref = x2_LM, 
                            color = "firebrick",
                            centered = F)

mesh_mshp_RM <- warpRefMesh(mesh = mesh_RM, 
                            mesh.coord = RM_templ, 
                            ref = mshp_RM, 
                            color = "cadetblue1",
                            centered = F)

mesh_x2_RM <- warpRefMesh(mesh = mesh_RM, 
                          mesh.coord = RM_templ, 
                          ref = x2_RM, 
                          color = "cadetblue3",
                          centered = F)




#-------------------------------------------------------------------------------
# Do the equivalent for the head, mirror mean shape then add or substract the 
# difference matrix N times, to show amplified differences
#-------------------------------------------------------------------------------

# When mirroring the mean head configuration, the landmarks must be reordered
# so the new config matches the original one

# To define the reordering of LMs viamanually via plotting
#reord_head <- locate.reorder(shape = mshp_head, along = 2)

# Alternatively, the order is the following:
reord_head <- c(1,2,12,11,9,8,10,6,5,7,4,3,13,15,14,17,16)

# Prepare reordered matrix#
mirror_mshp_head <- mshp_head[reord_head,] 

# Mirror along the medio lateral axis
mirror_mshp_head[,2] <- mirror_mshp_head[,2] * -1

# Compute a symmetrical head configuration
sym_head <- (mshp_head + mirror_mshp_head) / 2

# Compute the head configuration difference matrix
diff_head <- mshp_head - mirror_mshp_head

# Compute amplified differences configurations
x2_head <- sym_head + 2 * diff_head

x2_mirror_head <- sym_head - 2 * diff_head

# Warp meshes to the different configurations computed
mesh_mshp_head <- warpRefMesh(mesh = mesh_head, 
                              mesh.coord = head_templ, 
                              ref = mshp_head, 
                              color = "red1",
                              centered = F)

mesh_x2_head <- warpRefMesh(mesh = mesh_head, 
                            mesh.coord = head_templ, 
                            ref = x2_head, 
                            color = "firebrick",
                            centered = F)

mesh_mirror_mshp_head <- warpRefMesh(mesh = mesh_head, 
                                     mesh.coord = head_templ, 
                                     ref = mirror_mshp_head, 
                                     color = "cadetblue1",
                                     centered = F)

mesh_x2_mirror_head <- warpRefMesh(mesh = mesh_head, 
                                   mesh.coord = head_templ, 
                                   ref = x2_mirror_head, 
                                   color = "cadetblue3",
                                   centered = F)

mesh_sym_head <- warpRefMesh(mesh = mesh_head, 
                             mesh.coord = head_templ, 
                             ref = sym_head, 
                             color = "grey",
                             centered = F)

#-------------------------------------------------------------------------------
# Output 3D deformation plots: mandibles
#-------------------------------------------------------------------------------

# Define the viewpoints to standardize the position of 3D meshes on the plots

antview <- matrix(c(-0.1, 1, 0, 0,
                    0.4, 0.1, 0.9, 0,
                    0.9, 0.2, -0.4, 0,
                    0, 0, 0, 1), 
                  ncol = 4,
                  byrow = T)

dorsview <- matrix(c(0.06,  0.89, -0.44,  0.00,
                     -0.99, -0.01, -0.14,  0.00,
                     -0.13,  0.45,  0.89,  0.00,
                     0.00,  0.00,  0.00,  1.00), 
                   ncol = 4,
                   byrow = T)

postview <- matrix(c(-0.28, -0.24, -0.93, 0.00,
                     -0.95, -0.05, 0.30,  0.00,
                     -0.12,  0.97, -0.21,  0.00,
                     0.00,  0.00,  0.00,  1.00), 
                   ncol = 4,
                   byrow = F)

# Anterior view of deformed mandibles

par3d(windowRect = c(50, 50, 3000, 1000),
      userMatrix = antview)

newSubscene3d(newviewport = c(0, 0, 900, 900))
plot3d(mesh_x2_LM, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(x2_LM, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(600, 0, 900, 900))
plot3d(mesh_mshp_LM, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(mshp_LM, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1200, 0, 900, 900))
plot3d(mesh_mshp_RM, 
       box = F,
       axes = F, 
       xlab = "", 
       ylab = "", 
       zlab = "")

plot3d(mshp_RM, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1800, 0, 900, 900))
plot3d(mesh_x2_RM, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "", 
       zlab = "")

plot3d(x2_RM, 
       add = T,
       type = "s",
       size = 1.5)

rgl.snapshot(filename = paste(output_folder,
                              "mandi_asym_deform_anterior.png",
                              sep = "/"))

close3d()

# Dorsal view mandi deformation

par3d(windowRect = c(50, 50, 3000, 1000),
      userMatrix = dorsview)

newSubscene3d(newviewport = c(0, 0, 900, 900))
plot3d(mesh_x2_LM, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(x2_LM, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(600, 0, 900, 900))
plot3d(mesh_mshp_LM, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(mshp_LM, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1200, 0, 900, 900))
plot3d(mesh_mshp_RM, 
       box = F,
       axes = F, 
       xlab = "", 
       ylab = "", 
       zlab = "")

plot3d(mshp_RM, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1800, 0, 900, 900))
plot3d(mesh_x2_RM, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "", 
       zlab = "")

plot3d(x2_RM, 
       add = T,
       type = "s",
       size = 1.5)

rgl.snapshot(filename = paste(output_folder,
                              "mandi_asym_deform_dorsal.png",
                              sep = "/"))

close3d()

# Posterior view mandi deformation

par3d(windowRect = c(50, 50, 3000, 1000),
      userMatrix = round(postview,2))

newSubscene3d(newviewport = c(0, 0, 900, 900))
plot3d(mesh_x2_LM, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(x2_LM, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(600, 0, 900, 900))
plot3d(mesh_mshp_LM, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(mshp_LM, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1200, 0, 900, 900))
plot3d(mesh_mshp_RM, 
       box = F,
       axes = F, 
       xlab = "", 
       ylab = "", 
       zlab = "")

plot3d(mshp_RM, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1800, 0, 900, 900))
plot3d(mesh_x2_RM, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "", 
       zlab = "")

plot3d(x2_RM, 
       add = T,
       type = "s",
       size = 1.5)

rgl.snapshot(filename = paste(output_folder,
                              "mandi_asym_deform_post.png",
                              sep = "/"))

close3d()

#-------------------------------------------------------------------------------
# Output 3D deformation plots: head
#-------------------------------------------------------------------------------

# Anterior views of deformed head

antview2 <- matrix(c(0.03, 0.07,  0.99,  0.00,
                     1.02,  0.01,  0.09,  0.00,
                     0.01,  0.98, -0.07,  0.00,
                     0.00,  0.00,  0.00,  1.00),
                   ncol = 4)

par3d(windowRect = c(50, 50, 3000, 1000),
      userMatrix = antview2)

newSubscene3d(newviewport = c(-400, -500, 1500, 2000))
plot3d(mesh_x2_mirror_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "", 
       xlim = c(-0.5, 0.5),
       ylim = c(-0.5, 0.5),
       zlim = c(-0.5, 0.5))

plot3d(x2_mirror_head, 
       add = T,
       type = "s",
       size = 1)

newSubscene3d(newviewport = c(0, -500, 1500, 2000))
plot3d(mesh_mirror_mshp_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "", 
       xlim = c(-0.5, 0.5),
       ylim = c(-0.5, 0.5),
       zlim = c(-0.5, 0.5))

plot3d(mirror_mshp_head, 
       add = T,
       type = "s",
       size = 1)

newSubscene3d(newviewport = c(400, -500, 1500, 2000))
plot3d(mesh_sym_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "", 
       xlim = c(-0.5, 0.5),
       ylim = c(-0.5, 0.5),
       zlim = c(-0.5, 0.5))


plot3d(sym_head, 
       add = T,
       type = "s",
       size = 1)

newSubscene3d(newviewport = c(800, -500, 1500, 2000))
plot3d(mesh_mshp_head, 
       box = F,
       axes = F, 
       xlab = "", 
       ylab = "", 
       zlab = "", 
       xlim = c(-0.5, 0.5),
       ylim = c(-0.5, 0.5),
       zlim = c(-0.5, 0.5))

plot3d(mshp_head, 
       add = T,
       type = "s",
       size = 1)

newSubscene3d(newviewport = c(1200, -500, 1500, 2000))
plot3d(mesh_x2_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "", 
       zlab = "", 
       xlim = c(-0.5, 0.5),
       ylim = c(-0.5, 0.5),
       zlim = c(-0.5, 0.5))

plot3d(x2_head, 
       add = T,
       type = "s",
       size = 1)

rgl.snapshot(filename = paste(output_folder,
                              "head_asym_deform_anterior.png",
                              sep = "/"))

close3d()

# Dorsal views of head deformation

dorsview2 <- matrix(c(-0.04, -0.98, -0.18,  0.00,
                      1.00, -0.04,  0.00,  0.00,
                      0.01, -0.18,  0.99,  0.00,
                      0.00,  0.00,  0.00,  1.00),
                    ncol = 4)

par3d(windowRect = c(50, 50, 3000, 1000),
      userMatrix = dorsview2)

newSubscene3d(newviewport = c(-200, 0, 1000, 1000))
plot3d(mesh_x2_mirror_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(x2_mirror_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(300, 0, 1000, 1000))
plot3d(mesh_mirror_mshp_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(mirror_mshp_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(800, 0, 1000, 1000))
plot3d(mesh_sym_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(sym_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1300, 0, 1000, 1000))
plot3d(mesh_mshp_head, 
       box = F,
       axes = F, 
       xlab = "", 
       ylab = "", 
       zlab = "")

plot3d(mshp_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1800, 0, 1000, 1000))
plot3d(mesh_x2_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "", 
       zlab = "")

plot3d(x2_head, 
       add = T,
       type = "s",
       size = 1.5)

rgl.snapshot(filename = paste(output_folder,
                              "head_asym_deform_dorsal.png",
                              sep = "/"))

close3d()

# Posterior views of head deformation

postview2 <- matrix(c(0.10, -0.15, -0.99,  0.00,
                      -0.99,  0.02, -0.10,  0.00,
                      0.03,  0.99, -0.14,  0.00,
                      0.00, 0.00,  0.00, 1.00),
                    ncol = 4)

par3d(windowRect = c(50, 50, 3000, 1000),
      userMatrix = postview2)

newSubscene3d(newviewport = c(-300, -500, 1500, 2000))
plot3d(mesh_x2_mirror_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "", 
       xlim = c(-0.5, 0.5),
       ylim = c(-0.5, 0.5),
       zlim = c(-0.5, 0.5))

plot3d(x2_mirror_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(100, -500, 1500, 2000))
plot3d(mesh_mirror_mshp_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "", 
       xlim = c(-0.5, 0.5),
       ylim = c(-0.5, 0.5),
       zlim = c(-0.5, 0.5))

plot3d(mirror_mshp_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(500, -500, 1500, 2000))
plot3d(mesh_sym_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "", 
       xlim = c(-0.5, 0.5),
       ylim = c(-0.5, 0.5),
       zlim = c(-0.5, 0.5))

plot3d(sym_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(900, -500, 1500, 2000))
plot3d(mesh_mshp_head, 
       box = F,
       axes = F, 
       xlab = "", 
       ylab = "", 
       zlab = "", 
       xlim = c(-0.5, 0.5),
       ylim = c(-0.5, 0.5),
       zlim = c(-0.5, 0.5))

plot3d(mshp_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1300, -500, 1500, 2000))
plot3d(mesh_x2_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "", 
       zlab = "", 
       xlim = c(-0.5, 0.5),
       ylim = c(-0.5, 0.5),
       zlim = c(-0.5, 0.5))

plot3d(x2_head, 
       add = T,
       type = "s",
       size = 1.5)

rgl.snapshot(filename = paste(output_folder,
                              "head_asym_deform_posterior.png",
                              sep = "/"))

close3d()

# Ventral view of head asymmetric deformations

ventview <- matrix(c(0.03,  1.00, -0.14,  0.00,
                     0.99, -0.03,  0.01,  0.00,
                     0.01, -0.14, -0.99,  0.00,
                     0.00,  0.00,  0.00,  1.00),
                   ncol = 4)

par3d(windowRect = c(50, 50, 3000, 1000),
      userMatrix = ventview)

newSubscene3d(newviewport = c(-200, 0, 1000, 1000))
plot3d(mesh_x2_mirror_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(x2_mirror_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(300, 0, 1000, 1000))
plot3d(mesh_mirror_mshp_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(mirror_mshp_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(800, 0, 1000, 1000))
plot3d(mesh_sym_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "",
       zlab = "")

plot3d(sym_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1300, 0, 1000, 1000))
plot3d(mesh_mshp_head, 
       box = F,
       axes = F, 
       xlab = "", 
       ylab = "", 
       zlab = "")

plot3d(mshp_head, 
       add = T,
       type = "s",
       size = 1.5)

newSubscene3d(newviewport = c(1800, 0, 1000, 1000))
plot3d(mesh_x2_head, 
       box = F, 
       axes = F, 
       xlab = "",
       ylab = "", 
       zlab = "")

plot3d(x2_head, 
       add = T,
       type = "s",
       size = 1.5)

rgl.snapshot(filename = paste(output_folder,
                              "head_asym_deform_ventral.png",
                              sep = "/"))

close3d()
