# Modularity / integration analysis
# Check if patterns of asynnetry (directional) match modules
# Modules: Head capsule / mandibles (both sides = 1 module OR L / R separate)

# Required packages:
library(geomorph)
library(EMMLi)
library(abind)
library(scales)
library(paleomorph)
source("../asym_simulation/001_functions.R")
source("../asym_simulation/002_asym_components.R")
source("Rfunctions1.txt")
palette(palette.colors(palette = "Okabe-Ito"))

#-------------------------------------------------------------------------------
# Load landöark template
LM_template <- read.csv("data/LM_template.csv")

#Make modularity models dataframe

names_LM <- LM_template[,1]

no_modularity <- rep(1, 
                     length(names_LM))
# No modularity: all landmarks are in only one module

head_mand <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
               2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
               2, 2, 2, 2, 2, 2, 1, 1)
# Two modules: 
#The head = 1
#The mandibles (together as one module) = 2


head_mand_sens <- c(1, 1, 1, 1, 4, 4, 4, 1, 1, 1,
                    4, 4, 4, 1, 1, 4, 1, 1, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 4, 4)
# Three modules: 
#The head capsule = 1, 
#The mandibles = 2, 
#The sensory structures = 4

head_mand_asym <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 2, 3,
                    2, 2, 2, 3, 3, 3, 2, 3, 2, 2,
                    2, 2, 3, 3, 3, 3, 1, 1)
# Three modules:
# The head capsule = 1
# The right mandible = 2
# The left madible = 3


head_mand_asym_sens <- c(1, 1, 1, 1, 4, 4, 4, 1, 1, 1,
                         4, 4, 4, 1, 1, 4, 1, 1, 2, 3,
                         2, 2, 2, 3, 3, 3, 2, 3, 2, 2,
                         2, 2, 3, 3, 3, 3, 4, 4)
# Four modules:
# The head capsule = 1
# The right mandible = 2
# The left mandible = 3
# Sensory structures = 4

ventral_dorsal <- c(1, 1, 1, 1, 4, 4, 4, 1, 1, 1,
                    4, 4, 4, 1, 1, 4, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 4, 4)
# Two modules:
# The dorsal part of the head (eyes, ocelli, antennae...)
# The mandibles, associated with ventral part of the head (tentorium...)

half_half <- c(NA, NA, 1, 1, 1, 1, 1, 1, NA, 2,
               2, 2, 2, 2, 2, NA, 1, 2, 1, 2, 
               1, 1, 1, 2, 2, 2, 1, 2, 1, 1,
               1, 1, 2, 2, 2, 2, 1, 2)

# Two modules:
# Left side of head and mandibles
# Right side
# Midline landmarks excludede

mandi_only <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, NA, NA, NA, NA, NA, 2, 3,
                2, 2, 2, 3, 3, 3, 2, 3, 2, 2,
                2, 2, 3, 3, 3, 3, NA, NA)

# Two modules:
# Left mandible
# Right mandible


df_modul_models <- data.frame(name = names_LM,
                              no_modul = no_modularity,
                              head_mand = head_mand,
                              head_mand_sens = head_mand_sens,
                              head_mand_asym = head_mand_asym,
                              head_mand_asym_sens = head_mand_asym_sens,
                              ventral_dorsal = ventral_dorsal,
                              half_half = half_half,
                              mandi_only = mandi_only)

df_modul_models <- df_modul_models[-c(8:10),]

DF <- df_modul_models
#-------------------------------------------------------------------------------

input_folder <- "Figures/"

#Plot the different views of the landmarks configuration, with lines to show
# the main structures.

template <- shapes[,,1]
  
oculeyeR <- c(6, 7, 34)
oculeyeL <- c(8, 10, 35)

mandiR <- c(16, 18:20, 28, 29, 24, 26, 27)
mandiL <- c(17, 21:23, 32, 33, 25, 31, 30)

pdf(file = paste(input_folder, "LM_config.pdf"),
    width = 12, 
    height = 6)

layout(matrix(c(1:3), 
              ncol = 3,
              byrow = T))

par(mar = c(1, 2, 4, 1))

plot(template[,c(1,2)], 
     asp = 1, 
     pch = 21, 
     cex = 2,
     bg = "firebrick",
     main = "Dorsal view",
     axes = F)

text(template[,c(1,2)], 
     labels = c(1:35), 
     pos = c(rep(1,10), 
             3, 2, 1, 2, 2, 1, 1, 1, 3, 3,
             2, 3, 1, 4, 4, 4, 2, 3, 1, 4, 
             2, 1, 3, 1, 1))

#polygon(template[oculeyeL,
#                 c(1,2)])
#polygon(template[oculeyeR,
#                 c(1,2)])
#polygon(template[mandiL,
#                 c(1,2)])
#polygon(template[mandiR,
#                 c(1,2)])

plot(template[,c(2,3)], 
     asp = 1, 
     pch = 21, 
     cex = 2,
     bg = "forestgreen",
     main = "Anterior view",
     axes = F)

text(template[,c(2,3)], 
     labels = c(1:35), 
     pos = c(rep(1,18),
             4, 3, 1, 2, 2, 4, 2, 3, 1, 4, 1, 1, 3, 4, 1))

#polygon(template[oculeyeL,
#               c(2,3)])
#polygon(template[oculeyeR,
#                 c(2,3)])
#polygon(template[mandiL,
#                 c(2,3)])
#(template[mandiR,
#                 c(2,3)])


plot(template[,c(1,3)], 
     asp = 1, 
     pch = 21, 
     cex = 2,
     bg = "navyblue",
     main = "Lateral view",
     axes = F)

text(template[,c(1,3)], 
     labels = c(1:35), 
     pos = c(1,1,1,3,3,3,3,1,1,1,
             2,1,1,4,2,3,1,1,1,1,
             1,1,1,4,1,3,3,1,1,4,
             1,1,1,1,1))

dev.off()


# Plot the different models
palette(palette.colors(palette = "Okabe-Ito")[2:7])
#palette(c("bisque", "firebrick", "navyblue", "darkolivegreen2"))

pdf(file = paste(input_folder, "Module_Models.pdf"),
    width = 6, 
    height = 12)

layout(matrix(c(1:6), 
              ncol = 2,
              byrow = T))

par(mar = c(1, 2, 4, 1))

plot(x = -template[,2],
     y = template[,3],
     asp = 1, 
     pch = 21, 
     cex = 3,
     bg = DF[,3],
     main = "Head-Mandibles",
     axes = F)

plot(x = -template[,2],
     y = template[,3],
     asp = 1, 
     pch = 21, 
     cex = 3,
     bg = df_modul_models[,4],
     main = "Head-Mandibles-Sensory",
     axes = F)

plot(x = -template[,2],
     y = template[,3],
     asp = 1, 
     pch = 21, 
     cex = 3,
     bg = df_modul_models[,5],
     main = "Head-Mandibles-Asym",
     axes = F)

plot(x = -template[,2],
     y = template[,3],
     asp = 1, 
     pch = 21, 
     cex = 3,
     bg = df_modul_models[,6],
     main = "Head-Mandibles-Asym-Sensory",
     axes = F)

plot(x = -template[,2],
     y = template[,3],
     asp = 1, 
     pch = 21, 
     cex = 3,
     bg = df_modul_models[,7],
     main = "Ventral-Dorsal",
     axes = F)

plot(x = -template[,2],
     y = template[,3],
     asp = 1, 
     pch = 21, 
     cex = 3,
     bg = df_modul_models[,8],
     main = "Half-Half",
     axes = F)

#plot(x = -template[,2],
#     y = template[,3],
#     asp = 1, 
#     pch = 21, 
#     cex = 2,
#     bg = df_modul_models[,9],
#     main = "Mandibles-only",
#     axes = F)

dev.off()

#-------------------------------------------------------------------------------
#Test run EMMLi

dim(shapes)

A <- av_A #No need for replicates

M <- matrix(NA, 
            nrow = dim(A)[3],
            ncol = dim(A)[2] * dim(A)[1])

for (i in 1:dim(A)[3]) {M[i,] <- c(t(A[,,i]))}

DF_3 <- data.frame(apply(X = DF, 
                         MARGIN = 2, 
                         FUN = rep, 
                         each = 3))

DF_3[, 2:6] <- apply(X = DF_3[, 2:6], 
                     MARGIN = 2, 
                     FUN = as.numeric)

cor_M <- cor(M)

cong_M <- dotcorr(A)

f <- "output_EMMLi_individual_coord.csv"
f2 <- "output_EMMLi_LM_correl.csv"

EMMLi(corr = as.data.frame(cor_M),
      N_sample = nrow(M),
      mod = DF_3[,1:6],
      abs = T,
      saveAs = f)
#Test here with correlation matrix between individual coordinates


EMMLi(corr = as.data.frame(cong_M),
      N_sample = nrow(M),
      mod = DF[,1:6],
      abs = T,
      saveAs = f2)

#To produce correlation matrix between LMs -> use congruence coefficient R
#See Goswami & Polly 2010
#Not necessary, already a function available in paleomorph

#Rij = sum((li - ui) . (lj - uj)) / sqrt(sum(li^2) . sum(lj^2))

#congruence_coef <- function(A) {
  
#  N <- dim(A)[3]
#  P <- dim(A)[1]
#  K <- dim(A)[2]
#  
#  mshp <- mshape(A)
#  
#  mat <- matrix(NA,
#                ncol = P,
#                nrow = P)
#  
#  for (i in 1:P) {
#    for (j in 1:P) {
#      
#      ui <- mshp[i,]
#      uj <- mshp[j,]
#      
#      dprod <- rep(NA, N)
#      
#      for (n in 1:N) {
#       
#        li <- A[i,,n] 
#        lj <- A[j,,n]
#        
#        dprod[n] <- (li - ui) %*% (lj - uj)
#        
#      }
#      
#      num <- sum(dprod)
#      
#      Li2 <- A[i,,]^2
#      Lj2 <- A[j,,]^2
#      
#      denom <- sqrt(sum(Li2) * sum(Lj2))
#      
#      Rij <- num / denom
#      
#      mat[i, j] <- Rij
#      
#    }
#    
#  }
#return(mat)
#}


#M <- congruence_coef(av_A)

#-------------------------------------------------------------------------------
# Geomorph CR and Z score for modularity comparisons
# EMMLi tends to artificially over-support modularity partitions with more 
# modules. This is apparently not the case with the Adams 2016 approach.

# Head-nmandibles 2 modules partition
modul_test_1 <- modularity.test(A = A, 
                                partition.gp = DF[,3],
                                iter = 999,
                                CI = T)

# Ventral head structures, sensory (dorsal) and mandibles partition 3 modules
modul_test_2 <- modularity.test(A = A, 
                                partition.gp = DF[,4],
                                iter = 999,
                                CI = T)

# Head as a whole, mandibles separate partition (3 modules)
modul_test_3 <- modularity.test(A = A, 
                                partition.gp = DF[,5],
                                iter = 999,
                                CI = T)

# Head ventral, head sensory and mandibles separate (4 modules)
modul_test_4 <- modularity.test(A = A, 
                                partition.gp = DF[,6],
                                iter = 999,
                                CI = T)

# Mandibles AND head ventral together, vs head dorsal sensory (2 modules)
modul_test_5 <- modularity.test(A = A, 
                                partition.gp = DF[,7],
                                iter = 999,
                                CI = T)

# Left and right sides of the head as 2 modules. Midline LMs removed.
modul_test_6 <- modularity.test(A = A[-which(is.na(DF[,8])),,], 
                                partition.gp = na.omit(DF[,8]),
                                iter = 999,
                                CI = T)

# Left and right mandibles only (2 modules) all other LM removed
modul_test_7 <- modularity.test(A = A[-which(is.na(DF[,9])),,], 
                                partition.gp = na.omit(DF[,9]),
                                iter = 999,
                                CI = T)

ls_modul <- list(modul_test_1,
                 modul_test_2,
                 modul_test_3,
                 modul_test_4,
                 modul_test_5)

modul_compar <- compare.CR(modul_test_1, 
                           modul_test_2, 
                           modul_test_3, 
                           modul_test_4,
                           modul_test_5,
                           modul_test_6,
                           modul_test_7,
                           CR.null = T)

integ_test_1 <- integration.test(A = A, 
                                partition.gp = DF[,3],
                                iter = 999)

integ_test_2 <- integration.test(A = A, 
                                partition.gp = DF[,4],
                                iter = 999)

integ_test_3 <- integration.test(A = A, 
                                partition.gp = DF[,5],
                                iter = 999)

integ_test_4 <- integration.test(A = A, 
                                partition.gp = DF[,6],
                                iter = 999)

integ_test_5 <- integration.test(A = A, 
                                 partition.gp = DF[,7],
                                 iter = 999)

integ_test_6 <- integration.test(A[-which(is.na(DF[,8])),,], 
                                 partition.gp = na.omit(DF[,8]),
                                 iter = 999)

integ_test_7 <- integration.test(A[-which(is.na(DF[,9])),,], 
                                 partition.gp = na.omit(DF[,9]),
                                 iter = 999)


A1 <- A[which(DF[,5] == 2),,]
A2 <- A[which(DF[,5] == 3),,]

tbpls <- two.b.pls(A1 = A1,
                   A2 = A2)

integ_compar <- compare.pls(integ_test_1, 
                            integ_test_2, 
                            integ_test_3, 
                            integ_test_4,
                            integ_test_5,
                            integ_test_6,
                            integ_test_7)

# Problem with Procrustes for entire dataset (Cardini 2019), 
# test modul/integ with module-specific alignment.
# However, this removes the spatial and size relationship between "real" LMs.

# Make function to align module by modul and output an array of same dimensions

modul.intra.gpa <- function(A, partition, plot.compar = T) {
  
  mods <- unique(partition)
  
  A_intra_gpa <- array(data = NA,
                       dim = dim(A))
  
  for (i in 1:length(mods)) {
    
    A_mod <- A[which(partition == mods[i]),,]
    
    centroid_mod <- apply(A_mod, 2, mean)
    
    A_mod_gpa <- pgpa(A_mod)$rotated
    
    A_intra_gpa[which(partition == mods[i]),,] <- A_mod_gpa
      
    #for (j in 1:dim(A_mod_gpa)[3]) {
      
    #  A_intra_gpa[which(partition == mods[i]),,j] <- 
    #    t(t(A_mod_gpa[,,j]) + centroid_mod)
      
    #}
    
  }
  
  if (plot.compar) {
    
    layout(matrix(1:2, 
                  ncol =2))
    
    plot(x = A[,2,],
         y = A[,3,],
         pch = 19,
         cex = 0.5,
         col = partition, 
         asp = 1)
    
    plot(x = A_intra_gpa[,2,],
         y = A_intra_gpa[,3,],
         pch = 19,
         cex = 0.5,
         col = partition, 
         asp = 1)
  }
  
  mod_test <- modularity.test(A_intra_gpa,
                              partition.gp = partition,
                              CI = T)
  
  int_test <- integration.test(A_intra_gpa,
                               partition.gp = partition)
  
  return(list(A_intra_gpa, mod_test, int_test))

}

# Apply function to 5 "really possible" modular partitions.

ls_results <- list()

for (i in 1:5) {
  
  ls_results[[i]] <- modul.intra.gpa(A,
                                     partition = DF[,i+2],
                                     plot.compar = T)
}

names(ls_results) <- names(DF)[3:7]

modul_compar_intra_gpa <- compare.CR(ls_results$head_mand[[2]],
                                     ls_results$head_mand_sens[[2]],
                                     ls_results$head_mand_asym[[2]],
                                     ls_results$head_mand_asym_sens[[2]],
                                     ls_results$ventral_dorsal[[2]], 
                                     CR.null = T)

# Compare results from modularity analyses with global gpa, vs local gpa

intra_gpa_Z <- modul_compar_intra_gpa$sample.z

global_gpa_Z <- modul_compar$sample.z[1:6]

names(intra_gpa_Z) <- names(global_gpa_Z) <- names(DF)[2:7]

matZ <- rbind(intra_gpa_Z, 
              global_gpa_Z)

# Barplot for comparison
par(mar = c(11,3,2,3))

barplot(height = matZ, 
        beside = T,
        col = 1:2,
        space = c(0,2),
        las = 2)

# Barplot ordered by global gpa Z score

matZo <- matZ[, order(intra_gpa_Z)]

par(mar = c(11,3,2,3))

barplot(height = matZo, 
        beside = T,
        col = 1:2,
        space = c(0,2),
        las = 2)

# CR score comparison with CI

matCR_intra_gpa <- matrix(NA, ncol = 3, 
                          nrow = length(ls_results))

for (i in 1:length(ls_results)) {
  
  matCR_intra_gpa[i, 1] <- ls_results[[i]][[2]]$CR
  matCR_intra_gpa[i, 2:3] <- ls_results[[i]][[2]]$CInterval
  
}

rownames(matCR_intra_gpa) <- names(ls_results)

#---

matCR_global_gpa <- matrix(NA, ncol = 3, 
                           nrow = 5)

rownames(matCR_global_gpa) <- names(ls_results)

matCR_global_gpa[1,] <- c(modul_test_1$CR, 
                          modul_test_1$CInterval)

matCR_global_gpa[2,] <- c(modul_test_2$CR, 
                          modul_test_2$CInterval)

matCR_global_gpa[3,] <- c(modul_test_3$CR, 
                          modul_test_3$CInterval)

matCR_global_gpa[4,] <- c(modul_test_4$CR, 
                          modul_test_4$CInterval)

matCR_global_gpa[5,] <- c(modul_test_5$CR, 
                          modul_test_5$CInterval)

#---

plot(x = 1:5,
     y = matCR_global_gpa[,1],
     pch = 21,
     bg = 1,
     cex =2,
     ylim = c(0.4, 1))

points(x = 1:5 + 0.1,
     y = matCR_intra_gpa[,1],
     pch = 21,
     bg = 2,
     cex =2)

for (i in 1:5) {
  
  lines(x = rep(i, 2),
        y = matCR_global_gpa[i,2:3],
        lwd = 2,
        col = 1)
  
  lines(x = rep(i, 2) + 0.1,
        y = matCR_intra_gpa[i,2:3],
        lwd = 2,
        col = 2)
}

# Something appears wrong with the CI produced when using local alignment
# Maybe due to something in geomorph, but probably rather due to the 
# bootstrapping across modules which are not in their original positions?

plot(x = 1:5,
     y = matCR_global_gpa[,1],
     pch = 21,
     bg = 1,
     cex =2,
     ylim = c(0.4, 1))

points(x = 1:5 + 0.1,
       y = matCR_intra_gpa[,1],
       pch = 21,
       bg = 2,
       cex =2)

for (i in 1:5) {
  
  lines(x = rep(i, 2),
        y = c(min(ls_modul[[i]]$CR.boot),
              max(ls_modul[[i]]$CR.boot)),
        lwd = 2,
        col = 1)
  
  lines(x = rep(i, 2) + 0.1,
        y = c(min(ls_results[[i]][[2]]$CR.boot),
              max(ls_results[[i]][[2]]$CR.boot)),
        lwd = 2,
        col = 2)
}

#-------------------------------------------------------------------------------
# Bootstrapping to get 95CI for Z_CR values

bootZ <- function(A, partition.gp, it = 100) {
  
  Z <- list()
  
  i <- 1
  
  while (i < it + 1) {
    
    index <- sample(1:dim(A)[3], 
                    replace = T)
    
    Ab <- A[,,index]
    
    test <- modularity.test(A = Ab, 
                            partition.gp = partition.gp,
                            iter = 99,
                            CI = F)
    
    Z[[i]] <- test$Z
    
    i <- i + 1
  }

Z}

ls_bootZ <- lapply(1:length(ls_modul), 
                   function(i) bootZ(A = A, 
                                     partition.gp = DF[, i+2],
                                     it = 1000))

plot(1:5, matZ[2,2:6], ylim = c(0,-8))
for (i in 1:5){lines(rep(i,2), 
                       quantile(unlist(ls_bootZ[[i]]), 
                                c(0.025,0.975)))}

ls_bootZ_intra <- lapply(1:length(ls_results), 
                   function(i) bootZ(A = A, 
                                     partition.gp = DF[, i+2],
                                     it = 1000))

save.image()

#-------------------------------------------------------------------------------
# Integration measure between modules (different partitions but only global GPA)

rPLS1_head_mand <- integ_test_1$r.pls
rPLS2_head_mand_sens <- integ_test_2$r.pls.mat
rPLS3_head_mand_asym <- integ_test_3$r.pls.mat
rPLS4_head_mand_asym_sens <- integ_test_4$r.pls.mat
rPLS5_ventral_dorsal <- integ_test_5$r.pls
rPLS6_half_half <- integ_test_6$r.pls
rPLS7_mandi_only <- integ_test_7$r.pls

rPLS_all <- list(rPLS1_head_mand,
                 rPLS2_head_mand_sens,
                 rPLS3_head_mand_asym,
                 rPLS4_head_mand_asym_sens,
                 rPLS5_ventral_dorsal,
                 rPLS6_half_half,
                 rPLS7_mandi_only)


pdf(file = paste(input_folder, 
                 "pairwise_integration_PLS.pdf", 
                 sep = ""),
    height = 12,
    width = 5)

layout(matrix(c(1:6), 
              ncol = 2,
              byrow = T))

par(mar = c(1, 2, 4, 1))

for (i in 1:6) {
  
  plot(x = template[,2],
       y = template[,3],
       asp = 1, 
       pch = 21, 
       cex = 1.5,
       bg = DF[, i+2],
       main = colnames(DF)[i+2],
       axes = F)
  
  m_cent <- apply(template, 
                  2, 
                  tapply, 
                  INDEX = as.factor(DF[, i+2]), 
                  mean)[,2:3]
  
  cor_mat <- as.matrix(rPLS_all[[i]])
  
  if (length(cor_mat) == 1) {
    cor_mat <- matrix(c(0, cor_mat, cor_mat, 0), 
                      ncol = 2, byrow = T)
    
    colnames(cor_mat) <- rownames(m_cent)
    rownames(cor_mat) <- rownames(m_cent)
  } 
  
  for (j in rownames(cor_mat)) {
    for (k in colnames(cor_mat)) {
      
      m_line <- rbind(m_cent[j,],
                      m_cent[k,])
      
      lines(m_line,
            lwd = cor_mat[j,k] * 20,
            col = alpha("gray", 0.3))
      
      text(apply(m_line, 2, mean)[1],
           apply(m_line, 2, mean)[2],
           labels = round(cor_mat[j,k], 2),
           cex = 1.5)
      
    }
    
    points(m_cent,
           cex = 4,
           pch = 21,
           bg = rownames(m_cent))
    
  }
  
}

dev.off()

