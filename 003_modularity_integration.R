# Modularity / integration analysis
# Check if patterns of asynnetry (directional) match modules
# Modules: Head capsule / mandibles (both sides = 1 module OR L / R separate)

# Required packages:
library(geomorph)
library(EMMLi)
library(abind)
library(scales)
library(paleomorph)
source("001_functions.R")
source("002_asym_components.R")
source("Rfunctions1.txt")
palette(palette.colors(palette = "Okabe-Ito"))

#-------------------------------------------------------------------------------
# Load landmark template
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
# The left mandible = 2
# The right mandible = 3


head_mand_asym_sens <- c(1, 1, 1, 1, 4, 4, 4, 1, 1, 1,
                         4, 4, 4, 1, 1, 4, 1, 1, 2, 3,
                         2, 2, 2, 3, 3, 3, 2, 3, 2, 2,
                         2, 2, 3, 3, 3, 3, 4, 4)
# Four modules:
# The head capsule = 1
# The left mandible = 2
# The right mandible = 3
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

abbrev_LM <- LM_template$Name

df_modul_models <- data.frame(name = names_LM,
                              no_modul = no_modularity,
                              head_mand = head_mand,
                              head_mand_sens = head_mand_sens,
                              head_mand_asym = head_mand_asym,
                              head_mand_asym_sens = head_mand_asym_sens,
                              ventral_dorsal = ventral_dorsal,
                              half_half = half_half,
                              mandi_only = mandi_only,
                              abbrev = abbrev_LM)

df_modul_models <- df_modul_models[-c(8:10),]

DF <- df_modul_models

rownames(DF) <- 1:dim(DF)[1]
#-------------------------------------------------------------------------------

input_folder <- "Figures/"

#Plot the different views of the landmarks configuration, with lines to show
# the main structures.

template <- shapes[,,1]
  
oculeyeL <- c(6, 7, 34)
oculeyeR <- c(8, 10, 35)

mandiL <- c(16, 18:20, 28, 29, 24, 26, 27)
mandiR <- c(17, 21:23, 32, 33, 25, 31, 30)

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
# Run EMMLi analyses

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

integ_compar <- compare.pls(integ_test_1, 
                            integ_test_2, 
                            integ_test_3, 
                            integ_test_4,
                            integ_test_5,
                            integ_test_6,
                            integ_test_7)

# 2B-PLS analysis between mandibles, and between head / mandi and with BF.
A1 <- A[which(DF[,5] == 2),,] # Left mandible
A2 <- A[which(DF[,5] == 3),,] # Right mandible
A3 <- A[which(DF[,5] == 1),,] # Head

tbpls_mandis <- two.b.pls(A1, A2)
tbpls_LM_head <- two.b.pls(A1, A3)
tbpls_RM_head <- two.b.pls(A2, A3)

M1 <- M2 <- matrix(NA, 
                   nrow = dim(A1)[3],
                   ncol = dim(A1)[1]*dim(A1)[2])

M3 <- matrix(NA, 
             nrow = dim(A3)[3],
             ncol = dim(A3)[1]*dim(A3)[2])

for (i in 1:dim(A1)[3]) {
  M1[i,] <- c(t(A1[,,i]))
  M2[i,] <- c(t(A2[,,i]))
  M3[i,] <- c(t(A3[,,i]))
  }

tbpls_BF_LM <- two.b.pls(bf[-which(is.na(bf))], 
                         M1[-which(is.na(bf)),])
tbpls_BF_RM <- two.b.pls(bf[-which(is.na(bf))], 
                         M2[-which(is.na(bf)),])
tbpls_BF_head <- two.b.pls(bf[-which(is.na(bf))], 
                           M3[-which(is.na(bf)),])

tbpls_BF2_LM <- two.b.pls(bf2[-which(is.na(bf2))], 
                          M1[-which(is.na(bf2)),])
tbpls_BF2_RM <- two.b.pls(bf2[-which(is.na(bf2))], 
                          M2[-which(is.na(bf2)),])
tbpls_BF2_head <- two.b.pls(bf2[-which(is.na(bf2))], 
                            M3[-which(is.na(bf2)),])

#-------------------------------------------------------------------------------
# Graphical representation of covariance patterns between modules,
# à la Burns et al. 2023

# Make function

covar.modules <- function(A, 
                          partition, 
                          part.names = sort(unique(partition)), 
                          LM.names = 1:length(partition),
                          color = hcl.colors(20, 
                                             palette = "viridis",
                                             alpha = NULL, 
                                             rev = FALSE, 
                                             fixup = TRUE)) {
  
  P_fac <- as.factor(partition) # Making sure this is a factor

  # Define dimensions
  n <- dim(A)[3]
  p <- dim(A)[1]
  k <- dim(A)[2]
  
  #Empty matrix for coordinates
  cooM <- matrix(NA,
                  nrow = n,
                  ncol = p * k)
  
  # Fill matrix with coordinates (columns) for each individual (row)
  for (i in 1:dim(A)[3]) {cooM[i,] <- c(t(A[,,i]))}
  
  # Compute covariance and correlation matrices for all coordinates
  covM <- abs(cov(cooM))
  corM <- abs(cor(cooM))
  
  # Compute congruence matrix for LANDMARKS (combined coordinates)
  congruM <- abs(dotcorr(A))
  
  # Reorder matrices so that columns of the same module are together
  o <- order(P_fac)
  congruMo <- congruM[o, o]
  LM.names <- LM.names[o]
  
  oo <- c(t(cbind(o * 3 - 2, 
                  o * 3 - 1, 
                  o * 3)))
  
  covMo <- covM[oo, oo]
  corMo <- corM[oo, oo]

  # Number of columns belonging to each module
  l_part <- table(P_fac)
  l_part_coo <- l_part * 3
  
  # Start plotting
  par(mfrow = c(1, 2),
      mar = c(4, 4, 5, 4))
  
  # Plot full covar matrix (i.e. for coordinates)
  image(x = 1:dim(covMo)[1], 
        y = 1:dim(covMo)[1], 
        z = covMo, 
        col = color,
        main = "Coordinates covariance matrix",
        xlab = "",
        ylab = "",
        xaxt = "n",
        yaxt = "n")
  
  abline(v = cumsum(l_part_coo) + 0.5, 
         h = cumsum(l_part_coo) + 0.5,
         lwd = 3, 
         col = "grey") # Add lines separating the different modules
  
  axis(side = 1,
       at = seq(2, dim(covMo)[1], by = 3),
       labels = LM.names, 
       las = 2,
       cex.axis = 0.7)
  
  axis(side = 2,
       at = seq(2, dim(covMo)[1], by = 3),
       labels = LM.names, 
       las = 2,
       cex.axis = 0.7)
  
  legend(x = p * k + 1, 
         y = p * k + 0.5, 
         bty = "o",
         xpd = T,
         legend = round(seq(min(covMo), 
                            max(covMo), 
                            length = length(color)),
                        7), 
         col = color, 
         pch = 15, 
         cex = 0.7)
  
  cms <- c(0, cumsum(l_part_coo))
  vt <- cms[-length(cms)] + diff(cms) / 2
  
  text(x = vt, 
       y = p * k + 5, 
       labels = part.names,
       xpd = T,
       font = 2)
  
  text(y = vt, 
       x = p * k - 3, 
       labels = part.names, 
       xpd = F, 
       srt = 90, 
       font = 2,
       col ="white")

  # Plot correl matrix for landmarks
  image(x = 1:dim(congruMo)[1], 
        y = 1:dim(congruMo)[1], 
        z = congruMo, 
        col = color,
        main = "Landmarks correlation matrix",
        xlab = "",
        ylab = "",
        xaxt = "n",
        yaxt = "n")
  
  axis(side = 1,
       at = seq(1, length(LM.names), by = 1),
       labels = LM.names, 
       las = 2, 
       cex.axis = 0.7)
  
  axis(side = 2,
       at = seq(1, length(LM.names), by = 1),
       labels = LM.names, 
       las = 2,
       cex.axis = 0.7)
  
  abline(v = cumsum(l_part) + 0.5, 
         h = cumsum(l_part) + 0.5,
         lwd = 3, 
         col = "grey") # Add lines separating modules
  
  legend(x = p + 1,
         y = p + 0.5,
         xpd = T,
         legend = round(seq(min(congruMo), 
                            max(congruMo), 
                            length = length(color)),
                        2), 
         col = color, 
         pch = 15, 
         cex = 0.7)
  
  cms <- c(0, cumsum(l_part))
  vt <- cms[-length(cms)] + diff(cms) / 2
  
  text(x = vt, 
       y = p + 2, 
       labels = part.names,
       xpd = T,
       font = 2)
  
text(y = vt, 
     x = p - 0.5, 
     labels = part.names, 
     xpd = F, 
     font = 2, 
     srt = 90, 
     col = "white")
  
}

#-------------------------------------------------------------------------------
# Produce figures with combined matrices plots

pdf(file = "Figures/head_mandible.pdf", 
    width = 11, 
    height = 6)

covar.modules(A = av_A, 
              partition = DF[,3], 
              part.names = c("head", "mandibles"),
              LM.names = DF$abbrev)
dev.off()

pdf(file = "Figures/head_mand_sens.pdf", 
    width = 11, 
    height = 6)

covar.modules(A = av_A, 
              partition = DF[,4], 
              part.names = c("head", "mandibles", "sensory"),
              LM.names = DF$abbrev)
dev.off()

pdf(file = "Figures/head_mandi_asym.pdf", 
    width = 11, 
    height = 6)

covar.modules(A = av_A, 
              partition = DF[,5], 
              part.names = c("head", "mandi_L", "mandi_R"),
              LM.names = DF$abbrev)
dev.off()

pdf(file = "Figures/head_mandi_asym_sensory.pdf", 
    width = 11, 
    height = 6)

covar.modules(A = av_A, 
              partition = DF[,6], 
              part.names = c("head", "mandi_L", "mandi_R", "sensory"),
              LM.names = DF$abbrev)
dev.off()

pdf(file = "Figures/ventral_dorsal.pdf", 
    width = 11, 
    height = 6)

covar.modules(A = av_A, 
              partition = DF[,7], 
              part.names = c("ventral", "dorsal"),
              LM.names = DF$abbrev)
dev.off()

# Caption for the landmark names
plot(x = rep(3.5, dim(DF)[1]),
     y = 1:dim(DF)[1],
     xlim = c(0, 4),
     pch = "-",
     bty = "n",
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = "")
text(x = rep(3.4, dim(DF)[1]),
     y = 1:dim(DF)[1],
     labels = DF$name,
     cex = 0.8,
     adj = 1)
text(x = rep(3.6, dim(DF)[1]),
     y = 1:dim(DF)[1],
     labels = DF$abbrev,
     cex = 0.8,
     adj = 0)

#-------------------------------------------------------------------------------
# Possible problem with Procrustes for entire dataset (Cardini 2019), 
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

ls_results[[6]] <- modul.intra.gpa(A = A[-which(is.na(DF[,8])),,], 
                                    partition = na.omit(DF[,8]),
                                    plot.compar = T)

ls_results[[7]] <- modul.intra.gpa(A = A[-which(is.na(DF[,9])),,], 
                                   partition = na.omit(DF[,9]),
                                   plot.compar = T)

names(ls_results) <- names(DF)[3:9]

modul_compar_intra_gpa <- compare.CR(ls_results$head_mand[[2]],
                                     ls_results$head_mand_sens[[2]],
                                     ls_results$head_mand_asym[[2]],
                                     ls_results$head_mand_asym_sens[[2]],
                                     ls_results$ventral_dorsal[[2]], 
                                     ls_results$half_half[[2]],
                                     ls_results$mandi_only[[2]],
                                     CR.null = T)

# Compare results from modularity analyses with global gpa, vs local gpa

intra_gpa_Z <- modul_compar_intra_gpa$sample.z[1:7]

global_gpa_Z <- modul_compar$sample.z[1:7]

names(intra_gpa_Z) <- names(global_gpa_Z) <- names(DF)[2:8]

matZ <- rbind(intra_gpa_Z, 
              global_gpa_Z)

# Barplot for comparison

pdf(file = paste(input_folder, "Barplot_ZCR.pdf"),
    height = 8,
    width = 6)

layout(1)
par(mar = c(16,5,2,1))

barplot(height = matZ, 
        beside = T,
        col = c("black", "gray"),
        space = c(0,1),
        las = 2, 
        ylab = "Z_CR",
        names.arg = c("No modularity",
                      "Head-Mandibles",
                      "Head-Mandibles-Sensory",
                      "Head-Mandibles asymmetric",
                      "Head-Mandibles asymmetric-Sensory",
                      "Ventral-Dorsal",
                      "Half-Half"))

dev.off()

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

rownames(matCR_global_gpa) <- names(ls_results)[1:5]

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
     y = matCR_intra_gpa[1:5,1],
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
       y = matCR_intra_gpa[1:5,1],
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

#-------------------------------------------------------------------------------
# Integration after module by module superimposition

rPLS1_head_mand_mbm <- ls_results[[1]][[3]]$r.pls
rPLS2_head_mand_sens_mbm <- ls_results[[2]][[3]]$r.pls.mat
rPLS3_head_mand_asym_mbm <- ls_results[[3]][[3]]$r.pls.mat
rPLS4_head_mand_asym_sens_mbm <- ls_results[[4]][[3]]$r.pls.mat
rPLS5_ventral_dorsal_mbm <- ls_results[[5]][[3]]$r.pls
rPLS6_half_half_mbm <- ls_results[[6]][[3]]$r.pls
rPLS7_mandi_only_mbm <- ls_results[[7]][[3]]$r.pls

rPLS_all_mbm <- list(rPLS1_head_mand_mbm,
                 rPLS2_head_mand_sens_mbm,
                 rPLS3_head_mand_asym_mbm,
                 rPLS4_head_mand_asym_sens_mbm,
                 rPLS5_ventral_dorsal_mbm,
                 rPLS6_half_half_mbm,
                 rPLS7_mandi_only_mbm)


#-------------------------------------------------------------------------------
# Plot

palette(palette.colors(palette = "Okabe-Ito")[2:7])

titles <- c("A. Head-Mandibles", 
            "B. Head-Mandibles-Sensory",
            "C. Head-Mandibles asymmetric",
            "D. Head-Mandibles asymmetric-Sensory",
            "E. Ventral-Dorsal",
            "F. Half-Half")

pdf(file = paste(input_folder, 
                 "pairwise_integration_PLS.pdf", 
                 sep = ""),
    height = 13,
    width = 7)

layout(matrix(c(1:6), 
              ncol = 2,
              byrow = T))

par(mar = c(1, 1, 3, 1))

for (i in 1:6) {
  
  plot(x = template[,2],
       y = template[,3],
       asp = 1, 
       pch = 21, 
       cex = 2.5,
       bg = DF[, i+2],
       main = titles[i],
       axes = F)
  
  m_cent <- apply(template, 
                  2, 
                  tapply, 
                  INDEX = as.factor(DF[, i+2]), 
                  mean)[,2:3]
  
  cor_mat <- as.matrix(rPLS_all[[i]])
  cor_mat_mbm <- as.matrix(rPLS_all_mbm[[i]])
  
  if (length(cor_mat) == 1) {
    cor_mat <- matrix(c(0, cor_mat, cor_mat, 0), 
                      ncol = 2, byrow = T)
    cor_mat_mbm <- matrix(c(0, cor_mat_mbm, cor_mat_mbm, 0), 
                      ncol = 2, byrow = T)
    
    colnames(cor_mat_mbm) <- colnames(cor_mat) <- rownames(m_cent)
    rownames(cor_mat_mbm) <- rownames(cor_mat) <- rownames(m_cent)
  } 
  
  for (j in rownames(cor_mat)) {
    for (k in colnames(cor_mat)) {
      
      m_line <- rbind(m_cent[j,],
                      m_cent[k,])
      
      lines(m_line,
            lwd = cor_mat[j,k] * 20,
            col = alpha("gray", 0.3))
    }
  }
  
  points(m_cent,
         cex = 4,
         pch = 21,
         bg = rownames(m_cent)) 
  
  for (j in rownames(cor_mat)) {
    for (k in colnames(cor_mat)) {
      
      m_line <- rbind(m_cent[j,],
                      m_cent[k,])
      
      if (round(cor_mat[j,k], 2) == 0) txt_col <- alpha("white", alpha = 0)
      else txt_col <- "black"
        
      text(apply(m_line, 2, mean)[1],
           apply(m_line, 2, mean)[2],
           labels = round(cor_mat[j,k], 2),
           cex = 1.2,
           pos = 3,
           col = txt_col)
      
      text(apply(m_line, 2, mean)[1],
           apply(m_line, 2, mean)[2],
           labels = round(cor_mat_mbm[j,k], 2),
           cex = 0.9,
           col = txt_col)
    }
    
  }
  
}

dev.off()

#-------------------------------------------------------------------------------
# Plot relationship between rPLS with global superimposition vs module by module

P_values_pairwise <- c(integ_test_1$P.value,
                       integ_test_2$pairwise.P.values,
                       integ_test_3$pairwise.P.values,
                       integ_test_4$pairwise.P.values,
                       integ_test_5$P.value,
                       integ_test_6$P.value,
                       integ_test_7$P.value)

P_values_pairwise_mbm <- c(ls_results[[1]][[3]]$P.value,
                           ls_results[[2]][[3]]$pairwise.P.values,
                           ls_results[[3]][[3]]$pairwise.P.values,
                           ls_results[[4]][[3]]$pairwise.P.values,
                           ls_results[[5]][[3]]$P.value,
                           ls_results[[6]][[3]]$P.value,
                           ls_results[[7]][[3]]$P.value)

integ_compar_intra_gpa <- compare.pls(ls_results[[1]][[3]],
                                      ls_results[[2]][[3]],
                                      ls_results[[3]][[3]],
                                      ls_results[[4]][[3]],
                                      ls_results[[5]][[3]],
                                      ls_results[[6]][[3]],
                                      ls_results[[7]][[3]])

X <- unlist(rPLS_all)
Y <- unlist(rPLS_all_mbm)

mod_p_val <- lm(Y ~ X)

newx <- seq(min(X), 
            max(X), 
            length.out = 100)

preds <- predict(mod_p_val,
                 newdata = data.frame(X = newx),
                 interval = "confidence")

pdf(file = paste(input_folder, 
                 "rPLS_superimpositions_correlation", 
                 sep = ""),
    height = 6,
    width = 6)

plot(X, Y,  pch = NA,
     xlim = c(0.2,1),
     ylim = c(0.2,1),
     xlab = "Global superimposition pairwise rPLS",
     ylab = "Module by moudle superimposition pairwise rPLS",)

abline(0,1, 
       lty = 2, 
       lwd = 2, 
       col = "gray")

clip(x1 = min(X),
     x2 = max(X),
     y1 = 0,
     y2 = 1)

abline(mod_p_val,
       lwd = 2)

polygon(x = c(newx, rev(newx)), 
        y = c(preds[,2], rev(preds[,3])),
        border = NA,
        col = alpha("black", 0.2))

clip(x1 = 0,
     x2 = 2,
     y1 = 0,
     y2 = 2)

points(X, 
       Y,
       cex = 2,
       pch = as.numeric(P_values_pairwise < 0.05 
                        & P_values_pairwise_mbm < 0.05)+ 21,
       bg = as.numeric(P_values_pairwise < 0.05 
                       & P_values_pairwise_mbm < 0.05) + 1)

legend("topleft", 
       legend = c("Only significant with global superimposition",
                  "Significant with both superimposition approaches"),
       pch = c(21, 22),
       pt.bg = c(1, 2),
       bty = "n",
       pt.cex = 2)

dev.off()
