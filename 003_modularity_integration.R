# Modularity / integration analysis
# Check if patterns of asynnetry (directional) match modules
# Modules: Head capsule / mandibles (both sides = 1 module OR L / R separate)

# Required packages:
library(geomorph)
library(EMMLi)
library(abind)
library(scales)
source("../repo/001_functions.R")
source("../repo/002_asym_components.R")
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

head_mand <- c(1, 1, 2, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 2, 1, 1, 1, 2, 2,
               2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
               2, 2, 2, 2, 2, 2, 1, 1)
# Two modules: 
#The head = 1
#The mandibles (together as one module) = 2


head_mand_sens <- c(1, 1, 2, 1, 3, 3, 3, 1, 1, 1,
                    3, 3, 3, 1, 2, 3, 1, 1, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 3, 3)
# Three modules: 
#The head capsule = 1, 
#The mandibles = 2, 
#The sensory structures = 3

head_mand_asym <- c(1, 1, 2, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 3, 1, 1, 1, 2, 3,
                    2, 2, 2, 3, 3, 3, 2, 3, 2, 2,
                    2, 2, 3, 3, 3, 3, 1, 1)
# Three modules:
# The head capsule = 1
# The right mandible = 2
# The left madible = 3


head_mand_asym_sens <- c(1, 1, 2, 1, 4, 4, 4, 1, 1, 1,
                         4, 4, 4, 1, 3, 4, 1, 1, 2, 3,
                         2, 2, 2, 3, 3, 3, 2, 3, 2, 2,
                         2, 2, 3, 3, 3, 3, 4, 4)
# Four modules:
# The head capsule = 1
# The right mandible = 2
# The left mandible = 3
# Sensory structures = 4

df_modul_models <- data.frame(name = names_LM,
                              no_modul = no_modularity,
                              head_mand = head_mand,
                              head_mand_sens = head_mand_sens,
                              head_mand_asym = head_mand_asym,
                              head_mand_asym_sens = head_mand_asym_sens)

#-------------------------------------------------------------------------------
# Plot the different models

input_folder <- "Figures/"

pdf(file = paste(input_folder, "Module_Models.pdf"),
    width = 6, 
    height = 9)

layout(matrix(c(1:4), 
              ncol = 2,
              byrow = T))

par(mar = c(1, 2, 4, 1))

plot(shapes[,2:3,1], 
     asp = 1, 
     pch = 21, 
     cex = 2,
     bg = head_mand[-c(8:10)],
     main = "Head-Mandibles",
     axes = F)

plot(shapes[,2:3,1], 
     asp = 1, 
     pch = 21, 
     cex = 2,
     bg = head_mand_asym[-c(8:10)],
     main = "Head-Mandibles-Asymmetric",
     axes = F)

plot(shapes[,2:3,1], 
     asp = 1, 
     pch = 21, 
     cex = 2,
     bg = head_mand_sens[-c(8:10)],
     main = "Head-Mandibles-Sensory",
     axes = F)

plot(shapes[,2:3,1], 
     asp = 1, 
     pch = 21, 
     cex = 2,
     bg = head_mand_asym_sens[-c(8:10)],
     main = "Head-Mandibles-Asym-Sensory",
     axes = F)

dev.off()

#-------------------------------------------------------------------------------
#Test run EMMLi

dim(shapes)

A <- av_A #No need for replicates

M <- matrix(NA, 
            nrow = dim(A)[3],
            ncol = dim(A)[2] * dim(A)[1])

for (i in 1:dim(A)[3]) {M[i,] <- c(t(A[,,i]))}

DF <- df_modul_models[-c(8:10),]

DF_3 <- data.frame(apply(X = DF, 
                         MARGIN = 2, 
                         FUN = rep, 
                         each = 3))

DF_3[, 2:6] <- apply(X = DF_3[, 2:6], 
                     MARGIN = 2, 
                     FUN = as.numeric)

cor_M <- cor(M) 

f <- "output_EMMLi_individual_coord.csv"

EMMLi(corr = as.data.frame(cor_M),
      N_sample = nrow(M),
      mod = DF_3,
      abs = T,
      saveAs = f)
#Test here with correlation matrix between individual coordinates


EMMLi(corr = as.data.frame(cor_M),
      N_sample = nrow(M),
      mod = DF,
      abs = T)
#Does not work, correlation must be between LANDMARKS, 
#not individual coordinates

#To produce correlation matrix between LMs -> use congruence coefficient R
#See Goswami & Polly 2010

#Rij = sum((li - ui) . (lj - uj)) / sqrt(sum(li^2) . sum(lj^2))

congruence_coef <- function(A) {
  
  N <- dim(A)[3]
  P <- dim(A)[1]
  K <- dim(A)[2]
  
  mshp <- mshape(A)
  
  mat <- matrix(NA,
                ncol = P,
                nrow = P)
  
  for (i in 1:P) {
    for (j in 1:P) {
      
      ui <- mshp[i,]
      uj <- mshp[j,]
      
      dprod <- rep(NA, N)
      
      for (n in 1:N) {
       
        li <- A[i,,n] 
        lj <- A[j,,n]
        
        dprod[n] <- (li - ui) %*% (lj - uj)
        
      }
      
      num <- sum(dprod)
      
      Li2 <- A[i,,]^2
      Lj2 <- A[j,,]^2
      
      denom <- sqrt(sum(Li2) * sum(Lj2))
      
      Rij <- num / denom
      
      mat[i, j] <- Rij
      
    }
    
  }
return(mat)
}


M <- congruence_coef(av_A)

#Geomorph much simpler tests

modul_test_1 <- modularity.test(A = A, 
                                partition.gp = DF[,3],
                                iter = 999,
                                CI = T)

modul_test_2 <- modularity.test(A = A, 
                                partition.gp = DF[,4],
                                iter = 999,
                                CI = T)

modul_test_3 <- modularity.test(A = A, 
                                partition.gp = DF[,5],
                                iter = 999,
                                CI = T)

modul_test_4 <- modularity.test(A = A, 
                                partition.gp = DF[,6],
                                iter = 999,
                                CI = T)

modul_compar <- compare.CR(modul_test_1, 
                           modul_test_2, 
                           modul_test_3, 
                           modul_test_4)

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


newpart <- DF[,5]
newpart[which(newpart == 1)] <- NA

integ_test_5 <- integration.test(A = A, 
                                 partition.gp = newpart,
                                 iter = 999)

A1 <- A[which(DF[,5] == 2),,]
A2 <- A[which(DF[,5] == 3),,]

tbpls <- two.b.pls(A1 = A1,
                   A2 = A2)

integ_compar <- compare.pls(integ_test_1, 
                            integ_test_2, 
                            integ_test_3, 
                            integ_test_4)
