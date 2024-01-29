# Asymmetry and bite force study in Schistocerca gregaria
# Data loading, cleaning and exploration

# Required packages:
source(file = "Rfunctions1.txt") # Functions from Claude 2008
library(scales)
library(car)
library(rgl)
library(geomorph)
options(rgl.printRglwidget = TRUE)
#options(browser = "opera") #Can be modified (some bugs with firefox)

#-------------------------------------------------------------------------------

# Load data frame

dat <- read.csv(file = "data/data_subset.csv",
                header = T,
                sep = ",",
                dec = ".")

# Load landmark files

fil <- list.files(path = "data/",
                  pattern = ".tps",
                  full.names = T)

ls <- list()
  
for (i in 1:length(fil)) { # Import TPS coordinate data into a list
  
  sc <- scan(file = fil[i],
             what = "character")
  
  o <- sc[grep("LM=", sc)]
  
  NLM <- as.numeric(sub("LM=", 
                        replacement = "", 
                        o))
  
  if (length(sc) != NLM*3+2) {
    warning("Not the right number of landmarks or dimensions")
  }
  
  coo <- as.numeric(sc[-c(1, length(sc))])
  
  nam <- sub(pattern = "IMAGE=",
             replacement = "", 
             x = sc[length(sc)])
  
  ls[[i]] <- matrix(data = coo, 
                    nrow = NLM, 
                    ncol = 3, 
                    byrow = T)
  
  names(ls)[i] <- nam
  
}

name <- names(ls)

# Make coordinate list into array
arr <- array(data = NA,
             dim = c(dim(ls[[1]]), length(ls)))

for (i in 1:length(ls)) {arr[,,i] <- ls[[i]]}

arr <- arr[-c(8:10),,] 
# LM 9 was not digitizable properly in all specimens
# LM 8 and 10 are imprecise, removed

arr <- arr[,,-grep("0497",
                   x = name)] # Specimen 0497 was an outlier, from another species
dat <- dat[-which(dat$ID == "497"),]

name <- name[-grep("0497",
                   x = name)]

#-------------------------------------------------------------------------------

# Partial generalized Procrustes analysis
pA <- pgpa(A = arr)

# Extract the aligned shapes, and orthogonal projection to Euclidean space
shapes <- orp(pA$rotated)

# Extract centroid sizes
csiz <- pA$cent.size

#-------------------------------------------------------------------------------

# Turn shapes array into matrix
mat_shps <- matrix(data = NA,
                   nrow = dim(shapes)[3],
                   ncol = dim(shapes)[1] * dim(shapes)[2])

for (i in 1:dim(shapes)[3]) {
  
  mat_shps[i,] <- as.vector(t(shapes[,,i]))
  
}

# Define factor to identify shapes and their replicate
replic <- rep(0, 
              length(name))

replic[grep("REPLICATE",
            name)] <- 1

# Define other variables of interest
sex <- as.factor(dat$sex)

body_siz <- dat$body.l

head_w <- dat$head.w
head_l <- dat$head.l
head_h <- dat$head.h

bf <- dat$max_bf
bf2 <- dat$bf.max
mbf <- dat$mean_bf
mbf2 <- dat$bf.mean

#-------------------------------------------------------------------------------
# Shapes PCA with replicates

pca <- prcomp(mat_shps)

var_exp <- pca$sdev^2 / sum(pca$sdev^2) * 100

plot(pca$x[,1:2], 
     pch = c(15, 17)[as.factor(replic)],
     asp = 1)

for (i in 1:(length(replic)/2)) {
  
  lines(pca$x[c(i*2-1, i*2), 1:2])
  
} # Plot lines between replicates to see differences

#-------------------------------------------------------------------------------
# Checking for outliers

plot(shapes[,1,],
     shapes[,2,],
     asp = 1)

plot(shapes[,3,],
     shapes[,2,],
     asp = 1)

#plot3d(shapes[,,1])

#-------------------------------------------------------------------------------
# Average across replicates and re-analysis

name_split <- strsplit(x = name,
                       split = "_")

ID <- rep(NA, length(name))

for (i in 1:length(name)) {
  ID[i] <- name_split[[i]][1]
}

ID <- as.factor(ID)

av_A <- array(data = NA,
              dim = c(dim(shapes)[1:2], 
                      dim(shapes)[3]/length(unique(replic))))

# Average for each individual using mean shape of replicate superimposed 
# configurations 
for (i in 1:length(levels(ID))) {
  
  ind <- levels(ID)[i]
  av_A[,,i] <- mshape(shapes[,, which(ID == ind)])
  
}

# Turn array into a matrix
mat_av <- matrix(data = NA,
                   nrow = dim(av_A)[3],
                   ncol = dim(av_A)[1] * dim(av_A)[2])

for (i in 1:dim(av_A)[3]) {
  mat_av[i,] <- as.vector(t(av_A[,,i]))
}

# Average centroid size across repliactes for each individual.
csiz.av <- tapply(X = csiz, 
                  INDEX = ID,
                  FUN = mean)

# PCA with averaged individuals
pca_av <- prcomp(mat_av)

PC1.av <- pca_av$x[,1]

PC2.av <- pca_av$x[,2]

# Checking for correlations
cor.test(csiz.av, head_h)
cor.test(csiz.av, head_w)
cor.test(csiz.av, head_l)

cor.test(csiz.av, PC1.av)
cor.test(csiz.av, PC2.av)

plot(PC1.av,
     PC2.av,
     pch = c(15, 17)[sex],
     asp = 1)

cor.test(csiz.av, bf)
cor.test(PC1.av, bf)
cor.test(PC2.av, bf)

#-------------------------------------------------------------------------------
# Shape deformations along PC axes

# Make function to obtain max and min deformation along all PCs
# Deformation can be amplified by factor "ampli"
# Plots the first 3 PC deformations if plots = T (open in new device)
# If using 3D data, plots.dims allows to see deformation across orientations
PC.deform <- function(A = A,
                      pca = pca,
                      ampli = 1,
                      plots = T,
                      plots.dims = c(1:2)) {
  
  n <- dim(A)[3]
  p <- dim(A)[1]
  k <- dim(A)[2]
  
  mshp <- mshape(A = A) # Calculate mean shape
  
  v <- as.vector(mshp) # Mean shape coords as vector
  
  pcmax <- apply(X = pca$x, 
                 MARGIN = 2, 
                 FUN = max) # Maximum values along PCs
  
  pcmin <- apply(X = pca$x, 
                 MARGIN = 2, 
                 FUN = min) # Minimum values along PCs
  
  deformed <- list()
  
  for (i in 1:length(pcmax)) {
    
    vrot <- pca_av$rotation[,i]
    
    defmax <- matrix(v + pcmax[i] * vrot,
                     ncol = 3,
                     byrow = F)
    
    defmin <- matrix(v + pcmin[i] * vrot,
                     ncol = 3,
                     byrow = F)
    
    arr <- array(data = NA,
                 dim = c(p, k, 2))
    
    arr[,,1] <- defmax
    arr[,,2] <- defmin
    
    deformed[[i]] <- arr
  }
  
  names(deformed) <- names(pcmax)
  
  if (plots) {
    dev.new()
    plot(mshp[, plots.dims],
         asp = 1,
         main = "PC1 deformation",
         pch = 19,
         cex = 2)
    
    points(deformed$PC1[, plots.dims, 1], 
           pch = 19,
           col = "darkorange2")
    
    points(deformed$PC1[, plots.dims, 2], 
           pch = 19,
           col = "darkorchid")
    
    dev.new()
    plot(mshp[, plots.dims],
         asp = 1,
         main = "PC2 deformation",
         pch = 19, cex = 2)
    
    points(deformed$PC2[, plots.dims, 1], 
           pch = 19,
           col = "darkorange2")
    
    points(deformed$PC2[, plots.dims, 2], 
           pch = 19,
           col = "darkorchid")
    
    dev.new()
    plot(mshp[, plots.dims],
         asp = 1,
         main = "PC3 deformation",
         pch = 19, cex = 2)
    
    points(deformed$PC3[, plots.dims, 1], 
           pch = 19,
           col = "darkorange2")
    
    points(deformed$PC3[, plots.dims, 2], 
           pch = 19,
           col = "darkorchid")
  }
  
  return(deformed)
  
}

deform_av <- PC.deform(A = av_A,
                       pca = pca_av,
                       ampli = 1,
                       plots = T,
                       plots.dims = c(2, 3))

layout(mat = matrix(c(1,1,2,
                      1,1,3,
                      4,5,0), 
                    ncol = 3, 
                    nrow = 3, 
                    byrow = T))

plot(pca_av$x[,1:2], 
     asp = 1,
     pch = c(18,19)[sex],
     cex = 2,
     col = c("darkorange2", "darkorchid4")[sex])

abline(v = 0,
       h = 0,
       col = "grey")

par(mar = c(1,1,1,1))

plot(deform_av$PC2[,2:3,1], 
     pch = 19, 
     cex= 2, 
     asp = 1,
     xlab = "",
     ylab = "",
     bty = "n",
     xaxt = "n",
     yaxt = "n")

plot(deform_av$PC2[,2:3,2], 
     pch = 19, 
     cex= 2, 
     asp = 1,
     xlab = "",
     ylab = "",
     bty = "n",
     xaxt = "n",
     yaxt = "n")

plot(deform_av$PC1[,2:3,2], 
     pch = 19, 
     cex= 2, 
     asp = 1,
     xlab = "",
     ylab = "",
     bty = "n",
     xaxt = "n",
     yaxt = "n")

plot(deform_av$PC1[,2:3,1], 
     pch = 19, 
     cex= 2, 
     asp = 1,
     xlab = "",
     ylab = "",
     bty = "n",
     xaxt = "n",
     yaxt = "n")
