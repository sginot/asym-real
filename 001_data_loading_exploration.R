# Asymmetry and bite force study in Schistocerca gregaria
# Data loading, cleaning and exploration

# Required packages:
source(file = "Rfunctions1.txt")
library(scales)
library(car)
library(rgl)
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
  
for (i in 1:length(fil)) {
  
  sc <- scan(file = fil[i],
             what = "character")
  
  o <- sc[grep("LM=", sc)]
  
  NLM <- as.numeric(sub("LM=", 
                        replacement = "", 
                        o))
  
  if (length(sc) != NLM*3+2) {
    warning("Not the right number of landmarks or dimension")
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

arr <- array(data = NA,
             dim = c(dim(ls[[1]]), length(ls)))

for (i in 1:length(ls)) {arr[,,i] <- ls[[i]]}

arr <- arr[-c(8:10),,] 
# LM 9 was not digitizable properly in all specimens
# LM 8 and 10 are imprecise?

arr <- arr[,,-grep("0497",
                   x = name)] # Specimen 0497 was an outlier, probably another species
dat <- dat[-which(dat$ID == "497"),]

name <- name[-grep("0497",
                   x = name)]

#-------------------------------------------------------------------------------

# Procrustes
pA <- pgpa(A = arr)

#-------------------------------------------------------------------------------

shapes <- pA$rotated
csiz <- pA$cent.size

mat_shps <- matrix(data = NA,
                   nrow = dim(shapes)[3],
                   ncol = dim(shapes)[1] * dim(shapes)[2])

for (i in 1:dim(shapes)[3]) {
  
  mat_shps[i,] <- as.vector(t(shapes[,,i]))
  
}

replic <- rep(0, 
              length(name))

replic[grep("REPLICATE",
            name)] <- 1

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

# Shapes PCA

pca <- prcomp(mat_shps)

var_exp <- pca$sdev^2 / sum(pca$sdev^2) * 100

plot(pca$x[,1:2], 
     pch = c(15, 17)[as.factor(replic)],
     asp = 1)

for (i in 1:length(replic)) {
  
  lines(pca$x[c(i*2-1, i*2), 1:2])
  
}

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
# Checking for correlations

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

for (i in 1:length(levels(ID))) {
  
  ind <- levels(ID)[i]
  av_A[,,i] <- mshape(shapes[,, which(ID == ind)])
  
}

csiz.av <- tapply(X = csiz, 
                  INDEX = ID,
                  FUN = mean)

PC1.av <- tapply(X = pca$x[, 1], 
                 INDEX = ID,
                 FUN = mean)

PC2.av <- tapply(X = pca$x[, 2], 
                 INDEX = ID,
                 FUN = mean)

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

