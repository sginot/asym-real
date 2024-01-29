#-------------------------------------------------------------------------------
# This script gathers functions developed to run simulations of shapes with
# different amounts of Directional asymmetry, Fluctuating asymmetry, etc.
# Functions also include a way to relabel landmarks semi-automatically and a
# function to produce a multivariate asymmetry decomposition based on a method
# proposed by Neubauer et al. 2020.

#-------------------------------------------------------------------------------
# The general aim is to test a new way of computing fluctuating and directional
# asymmetry at the individual level. One cannot actually compute FA and DA at  
# the individual level since the calculation is based on variance decomposition.
# The idea here is to compute an approximation of individual asymmetries by
# 1) computing FA and DA for individual landmarks at the population level,
# 2) computing individual total asymmetry (ie the difference in landmark 
# position) at each landmark,
# Decomposing that total asymmetry based on proportions of FA and DA at the
# respective landmarks

# The validation and assessment of that process can be achieved by simulation.

# This script must be run before all others in the repository

#Required packages:
library(abind)
source("Rfunctions1.txt")

#-------------------------------------------------------------------------------
# Function simul.LM.asym simulates bias (directional asymmetry) and noise
# (fluctuating asymmetry and error) for ONE landmark.

simul.LM.asym <- function(coords.LM = c(0, 0, 0), 
                          # x y (z) coordinates for landmark
                          bias = rep(0.2, length(coords.LM)),
                          # bias applied to the coordinates
                          sd.bias = abs(0.1 * bias),
                          noise = abs(0.2 * bias)) {
  #noise by default 20% of bias, can also be
  #arbitrarily defined. This is passed to the 
  #sd argument in rnorm function.
  
  k <- length(coords.LM)
  
  bias_LM <- rep(NA, k)
    
  for (i in 1:k) {
    
    bias_LM[i] <- coords.LM[i] + rnorm(n = 1,
                                       mean = bias[i],
                                       sd = sd.bias[i])
    
  } 
  # Add the defined biases to coordinates of the targeted landmark
  
  new_LM <- rep(NA, k)
  # The vector for the new landmark is here defined.
  
  for (i in 1:k) {
  
    new_LM[i] <- rnorm(n = 1, 
                       mean = bias_LM[i], 
                       sd = noise[i])
    # Add some noise around the biased x and y coordinates to simulate FA
    # A normal distribution is therefore assumed for FA, which may be in the
    # same direction or a different direction as the bias.  
  }
  
  return(list(bias_LM, new_LM))
  # The value returned by the function is a list containing two elements:
  # The biased LM, and the biased + noisy LM
}


#-------------------------------------------------------------------------------
# Function simulsym simulates a given number of shapes by applying recursively
# the simul.LM.asym.

#-------------------------------------------------------------------------------
# End results wanted: 
# Array with N simulated shapes combining DA and FA
# Matrix with exact amount of DA, FA and error for each individual shape

# Parameters should be:
# Basic shape
# Number of simulated shapes to produce
# Number of simulated replicates
# Average bias in x and y for each landmark (ie DA)
# Standard deviation in the bias
# Amount of noise (FA) to add, in proportion of the bias
# Amount of noise (error) to add to replicates

#-------------------------------------------------------------------------------
# Define shape simulation function 

simulsym <- function(shape = stop("Basic symmetric shape not defined"),
                     N = 100,
                     # Number of shapes to be simulated
                     Nrep = 1,
                     # Number of replicate per shape
                     m.bias = stop("Matrix of biases must be defined"), 
                     # matrix of biases applied to the coordinates
                     m.sd.bias = abs(0.1 * m.bias),
                     # matrix of variation in bias, by default 10% of m.bias
                     m.noise = abs(0.2 * m.bias),
                     # matrix of amount of noise, by default 20% of m.bias
                     error.replic = abs(0.05 * mean(m.bias))
                     # error applied to all landmarks in all dimensions
                     ) {
  
  p <- dim(shape)[1]
  
  k <- dim(shape)[2]
  
  Niter <- N
  # Define the number of iterations of the simulation of landmarks. 
  # This number corresponds to the desired number of shapes in the end.
  
  arr_simul_DA <- array(NA, 
                        dim = c(p, 
                                k, 
                                Niter))
  arr_simul_DA_FA <- array(NA, 
                           dim = c(p,
                                   k, 
                                   Niter))
  arr_simul_err <- array(NA, 
                         dim = c(p,
                                 k,
                                 Niter * Nrep))
  # Initiate empty arrays to gather simulated shapes.
  
  it <- 1
  # Initialize i value
  
  while (it <= Niter) {
    # This will run 100 times by default
    
    sim_shp_DA <- sim_shp_FA <- shape
    # Define a matrix of the same dimensions as the original shape to gather 
    # simulated values for all landmarks of ONE simulated shape
    
    for (i in 1:p) {
      # Loop to apply simulation function to all landmarks, creating a full shape
      
      sim <- simul.LM.asym(coords.LM = shape[i, ],
                           bias = m.bias[i, ],
                           sd.bias = m.sd.bias[i, ], 
                           noise = m.noise[i, ])
      sim_shp_DA[i, ] <- sim[[1]]
      sim_shp_FA[i, ] <- sim[[2]]
      # Runs function, based on biases values defined previously.
      
    }
    
    arr_simul_DA[,, it] <- sim_shp_DA
    arr_simul_DA_FA[,, it] <- sim_shp_FA
    
    for (j in 1:Nrep) {
      
      arr_simul_err[,, (Niter * j) - Niter + it] <- 
        sim_shp_FA + rnorm(n = length(sim_shp_FA),
                           mean = 0,
                           sd = error.replic)
      # At the moment should not work for Nrep > 1
      # Assign iterated shape to the array
      
    }
    
    it <- it + 1
    
  }
  
  M <- matrix(NA, 
              nrow = N,
              ncol = 6)
  
  colnames(M) <- c("DAi", 
                   "FAi", 
                   "Ei", 
                   "distDAi", 
                   "distFAi", 
                   "distEi")
  
  m <- matrix(NA,
              nrow = N,
              ncol = p * k)
  
  for (i in 1:N) {
    
    m[i,] <- c(arr_simul_DA_FA[,, i] - shape)
    
    M[i, 1] <- sum(arr_simul_DA[,, i] - shape)
    M[i, 2] <- sum(arr_simul_DA_FA[,, i] - arr_simul_DA[,, i])
    M[i, 3] <- sum(arr_simul_err[,, i] - arr_simul_DA_FA[,, i])
    # At the moment should not work for Nrep > 1
    M[i, 4] <- sqrt(sum((arr_simul_DA[,, i] - shape)^2))
    M[i, 5] <- sqrt(sum((arr_simul_DA_FA[,, i] - arr_simul_DA[,, i])^2))
    M[i, 6] <- sqrt(sum((arr_simul_err[,, i] - arr_simul_DA_FA[,, i])^2))
    
  }
  
  return(list(A = abind(arr_simul_DA_FA, arr_simul_err),
              DA = arr_simul_DA,
              M = data.frame(M),
              m = m))
}

#-------------------------------------------------------------------------------
# Define function to obtain graphically the reordering of points

locate.reorder <- function(shape = stop("Original shape must be defined"),
                           along = 1) {
  
  k <- ncol(shape)
  
  shape.mirror <- shape
  shape.mirror[, along] <- shape[, along] * -1
  # This mirrors the shape along the given axis. Therefore it assumes that the
  # the given axis is the axis of bilateral symmetry
  
  if (k == 2) {
    
    plot(shape,
         asp = 1,
         xlim = c(min(c(shape[, 1],
                        shape.mirror[, 1])), 
                  max(c(shape[, 1],
                        shape.mirror[, 1]))),
         ylim = c(min(c(shape[, 2],
                        shape.mirror[, 2])), 
                  max(c(shape[, 2],
                        shape.mirror[, 2]))),
         type = "n")
    
    text(shape, 
         labels = 1:nrow(shape),
         cex = 2,
         pos = 1)
    
    points(shape.mirror,
           pch = 21,
           bg = "red")
    
    print("Click on the red points in the order indicated by numbers")
    
    ids <- identify(shape.mirror,
                    n = nrow(shape.mirror),
                    order = T)
    
    return(ids$order)
    
  }
  
  if (k == 3) {

    plot3d(shape,
           type = "p")
    
    
    text3d(shape, 
           texts = 1:nrow(shape),
           cex = 2,
           pos = 1)
    
    plot3d(shape.mirror,
           add = T,
           type = "s",
           col = "red",
           radius = 0.001)
    
    print("Click on the red points in the order indicated by numbers")
    
    ids <- identify3d(x = shape.mirror,
                      n = nrow(shape.mirror),
                      plot = T,
                      adj = c(0, 1),
                      buttons = c("right", "middle"))
    
    return(ids)
  }
  
}

#-------------------------------------------------------------------------------
# mv.asym function applies multivariate asym decomposition à la Neubauer et al.

# End result wanted:
# Matrix with individual FA, DA and error
# PCA of symmetric variation
# PCA of asymmetric variation

# Parameters should be:
# Array containing shapes and their replicates
# Vector defining reordering of landmarks when mirrored
# Factor defining individual replicates

mv.asym <- function(A = stop("Array with shapes must be defined"),
                    reorder.LM = stop("Vector reodering landmarks is needed"),
                    Nrep = 1,
                    along = 1,
                    indiv.fac = as.factor(c(1:(dim(A)[3]/(Nrep + 1)), 
                                            1:(dim(A)[3]/(Nrep + 1))))
                    # By default, indiv.fac assumes that every individual in the 
                    # sample has been replicated, and in the same order as the 
                    # first series of landmarking. Otherwise it can be manually
                    # defined.
                    ) {
  
  p <- dim(A)[1]
  k <- dim(A)[2]
  N <- dim(A)[3]/(Nrep + 1)
  
  Am <- A
  Am[, along, ] <- A[, along, ] * -1
  # Make mirrored shapes array
  
  Am <- Am[reorder.LM, , ]
  # Relabel landmarks so they are in the same order as original shape.
  
  Ac <- abind(A, Am)
  # Join originally oriented and mirrored arrays
  
  mirror.fac <- as.factor(c(rep("ORIGINAL", dim(A)[3]), 
                            rep("MIRROR", dim(A)[3])))
  # Define orientation factor for the global array
  
  indiv.fac.mirror <- rep(indiv.fac,
                          (Nrep + 1))
  # Define individual factor for the global array
  
  pAc <- pgpa(Ac)
  # Partial Procrustes analysis from Claude (2008)
  
  sym_shapes <- array(NA,
                      dim = c(p, k, N))
  # Initiate empty array
  
  for (i in 1:N) {
    
    index <- which(indiv.fac.mirror == levels(indiv.fac.mirror)[i])
    
    sym_shapes[ , , i] <- mshape(orp(pAc$rotated)[ , , index])
    
  } # Fill array with mean shape per individual across orientations and replicates

  mat_symshapes <- matrix(data = sym_shapes, 
                          nrow = N, 
                          ncol = p * k,
                          byrow = T)

  pca_sym <- prcomp(mat_symshapes)
  # Make PCA representing symmetrical variation (shapes and reflections averaged)
  
  mat_all <- matrix(data = orp(pAc$rotated), 
                    nrow = dim(pAc$rotated)[3], 
                    ncol = p * k,
                    byrow = T)
  
  mat_avg_rep <-  apply(mat_all, 
                        2, 
                        tapply, 
                        mirror.fac:indiv.fac.mirror, 
                        mean)
  # Matrix with replicate-average shapes 
  
  iMIR <- grep("MIRROR",
               rownames(mat_avg_rep))
  
  iORIG <- grep("ORIGINAL", 
                rownames(mat_avg_rep))
  
  mat_asym <- mat_avg_rep[iORIG, ] - mat_avg_rep[iMIR, ]
  # Matrix of asymmetric differences (averaged replicates)
  
  mat_asym2 <-  mat_all[grep("ORIG", 
                             mirror.fac), ] - mat_all[grep("MIRROR", 
                                                           mirror.fac), ]
  # Matrix of asymmetric differences (replicates separate)
  
  pca_asym <- prcomp(mat_asym, 
                     center = F, 
                     scale. = F)
  # PCA of asymmetric differences
  
  pcs <- pca_asym$x
  # PC scores
  
  eig <- pca_asym$sdev^2
  # eigenvalues
  
  preds <- predict(pca_asym, 
                   newdata = mat_asym2)
  # Projects the asymmetric differences for replicates onto the PCA space of
  # asymmetric differences with averaged replicates
  
  avgs <- apply(preds,
                2,
                mean)
  # Compute the average PC score across asymmetric PCs for all replicates
  
  M <- matrix(NA, ncol = 4, nrow = N)
  # Prepare empty matrix for individual values of DA FA and error indices.
  
  colnames(M) <- c("FAi", "FAi2", "DAi", "Ei")
  
  M[, 1] <- apply(pcs[, -1], 
                  1, 
                  sum)
  
  M[, 2] <- apply(t(t(pcs[, -1]) - avgs[-1]), 
                  1, 
                  sum)
  
  M[, 3] <- pcs[, 1]
  
  for (i in 1:N) {
    
    av <- pcs[i, -1]
    
    index <- which(indiv.fac == levels(indiv.fac)[i])
    
    reps <- preds[index, -1]

    diffs <- abs(t(t(reps) - av))

    M[i, 4] <- sum(diffs) / length(diffs)
    
  }
  
  return(list(PCA.sym = pca_sym, 
              PCA.asym = pca_asym, 
              M = data.frame(M),
              matshp = mat_all,
              Procrustes = pAc,
              mirror.fac = mirror.fac,
              indiv.fac.mirror = indiv.fac.mirror))
  
}
