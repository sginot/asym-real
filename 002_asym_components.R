#-------------------------------------------------------------------------------
# Computing estimates of inter-individual variation (IV), total asymmetry (TA)
# directional asymmetry (DA), fluctuating asymmetry (FA), and error in a sample
# (here simulated, but applicable to real samples as well).

# Script 001 must have been run first

#-------------------------------------------------------------------------------
# First, at the individual coordinate level (ie univariate case)

# Define function to do the decomposition of asymmetry
asym.decomp.univ <- function(x = stop("No variable defined"), 
                             replic = stop("No replicate factor defined"),
                             side = stop("No side factor defined")) {
  
  avg_replic <- tapply(x, 
                       replic,
                       mean)
  
  IV <- var(avg_replic) # Inter-individual variance
  
  avg_side <- tapply(x, 
                     side, 
                     mean)
  
  DA <- abs(diff(avg_side)) # Directional asymmetry, ie right average - left average
  
  avg_repside <- tapply(x, 
                        replic:side, 
                        mean)
  
  m_ar <- matrix(avg_repside, 
                 ncol = 2, 
                 byrow = T)
  
  TA <-  abs(apply(X = m_ar,
                   FUN = diff,
                   MARGIN = 1)) # Total asymmetry per individual
  
  FA <- abs(TA - DA) # Individual FA per individual
  
  avg_FA <- mean(FA) # Average FA
  
  error <- tapply(x, 
                  replic:side, 
                  diff) # Individual differences between replicates = error
  
  avg_error <- mean(abs(error))
  
  return(list(IV = IV,
              DA = DA,
              avg_FA = avg_FA,
              avg_error = avg_error,
              TA = TA,
              FA = FA,
              error = error))
}

#-------------------------------------------------------------------------------
# Wrapping function to apply asym.decomp.univ to all coordinates

asym.decomp.all <- function(mat = stop("A matrix containing coordinates must be defined"),
                            replic = stop("No replicate factor defined"),
                            side = stop("No side factor defined")) {
  
  mat_res <- matrix(NA,
                    ncol = 4,
                    nrow = ncol(mat))
  colnames(mat_res) <- c("IV", "DA", "FA", "error")
  
  for (i in 1:ncol(mat)) {
    
    res <- asym.decomp.univ(mat[,i],
                            replic = replic,
                            side = side)
    
    mat_res[i, ] <- unlist(res[1:4])
    
  }
  
  return(mat_res)
  
}

#-------------------------------------------------------------------------------
# Define function to do the decomposition of asymmetry for individual landmarks
# The difference with univariate case is that DISTANCES rather than differences
# must be used.
# Should apply to 2D landmarks
# Warning: do not use the function for the entire shapes, as it will only use
# the first 2 columns


asym.decomp.LM <- function(x = stop("No coordinates matrix defined"),
                             replic = stop("No replicate factor defined"),
                             side = stop("No side factor defined")) {
  
  avg_replic <- cbind(tapply(x[, 1], 
                             replic,
                             mean),
                      tapply(x[, 2], 
                             replic,
                             mean))
  
  IV <- var(avg_replic) # Inter-individual variance-covariance matrix
  
  avg_side <- cbind(tapply(x[, 1], 
                           side,
                           mean),
                    tapply(x[, 2], 
                           side,
                           mean))
  
  DA <- dist(avg_side) # Directional asymmetry, ie dist(right avg, left avg)
  
  avg_repside <- cbind(tapply(x[, 1], 
                              replic:side,
                              mean),
                       tapply(x[, 2], 
                              replic:side,
                              mean))
  
  
  TA <- rep(NA, nrow(avg_repside)/2)
  
  for (i in 1:nrow(avg_repside)/2) {
    
    TA[i] <- dist(rbind(avg_repside[i * 2, ],
                        avg_repside[(i * 2) - 1, ]))
    
  }
  # Total asymmetry per individual
  
  FA <- abs(TA - DA) # Individual FA per individual
  
  avg_FA <- mean(FA) # Average FA
  
  error <- rep(NA, nrow(x)/2)
  
  for (i in 1:nrow(x)/2) {
    
    error[i] <- dist(x[which(replic:side == levels(replic:side)[i]),])
    
  }
  # Individual differences between replicates = error
  
  avg_error <- mean(abs(error))
  
  return(list(IV = IV,
              DA = DA,
              avg_FA = avg_FA,
              avg_error = avg_error,
              TA = TA,
              FA = FA,
              error = error))
}


LM_DAs <- rep(NA, 6)
LM_FAs <- rep(NA, 6)
LM_error <- rep(NA, 6)

for (i in 1:6) {
  
  res <- asym.decomp.LM(x = mat_global[, c(i, i + 6)],
                        side = side_factor,
                        replic = replicat_factor)
  LM_DAs[i] <- res$DA
  LM_FAs[i] <- res$avg_FA
  LM_error[i] <- res$avg_error
}

#-------------------------------------------------------------------------------
# However, in general GMM datasets will be recorded as arrays of coordinates
# Therefore the following function will apply to an entire array, with user
# defined landmarks of interest,

# First problem is whether the array contains a structure that includes both
# sides or not (e.g. skull landmakrs vs hemi-mandible landmarks)

asym.decomp <- function(A = stop("No coordinates array defined"),
                        bilat = F,
                        replic = stop("No replicate factor defined"),
                        side = if (!bilat){stop("No side factor defined")},
                        symaxis = c(2, 8, 11)) {
  
  N <- dim(A)[3]
  P <- dim(A)[1]
  K <- dim(A)[2]
  
  if (length(replic) != N | length(side) != N) {
    stop("Replicate and side factors are not of same length as array")
  }
  
  print(c(N, P, K))
  
  M <- matrix(NA,
              ncol = P * K,
              nrow = N)
  
  for (i in 1:N) {M[i,] <- c(A[, , i])}
 
  avg_shp <- apply(M, 
                   2, 
                   mean)
  
  avg_sides <- apply(M, 
                     2, 
                     tapply, 
                     side, 
                     mean)
  
  avg_L <- avg_sides[1,]
  avg_R <- avg_sides[2,]
  
  avg_ind <- apply(M, 
                   2, 
                   tapply, 
                   replic, 
                   mean)
  
  avg_repside <- apply(M, 
                       2, 
                       tapply, 
                       replic:side, 
                       mean)
  
  dist_shp_repside <- rep(NA, N/4)
  
  for (i in 1 : (N/4)) {
    dist_shp_repside[i] <- dist(avg_repside[(i * 2 - 1) : (i * 2), ])
  }
  #dist between side averages should include only DA because FA and error should
  #have mean = 0.
  #However the mean distance between sides within individuals includes FA, there
  #fore should be higher than dist of side averages
  
  dist_shp_rep <- rep(NA, N/2)
  rs <- replic:side
  
  for (i in 1 : (N/2)) {
    dist_shp_rep[i] <- dist(M[which(rs == rs[i]), ])
  }
  
}

