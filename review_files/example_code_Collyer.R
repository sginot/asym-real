library(geomorph)
data(pupfish)

head <- c(4, 10:17, 39:56) # head landmarks
tail <- setdiff(1:56, head) # tail landmarks

# First, performing modularity and integration tests, using the full configuration.

partition <- rep(0, 56)
partition[tail] <- 1
MT <- modularity.test(pupfish$partition.gp = partition)
MT

IT <- integration.test(pupfish$partition.gp = partition)
IT

# Next, perform Procrustes superimposition on modules, separately.  Additionally, rotate the Procrustes coordinates from one module 90 degrees, just to show how arbitrary rotation can affect results.

GPA.h1 <- gpagen(pupfish$coords[head, , ]) # just head
GPA.h2 <- GPA.h1
rot.angle <- 90 # degrees
rot.rad <- rot.angle * pi / 180
rot.mat <- matrix(c(cos(rot.rad), sin(rot.rad), cos(rot.rad), -sin(rot.rad)), 2, 2)
GPA.h2$coords <- simplify2array(lapply(1:54, function(j) GPA.h2$coord[,,j] %*% rot.mat)) # just head landmarks, rotated 90 degrees
GPA.t <- gpagen(pupfish$coords[tail, , ]) # just tail

# These can be combined and analyses run.

coords1 <- combine.subsets(GPA.h1, GPA.t, gpa = FALSE)$coords
coords2 <- combine.subsets(GPA.h2, GPA.t, gpa = FALSE)$coords
gp <- factor(c(rep(1, 27), rep(2, 29)))
MT1 <- modularity.test(coords1, partition.gp = gp)
MT2 <- modularity.test(coords2, partition.gp = gp)
MT1

MT2

IT1 <- integration.test(coords1, partition.gp = gp)
IT2 <- integration.test(coords2, partition.gp = gp)
IT1

IT2

# So, as it can be seen, the arbitrary rotation of modules will change the results. As the authors noted, there is a problem of scaling, because each configuration is transformed to unit size. This can be avoided by combining landmark configurations and scaling them relative to their portions of the original centroid size.

coords3 <- combine.subsets(GPA.h1, GPA.t, gpa = FALSE, CS.sets = list(GPA.h1$Csize, GPA.t$Csize))$coords
coords4 <- combine.subsets(GPA.h2, GPA.t, gpa = FALSE, CS.sets = list(GPA.h2$Csize, GPA.t$Csize))$coords
MT3 <- modularity.test(coords3, partition.gp = gp)
MT4 <- modularity.test(coords4, partition.gp = gp)
MT3

MT4

IT3 <- integration.test(coords3,
partition.gp = gp)
IT4 <- integration.test(coords4, partition.gp = gp)
IT3

IT4

