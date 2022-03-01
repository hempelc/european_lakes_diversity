################################################################################
### CHAPTER 3: ASSOCIATION MEASURES
### Updated by F. Gillet on 25.08.2012
###
### Online supporting material for: 
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
################################################################################

# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(dplyr)
library(data.table)


# chunk1 variables 
workdir <- "M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen"
filename <- "kingdom_Euk_without_Metazoa_and_Embryophyta_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2"
winsys <- T
# end variables


# chunk2 data entry
data.entry(filename,workdir,winsys)
setwd(workdir)
source("coldiss (für Skript heatmaps_lakes).R")
# end data entry


# chunk3: read data
spefileprn <- paste(filename,".csv",sep="")
spe <- fread(spefileprn,dec=",", data.table=F)
spe <- as.data.frame(apply(spe, 2, FUN=function (x) as.numeric(as.character(x))))
# end read data



# Q-mode dissimilarity and distance measures for (semi-)quantitative data
# ***********************************************************************
# Bray-Curtis dissimilarity matrix on raw species data
spe.db <- vegdist(spe)	# Bray-Curtis dissimilarity (default)
head(spe.db)
# Bray-Curtis dissimilarity matrix on log-transformed abundances
spe.dbln <- vegdist(log1p(spe))
head(spe.dbln)
# Chord distance matrix
spe.norm <- decostand(spe, "nor")
spe.dc <- dist(spe.norm)
head(spe.dc)
# Hellinger distance matrix
spe.hel <- decostand(spe, "hel")
spe.dh <- dist(spe.hel)
head(spe.dh)


# Q-mode dissimilarity measures for binary data
# *********************************************

# Jaccard dissimilarity matrix using function vegdist()
spe.dj <- vegdist(spe, "jac", binary=TRUE)
head(spe.dj)
head(sqrt(spe.dj))
# Jaccard dissimilarity matrix using function dist()
spe.dj2 <- dist(spe, "binary")
head(spe.dj2)
# Jaccard dissimilarity matrix using function dist.binary()
spe.dj3 <- dist.binary(spe, method=1)
head(spe.dj3)
# Sorensen dissimilarity matrix using function dist.binary()
spe.ds <- dist.binary(spe, method=5)
head(spe.ds)
# Sorensen dissimilarity matrix using function vegdist()
spe.ds2 <- vegdist(spe, binary=TRUE)
head(spe.ds2)
head(sqrt(spe.ds2))
# Ochiai dissimilarity matrix
spe.och <- dist.binary(spe, method=7)
head(spe.och)


# Graphical display of association matrices
# Colour plots (also called heat maps, or trellis diagrams in the data 
# analysis literature) using the coldiss() function
# ********************************************************************

# Usage:
# coldiss(D = dissimilarity.matrix, nc = 4, byrank = TRUE, diag = FALSE)
# If D is not a dissimilarity matrix (max(D) > 1), then D is divided by max(D)
# nc 							number of colours (classes)
# byrank= TRUE		equal-sized classes
# byrank= FALSE		equal-length intervals
# diag = TRUE			print object labels also on the diagonal

# Compare dissimilarity and distance matrices obtained from the species data
# 4 colours with equal-length intervals
# --------------------------------------------------------------------------

# Bray-Curtis dissimilarity matrix on raw species abundance data
if (winsys) windows(title="Bray-Curtis (raw data)", 10, 5)
coldiss(spe.db, byrank=FALSE, diag=TRUE)

# Same but on log-transformed data
if (winsys) windows(title="Bray-Curtis [ln(y+1) data]", 10, 5)
coldiss(spe.dbln, byrank=FALSE, diag=TRUE)

# Chord distance matrix
if (winsys) windows(title="Chord", 10, 5)
coldiss(spe.dc, byrank=FALSE, diag=TRUE)

# Hellinger distance matrix
if (winsys) windows(title="Hellinger", 10, 5)
coldiss(spe.dh, byrank=FALSE, diag=TRUE)

# Jaccard distance matrix
if (winsys) windows(title="Jaccard", 10, 5)
coldiss(spe.dj, byrank=FALSE, diag=TRUE)

# Simple matching dissimilarity
# (called the Sokal and Michener index in ade4)
spe.s1 <- dist.binary(spe, method=2)
if (winsys) windows(title="S1 on species data", 10, 5) 
coldiss(spe.s1^2, byrank=FALSE, diag=TRUE)


# Compare distance matrices from environmental, species and spatial data
# 16 colours with equal-size classes
# ----------------------------------------------------------------------

# Hellinger distance matrix of the species data (equal-sized classes)
if (winsys) windows(title="Species Hellinger equal-sized", 10, 5)
coldiss(spe.dh, nc=16, diag=TRUE)







# R-mode dissimilarity matrices
# *****************************

# Transpose matrix of species abundances
spe.t <- t(spe)

# Chi-square pre-transformation followed by Euclidean distance
  #spe.t.chi <- decostand(spe.t, "chi.square")
  #spe.t.D16 <- dist(spe.t.chi)
  #if (winsys) windows(title="D16 on fish species (R-mode)", 10, 5)
  #coldiss(spe.t.D16, diag=F)

# Jaccard index on fish presence-absence
spe.t.S7 <- vegdist(spe.t, "jaccard", binary=TRUE)
  #if (winsys) windows(title="S7 on fish species (R-mode)", 10, 5) 
  #coldiss(spe.t.S7, diag=F)


