# PCA and ENVIRONMENT on rarefied data + Hellinger
# V3.3 Author M.Jensen, May 30, 2016
# Biodiversity Group University of Duisburg-Essen, Germany
# based on Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
#**************************************************************************************


# chunk0: default filenames etc
library(vegan)
#library(vegan3d) 3d plot possible with ordiplot3d, see line 244, but not very beautiful
library(gclus)
library(dplyr)
library(data.table)
library(ellipse)
workdir <- "M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen"
spefileprn <- "kingdom_Euk_without_Metazoa_and_Embryophyta_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2.agg"
envname <- "lakes_abiotic_data_just_numbers_to_work_with.csv" # name of environmental file
# end chunk 0 (files must be in the working directory)


# chunk1: some variables
# flags
#rrare <- T  # rarefaction without info loss
hel   <- T  # default: Hellinger/euclidic distance, if FALSE: Bray-Curtis distance
pa_flag <- F # presence-absence flag: all data are transformed to 0,1 and no rarefaction
minClust <- 1 # if False, the large cluster no is chosen according to both optimum criteria
winsys <- T # in case of Linux please set to FALSE
dendro <- T # if FALSE, no integrated dendrogram within NMDS plot
rlev_man <- 0
# variables for plotting
legend_pos <- "topright"
clust_ellip <- F
ellip_rad <- 0.95
# symbol colors:
colorvec <- c("blue","red","green","brown","cyan",
              "darkmagenta","orange","forestgreen","gold","steelblue1","gray55","hotpink",
              "darkorange1","darkgreen","brown3","pink","darkolivegreen1",
              "lightcyan","aquamarine3","skyblue","lightblue","lightgrey",
              "magenta","pink2","lightyellow",
              "darkgrey","beige","white","turquoise","gray","purple","bisque")
pchvec <- rep(c(22,21,24,23,25),10)      # symbol type
dim_num <- 2   # No of dimensions (axes) for calculating NMDS
axes <- c(1,2) # axes to be shown in the plot, could be also 3,2
envexists <- T # if FALSE no envfit takes place
decpoint <- ","  # German excel files with decimal comma, switch to "." in data.entry if needed
y_orient <- F
# end variables



# chunk2: data entry: modification of default values
# first: file-names and work-directory
data.entry(spefileprn,workdir,envexists,envname,decpoint)
# second: variables
#if (envexists) dendro <- F
data.entry(pa_flag,rlev_man,minClust)
data.entry(legend_pos,dendro,clust_ellip,ellip_rad,winsys)
# end data entry



# chunk3: read abundance data
setwd(workdir)
TI <- Sys.time()
spe <- fread(spefileprn,dec=decpoint, data.table=F)
lakenames <- spe[,1]
spe <- spe[,-1]
spe <- as.data.frame(apply(spe, 2, FUN=function (x) as.numeric(as.character(x))))
rownames(spe) <- lakenames
spe <- t(spe)
# end read data


# chunk4: data structure: how many OTUs only in 1,2..10 lakes
spepa <- decostand(spe,"pa")
spepasum <- rowSums(spepa)      # number of sites in which the OTUs (rows) occur...
sitecounts <- function(i,spepasum){sum(as.integer(spepasum==i))}
sitecounts10 <- sapply(1:10,sitecounts,spepasum=spepasum)
cat("community structure:\nhow many Groups restricted to 1..10 sites: ",sitecounts10,"\n")
# end structure info


# chunk5: transpose data (for vegan: samples in rows and OTUs in columns) and pa transform if needed
spe <- t(spe)
cat("number of sites: ",nrow(spe),"\n")
if (pa_flag) spe <- decostand(spe,"pa")
# now vegan functions can be used



# chunk6: rarefaction of sites with max possible level
#source("rrare100.R")
#if (rlev_man==0) rlev <- min(rowSums(spe)) else rlev <- rlev_man
#if (pa_flag==F & rrare) sperr <- sperrare <- rrarefy100(round(spe, digits=0),rlev) else sperr <- spe
# end chunk 6 further processing with sperr
# .agg data is already rarefied
sperr <- spe



# chunk7: generate dissimilarity matrix based on hellinger.bray_curtis.ochiai.sÃ¶rensen distances
if (hel) 
{spe.h <- decostand(sperr,"hel")
if (pa_flag) dist_text <- "Ochiai" else dist_text <- "Hel" # Hellinger => Ochiai for pa
dist_method <- "euc"
}else 
{if (pa_flag) spe.h <- sperr else 
  spe.h <- decostand(sperr,"total")
if (pa_flag) dist_text <- "Soerensen" else dist_text <- "Bray" # Bray-Curtis => Soerensen for pa
dist_method <- "bray"
}
spe.dist <- vegdist(spe.h,dist_method)
# end chunk7


# chunk11: Ward clustering of dissimilarity matrix and extraction of groups
spe.dist.ward <- hclust(spe.dist, "ward.D2")
# end


# chunk12: optimal number of clusters according to silhouette widths
# (Rousseeuw quality index)
# First, create an empty vector in which the asw values will be written
asw <- numeric(nrow(spe))
for (k in 2:(nrow(spe)-1)) {
  sil <- silhouette(cutree(spe.dist.ward, k=k), spe.dist)
  asw[k] <- summary(sil)$avg.width
}
k.best <- which.max(asw)
# end optimal number




# chunk13: Plot average silhouette widths (using Ward clustering) for all partitions 
# except for the trivial partition in a single group (k=1)
if (winsys) windows(title="Numbers: Silhouettes - Ward - k = 2 to n-1",12,12)
plot(1:nrow(spe), asw, type="h", main=paste("Silhouette-optimal number of clusters, Ward/",dist_text), xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")
# end plot



# chunk14: optimal number of clusters according to Mantel statistic (Pearson)
# Function to compute a binary distance matrix from groups
grpdist <- function(X)
{ require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  return(distgr) }
# Run based on the Ward clustering
kt <- data.frame(k=1:nrow(spe), r=0)
for (i in 2:(nrow(spe)-1)) {
  gr <- cutree(spe.dist.ward, i)
  distgr <- grpdist(gr)
  mt <- cor(spe.dist, distgr, method="pearson")
  kt[i,2] <- mt
}
kt
k.bestMant <- which.max(kt$r)
# end chunk 14 optimal number: loop could be optimized by e.g. compiler




# chunk15: plot Mantel-optimal number
if (winsys) windows(title="Optimal number of clusters - Mantel",12,12)
plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - Ward", 
     xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.bestMant, paste("optimum", k.bestMant, sep="\n"), col="red", font=2,
     col.axis="red")
points(k.bestMant, max(kt$r), pch=16, col="red", cex=1.5)
cat("", "Mantel-optimal number of clusters k =", k.bestMant, "\n", 
    "with a matrix linear correlation of", max(kt$r), "\n")
# end chunk15



# chunk16: decision about optimal number of clusters
if (minClust==1) k.best <- min(k.best,k.bestMant) else k.best <- max(k.best,k.bestMant)
if (minClust>1) k.best <- minClust # manual choice of cluster number
# end



# chunk17: dendrogramm
source("hcoplot (für Skript NMDSClust und PCAEnv&Clust).R")
if (winsys) windows(title="Dendrogramm - Ward",16,12)
par(mar=c(7,4,4,6)+0.1)
hcoplot(spe.dist.ward, spe.dist, k=k.best)
# end dendrogramm




# chunk18: Silhouette plot
cutg <- cutree(spe.dist.ward, k=k.best)
sil <- silhouette(cutg, spe.dist)
rownames(sil) <- row.names(spe)
maxlen <- max(length(rownames(sil)))
if (winsys) windows(title="Silhouette plot - Ward",13,12)
if (winsys) par(pin=c(7,7))
maintit <- paste("Silhouette Plot,",dist_text) #,spefileprn)
lettersize <-  50 / nrow(spe)/2
plot(sil, max.strlen = maxlen, main=maintit,
     cex.names=lettersize, col=colorvec[1:k.best], nmax=250)
# end chunk18 Silhouette plot




#*************************
# chunk 21 calculate PCA
spe.h.pca <- rda(spe.h,k=dim_num)
#*************************


# chunk 22 plot eigenvalues and % of variance for each axis
source("evplot (für Skript PCAEnv&Clust und PCAOTU_EnvPosthoc&Clust).R")
ev <- spe.h.pca$CA$eig
if (winsys) windows(title="PCA eigenvalues")
evplot(ev)
cat("Percentage/PCA-Axis","\n")
print(100*ev/sum(ev))
# end chunk 22




# chunk 23: cleanplot (biplot with OTUs)
if (winsys) windows(title="PCA(HELLINGER) on rrarefied OTUs", 16, 11)
source("cleanplot3.pca (für Skript PCAEnv&Clust und PCAOTU_EnvPosthoc&Clust).R")
gr <- cutree(spe.dist.ward, k=k.best)
grl <- levels(factor(gr))
cleanplot.pca(spe.h.pca, ax1=1, ax2=2, point=TRUE, ahead=0)
legend(legend_pos, paste("Cluster", c(1:length(grl))), pch=pchvec[1:length(grl)], 
       pt.bg=colorvec[c(1:length(grl))], pt.cex=1.2,cex=0.85)
# end chunk 23: cleanplot with OTUs




# chunk 24: plot the sites with cluster symbols and colours (scaling 1)
if (winsys) windows(title="Ordination and clustering",16,16)
sit.sc1 <- scores(spe.h.pca, display="wa", scaling=1)
pcatitle <- paste("PCA(OTU-HELLINGER) + WARD-Clust scaling 1",
                  "; rlevel:",as.character(rlev)) #,spefileprn 
yrange <- range(sit.sc1[,2])
if (y_orient==F) yrange <- rev(yrange)
p <- plot(spe.h.pca, choices=axes, display="wa", scaling=1, type="n", main=pcatitle, ylim=yrange)

abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
spe.bw.groups <- cutree(spe.dist.ward, k=k.best)
grp.lev <- levels(factor(spe.bw.groups))
for (i in 1:length(grp.lev))
{ points(sit.sc1[spe.bw.groups==i,axes[1]],sit.sc1[spe.bw.groups==i,axes[2]], 
         pch=pchvec[i], cex=1.2, col="black",bg=colorvec[i]) }
text(sit.sc1[ ,axes[1]], sit.sc1[ ,axes[2]],row.names(spe), pos=4, cex=0.7, col="black")
# Add the dendrogram if dendro=TRUE
if (dendro) ordicluster(p, spe.dist.ward, col="dark grey")
legend(legend_pos, paste("Cluster", c(1:length(grp.lev))), pch=pchvec[1:length(grp.lev)], 
       pt.bg=colorvec[1:length(grp.lev)], pt.cex=1.2,col="black",cex=0.85)
# end chunk 24 plot the sites with cluster symbols


# chunk24a: draw ellipses around clusters
if (clust_ellip)
  for(i in 1:length(grp.lev))
  {  
  cov <- cov(sit.sc1[spe.bw.groups==i, ], centre <- apply(sit.sc1[spe.bw.groups==i, ], 2, mean), lines(ellipse(cov, centre=centre, level=ellip_rad),col=colorvec[i]))
  }
# end chunk24a




# chunk25: a posteriori projection of environmental variables in PCA
# log nolog etc have to be modified depending on the env-file
if (envexists)
{ env <- fread(envname,dec=decpoint, data.table=F)
row.names(env) <- env[,1]
env <- select(env, -Seekürzel)
source("adjustenv (für Skript NMDSClust und PCAOTU_EnvPosthoc&Clust).R")
env <- adjustenv(spe,env)
#if (pa_flag!=T) env$spec_rich <- rowSums(sperr)
env <- scale(env) # z-transformations are necessary here ...
spe.nmds.env <- envfit(spe.h.pca, env, choices=axes, na.rm=T, permutations=9999)
print(spe.nmds.env)
plot(spe.nmds.env,col="cyan")
# Plot significant variables with a different colour
plot(spe.nmds.env, p.max=0.001, col="blue")
} #end chunk25 env posthoc


cat("runtime: ")
print(Sys.time() - TI)

#******************** end of script ******************************************************