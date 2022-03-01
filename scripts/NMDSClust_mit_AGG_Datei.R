# NMDS script by Manfred Jensen, V3.7, 1.6.2016, Biodiversity group, University Duisburg-Essen
# mainly based on Borcard,Gillet,Legendre 2011: Numerical Ecology with R (Springer:Berlin/Heidelberg)


# chunk0: default filenames etc
library(vegan)
#library(vegan3d) 3d plot possible with ordiplot3d, see line 244, but not very beautiful
library(gclus)
library(dplyr)
library(data.table)
library(ellipse)
workdir <- "M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen"
spefileprn <- "kingdom_Euk_without_Metazoa_and_Embryophyta_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2.agg"
envname <- "lakes_abiotic_data_just_numbers_to_work_with_log.csv" # name of environmental file
# end chunk 0 (files must be in the working directory)


# chunk1: some variables
# flags
rrare <- T  # rarefaction without info loss
rrare <- F  # no classical rarefaction
hel   <- T  # default: Hellinger/euclidic distance, if FALSE: Bray-Curtis distance
pa_flag <- F # presence-absence flag: all data are transformed to 0,1 and no rarefaction
minClust <- F # if False, the large cluster no is chosen according to both optimum criteria
winsys <- T # in case of Linux please set to FALSE
rlev_man <- 0
# variables for plotting
legend_later <- F
legend_pos <- "topright"
clust_ellip <- F
ellip_rad <- 0.95
dendro <- T # if FALSE, no integrated dendrogram within NMDS plot
x_orient <- T
y_orient <- T
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
# end variables



# chunk2: data entry: modification of default values
# first: file-names and work-directory
data.entry(spefileprn,workdir)
data.entry(envexists,envname,decpoint)
# second: variables
data.entry(pa_flag,hel,rlev_man,minClust,dim_num,axes,x_orient,y_orient)
data.entry(legend_pos,legend_later,dendro,clust_ellip,ellip_rad,winsys)
# end data entry


# chunk3: read abundance data
setwd(workdir)
TI <- Sys.time()
spe <- fread(spefileprn,dec=decpoint, data.table=F)
lakenames <- spe[,1]
spe <- as.data.frame(t(apply(spe[,-1], 2, FUN=function (x) as.numeric(as.character(x)))))
colnames(spe) <- lakenames
# end read data



# chunk4: data structure: how many OTUs only in 1,2..10 lakes
spepa <- decostand(spe,"pa")
spepasum <- rowSums(spepa)      # number of sites in which the OTUs (rows) occur...
sitecounts <- function(i,spepasum){sum(as.integer(spepasum==i))}
sitecounts10 <- sapply(1:10,sitecounts,spepasum=spepasum)
cat("community structure:\nhow many OTUs restricted to 1..10 sites: ",sitecounts10,"\n")
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
# further processing with sperr
#left out because .agg is already rarefied
sperr <- spe


# chunk7: generate dissimilarity matrix based on hellinger.bray_curtis.ochiai.sÃ¶rensen distances
if (hel) 
{spe.h <- decostand(sperr,"hel")
if (pa_flag) dist_text <- "Ochiai" else dist_text <- "Hel" # Hellinger => Ochiai for pa
dist_method <- "euc"
}else 
{if (pa_flag) spe.h <- sperr else 
  spe.h <- decostand(sperr,"total")
if (pa_flag) dist_text <- "SÃ¶rensen" else dist_text <- "Bray" # Bray-Curtis => SÃ¶rensen for pa
dist_method <- "bray"
}
spe.dist <- vegdist(spe.h,dist_method)
# end 


# ******************************************************************************************
# chunk8: calculate nmds and the stress variable as as measure of goodness
spe.nmds <- metaMDS(spe.h, distance=dist_method, k=dim_num, autotransform = FALSE,trymax=50)
stress <- round(spe.nmds$stress,digits=3)
cat("\n Stress = ",spe.nmds$stress,"\n")
#cat(" rarefaction level = ",rlev,"\n")
# end calculation
# ******************************************************************************************


# chunk9 (not active): plot nmds without cluster analysis
# if (winsys) windows(title=paste("NMDS on OTUs - ",dist_text))
# plot(spe.nmds, type="t", main=paste(dist_text,"NMDS - Stress =",round(spe.nmds$stress,3)))
# end plot



# chunk10: Shepard plot and goodness of fit
if (winsys) windows(title="NMDS - Shepard plot", 16, 9)
par(mfrow=c(1,2))
stressplot(spe.nmds, main="Shepard plot")
gof <- goodness(spe.nmds)
#plot(spe.nmds, type="t", main="Goodness of fit")
plot(spe.nmds, , main="Goodness of fit", display=c("sites"), cex=gof*300)



# Add colours from WARD clustering to an NMDS plot


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
plot(1:nrow(spe), asw, type="h", 
     main=paste("Silhouette-optimal number of clusters, Ward/",dist_text), 
     xlab="k (number of groups)", ylab="Average silhouette width")
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
par(mar=c(7,4,4,8)+0.1)
hcoplot(spe.dist.ward, spe.dist, k=k.best)
# end dendrogramm


# chunk18: Silhouette plot
cutg <- cutree(spe.dist.ward, k=k.best)
sil <- silhouette(cutg, spe.dist)
rownames(sil) <- row.names(spe)
maxlen <- max(length(rownames(sil)))
if (winsys) windows(title="Silhouette plot - Ward",13,12)
if (winsys) par(pin=c(7,7))
maintit <- paste("Silhouette Plot,",dist_text)#,spefileprn)
lettersize <-  50 / nrow(spe)/2
plot(sil, max.strlen = maxlen, main=maintit,
     cex.names=lettersize, col=colorvec[1:k.best], nmax=250)
# end chunk18 Silhouette plot




# chunk19: Combination of Clustering and NMDS result
spe.bw.groups <- cutree(spe.dist.ward, k=k.best)
grp.lev <- levels(factor(spe.bw.groups))
sit.sc <- scores(spe.nmds)
if (winsys) windows(title="NMDS plot with cluster colors",12,12)
titletext <- paste("NMDS(stress:",stress,", dim:",dim_num,") + Ward-clusters "," (",dist_text," dist.)",sep="")  #spefileprn,
xrange <- range(sit.sc[,axes[1]])
yrange <- range(sit.sc[,axes[2]])
if (x_orient==F) xrange <- rev(xrange)
if (y_orient==F) yrange <- rev(yrange)
p <- ordiplot(sit.sc, choices=axes,type="none",main=titletext,xlim=xrange,ylim=yrange)
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
for (i in 1:length(grp.lev))
{ points(sit.sc[spe.bw.groups==i,axes[1]],sit.sc[spe.bw.groups==i,axes[2]], 
         pch=pchvec[i], cex=1.2, col="black",bg=colorvec[i]) }
text(sit.sc[ ,axes[1]], sit.sc[ ,axes[2]],row.names(spe), pos=4, cex=0.7, col="black")
# Add the dendrogram if dendro=TRUE
if (dendro) ordicluster(p, spe.dist.ward, col="dark grey")
if (legend_later==F)
  legend(legend_pos, paste("Cluster", c(1:length(grp.lev))), pch=pchvec[1:length(grp.lev)], 
         pt.bg=colorvec[1:length(grp.lev)], pt.cex=1.2,col="black",cex=0.85)
# end chunk19 plot results
#p <- ordiplot3d(sit.sc,display="sites", choices=1:3, main=titletext) (3d possible)

# chunk19a: draw ellipses around clusters
if (clust_ellip)
{sit_2d <- cbind(sit.sc[,axes[1]],sit.sc[,axes[2]])
for(i in 1:length(grp.lev))
{	sit_2di <- sit_2d[spe.bw.groups==i,]
cov <- cov(sit_2di)
centre <- apply(sit_2di, 2, mean)
lines(ellipse(cov, centre=centre, level=ellip_rad),col=colorvec[i]) }}
# end chunk19a


# chunk20: a posteriori projection of environmental variables in NMDS (caution)
# log nolog etc have to be modified depending on the env-file
if (envexists)
{ # subchunk 1: read env data
  env <- fread(envname,dec=decpoint, data.table=F)
  row.names(env) <- env[,1]
  env <- select(env, -Seekürzel)
  # end subchunk 1  
  
  
  # subchunk3: adjust env data (only sites occurring in spe)
  source("adjustenv (für Skript NMDSClust und PCAOTU_EnvPosthoc&Clust).R")
  env <- adjustenv(spe,env)
  #env$spec_rich <- rowSums(sperr)
  env <- scale(env) # z-transformations are necessary here ...
  # end subchunk 3
  
  # subchunk 4: a posteriori (passive) fit env data to spe ordination
  if (dim_num==2) 
    spe.nmds.env <- envfit(spe.nmds, env, choices=axes, na.rm=T, permutations = 11000) else
      spe.nmds.env <- envfit(spe.nmds, env, na.rm=T, permutations = 11000)   
  print(spe.nmds.env)
  plot(spe.nmds.env,col="cyan")
  # Plot significant variables with a different colour
  plot(spe.nmds.env, p.max = 0.001, col="blue")
  # end subchunk 4
}#end chunk20 env posthoc

cat("runtime: ")
print(Sys.time() - TI)

if (legend_later)
  legend(locator(1), paste("Cluster", c(1:length(grp.lev))), pch=pchvec[1:length(grp.lev)], 
         pt.bg=colorvec[1:length(grp.lev)], pt.cex=1.2,col="black",cex=0.85)

#******************** end of script ******************************************************