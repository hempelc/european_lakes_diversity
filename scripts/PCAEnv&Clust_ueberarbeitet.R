# PCA and ENVIRONMENT
# V3.2 Author M.Jensen, May 12, 2016
# Biodiversity Group University of Duisburg-Essen, Germany
# based on Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
#**************************************************************************************


# chunk0: default filenames etc
library(vegan)
library(gclus)
library(dplyr)
library(data.table)
workdir <- "M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen"
envname <- "lakes_abiotic_data_just_numbers_to_work_with.csv" # name of environmental file
# end chunk 0 (files must be in the working directory)


# chunk1: some variables
# flags
minClust <- 1 # if False, the large cluster no is chosen according to both optimum criteria
winsys <- T # in case of Linux please set to FALSE
dendro <- T # if FALSE, no integrated dendrogram within PCA plot
log_env <- T
# variables for plotting
legend_pos <- "topright"
# symbol colors:
colorvec <- c("blue","red","green","brown","cyan",
              "darkmagenta","orange","forestgreen","gold","steelblue1","gray55","hotpink",
              "darkorange1","darkgreen","brown3","pink","darkolivegreen1",
              "lightcyan","aquamarine3","skyblue","lightblue","lightgrey",
              "magenta","pink2","lightyellow",
              "darkgrey","beige","white","turquoise","gray","purple","bisque")
pchvec <- rep(c(22,21,24,23,25),10)      # symbol type
envexists <- T # if FALSE no envfit takes place
decpoint <- ","  # German excel files with decimal comma, switch to "." in data.entry if needed
hoehe_einbeziehen <- T
axes <- c(1,2) # axes to be shown in the plot, could be also 3,2

# end variables



# chunk2: data entry: modification of default values
# first: file-names and work-directory
data.entry(envname,workdir,decpoint,log_env,minClust, hoehe_einbeziehen, legend_pos,winsys)
# second: variables
#if (envexists) dendro <- F
# end data entry



# chunk3: read environment data
setwd(workdir)
TI <- Sys.time()
env <- fread(envname,dec=decpoint, data.table=F)
rownames(env) <- env[,1]
env <- env[,-1]
colnames(env)[c(1:4, 34)]<- c("Altitude [m]", "pH", "Temp", "Cond", "GDP [mio]")
if (hoehe_einbeziehen==F) {env <- env[,-1]}
#env <- env[,c(1:3, 5, 15, 21, 23, 34)]
# end chunk3 read data


# chunk 4: log and select data; to be modified!!! no missing values allowed
if (log_env)
{ if (hoehe_einbeziehen==T) {nologs <- c("Altitude [m]", "pH", "Temp")} else {nologs <- c("pH", "Temp")}
cat("entry of env no_log_variables: no logarithm\n")
data.entry(nologs)
envnolog <- env[,nologs]
inds <- match(nologs,colnames(env),nomatch=0)
envlog <- env[, -inds]
envlog <- log1p(envlog)
env   <- cbind(envnolog,envlog)
#env <- select(env,Temp,pH,Cond,TP)
}  
# end log env data


# chunk 5: calculate PCA
env.pca <- rda(env, scale=TRUE)
# end chunk 5



# chunk 6 plot eigenvalues and % of variance for each axis
source("evplot (für Skript PCAEnv&Clust und PCAOTU_EnvPosthoc&Clust).R")
ev <- env.pca$CA$eig
if (winsys) windows(title="PCA eigenvalues")
evplot(ev)
cat("Percentage/PCA-Axis","\n")
print(100*ev/sum(ev))
# end chunk 6




# chunk7: Ward clustering of dissimilarity matrix and extraction of groups
dist_env <- dist(scale(env))
env.ward <- hclust(dist_env, "ward.D2")
# end


# chunk8: optimal number of clusters according to silhouette widths
# (Rousseeuw quality index)
# First, create an empty vector in which the asw values will be written
asw <- numeric(nrow(env))
for (k in 2:(nrow(env)-1)) {
  sil <- silhouette(cutree(env.ward, k=k), dist_env)
  asw[k] <- summary(sil)$avg.width
}
k.best <- which.max(asw)
# end chunk 8 optimal number



# chunk9: Plot average silhouette widths (using Ward clustering) for all partitions 
# except for the trivial partition in a single group (k=1)
if (winsys) windows(title="Numbers: Silhouettes - Ward - k = 2 to n-1",12,12)
plot(1:nrow(env), asw, type="h", 
     main=paste("Silhouette-optimal number of clusters, Ward/"), 
     xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")
# end plot



# chunk10: optimal number of clusters according to Mantel statistic (Pearson)
# Function to compute a binary distance matrix from groups
grpdist <- function(X)
{ require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  return(distgr) }
# Run based on the Ward clustering
kt <- data.frame(k=1:nrow(env), r=0)
for (i in 2:(nrow(env)-1)) {
  gr <- cutree(env.ward, i)
  distgr <- grpdist(gr)
  mt <- cor(dist_env, distgr, method="pearson")
  kt[i,2] <- mt
}
kt
k.bestMant <- which.max(kt$r)
# end chunk10 optimal number: loop could be optimized by e.g. compiler




# chunk11: plot Mantel-optimal number
if (winsys) windows(title="Optimal number of clusters - Mantel",12,12)
plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - Ward", 
     xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.bestMant, paste("optimum", k.bestMant, sep="\n"), col="red", font=2,
     col.axis="red")
points(k.bestMant, max(kt$r), pch=16, col="red", cex=1.5)
cat("", "Mantel-optimal number of clusters k =", k.bestMant, "\n", 
    "with a matrix linear correlation of", max(kt$r), "\n")
# end chunk11



# chunk12: decision about optimal number of clusters
if (minClust==1) k.best <- min(k.best,k.bestMant) else k.best <- max(k.best,k.bestMant)
if (minClust>1) k.best <- minClust # manual choice of cluster number
# end



# chunk13: dendrogramm
source("hcoplot (für Skript NMDSClust und PCAEnv&Clust).R")
if (winsys) windows(title="Dendrogramm - Ward",16,12)
par(mar=c(7,4,4,6)+0.1)
hcoplot(env.ward, dist_env, k=k.best)
# end dendrogramm




# chunk14: Silhouette plot
cutg <- cutree(env.ward, k=k.best)
sil <- silhouette(cutg, dist_env)
rownames(sil) <- row.names(env)
maxlen <- max(length(rownames(sil)))
if (winsys) windows(title="Silhouette plot - Ward",13,12)
if (winsys) par(pin=c(7,7))
maintit <- paste("Silhouette Plot,",envname)
lettersize <-  50 / nrow(env)/2
plot(sil, max.strlen = maxlen, main=maintit,
     cex.names=lettersize, col=colorvec[1:k.best], nmax=250)
# end chunk18 Silhouette plot



# chunk 19: cleanplot (biplot with env_data)
if (winsys) windows(title="PCA", 12, 8)
source("cleanplot3.pca (für Skript PCAEnv&Clust und PCAOTU_EnvPosthoc&Clust).R")
gr <- cutree(env.ward, k=k.best)
grl <- levels(factor(gr))
cleanplot.pca(env.pca, ax1=1, ax2=2, point=T, ahead=0)
legend(legend_pos, paste("Cluster", c(1:length(grl))), pch=pchvec[1:length(grl)], 
       pt.bg=colorvec[c(1:length(grl))], pt.cex=1.2,cex=0.85)
# end chunk 19 : cleanplot with env_data



# chunk 20: plot the sites with cluster symbols and colours (scaling 1)
if (winsys) windows(title="Ordination and clustering",16,16)
sit.sc1 <- scores(env.pca, display="wa", scaling=1)
pcatitle <- paste("PCA + WARD-Clust scaling 1;",envname)
p <- plot(env.pca, choices=axes, display="wa", scaling=1, type="n", main=pcatitle)
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
env.bw.groups <- cutree(env.ward, k=k.best)
grp.lev <- levels(factor(env.bw.groups))
for (i in 1:length(grp.lev))
{ points(sit.sc1[env.bw.groups==i,axes[1]],sit.sc1[env.bw.groups==i,axes[2]], 
         pch=pchvec[i], cex=1.2, col="black",bg=colorvec[i]) }
text(sit.sc1[ ,axes[1]], sit.sc1[ ,axes[2]],row.names(env), pos=4, cex=0.7, col="black")
# Add the dendrogram if dendro=TRUE
if (dendro) ordicluster(p, env.ward, col="dark grey")
legend(legend_pos, paste("Cluster", c(1:length(grp.lev))), pch=pchvec[1:length(grp.lev)], 
       pt.bg=colorvec[1:length(grp.lev)], pt.cex=1.2,col="black",cex=0.85)
# end chunk 20 plot the sites with cluster symbols



# end of script **********************************************************************************