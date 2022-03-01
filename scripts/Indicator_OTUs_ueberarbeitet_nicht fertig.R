# Constrained Cluster Analysis and Indicator OTUs on rarefied data distance Hellinger or Bray ...
# V2.9 Author M.Jensen, June 29, 2016
# Biodiversity Group University of Duisburg-Essen, Germany
# based on Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
# indicspecies by Miquel De Caceres & Florian Jansen plus rioja by Stephen Juggins
#**************************************************************************************

# chunk0: default filenames etc
library(gclus)
library(vegan)
library(rioja)
library(indicspecies)
library(dplyr)
library(data.table)


workdir <- "M:/Allgemeine Daten/Statistik_Kolloq"
spefileind <- "Alpine_new.ind"
envname <- "AlpineSpatialEnv_LGMJ.csv" # name of environmental file
# end chunk 0 (files must be in the working directory)


# chunk1: some variables
# flags
drare <- T  # rarefaction without info loss
rrare <- F  # no classical rarefaction
hel   <- F  # default: Hellinger/euclidic distance, if FALSE: Bray-Curtis distance
pa_flag <- F # presence-absence flag: all data are transformed to 0,1 and no rarefaction
minClust <- 1 # if False, the large cluster no is chosen according to both optimum criteria
winsys <- T # in case of Linux please set to FALSE
dendro <- T # if FALSE, no integrated dendrogram within NMDS plot
print_divs <- F
rlev_man <- 0
add_spec <- F
indyesno <- 1
maxorder <- 3
pval_lev <- 0.05
permnum <- 999
txt_file <- F
cov_plot <- F
A_t <- 0.95
# variables for plotting
legend_later <- F
clust_ellip <- T
ellip_rad <- 0.95
legend_pos <- "bottomleft"
x_orient <- F
y_orient <- T
# symbol colors:
colorvec <- c("red","green","cyan","darkblue","khaki",
              "magenta","yellow","lightgreen","orange","steelblue1","grey","black",
              "darkorange1","lightyellow","brown","pink","darkolivegreen1",
              "lightcyan","aquamarine3","blue","lightblue","lightgrey",
              "darkmagenta","pink2","white",
              "darkgrey","beige","darkgreen","turquoise","gray","white","purple","bisque")
pchvec <- rep(c(22,21,24,23,25),10)      # symbol type
dim_num <- 2   # No of dimensions (axes) for calculating NMDS
axes <- c(1,2) # axes to be shown in the plot, could be also 3,2
envexists <- T # if FALSE no envfit takes place
decpoint <- ","  # German excel files with decimal comma, switch to "." in data.entry if needed
site <- 3
constraint_var <- "altitude"
# end variables



# chunk2: data entry: modification of default values
# first: file-names and work-directory
data.entry(spefileind,workdir,site,winsys)
# second: variables
data.entry(hel,pa_flag,minClust,pval_lev,permnum,txt_file,cov_plot)
data.entry(envexists,envname,decpoint,constraint_var,x_orient)
# end data entry



# chunk3: read abundance data
setwd(workdir)
TI <- Sys.time()
spe <- fread(spefileind,dec=decpoint,data.table = F)
rownames(spe) <- spe[,1]
spe <- spe[,-1]
if (envexists) {speinds <- grep(site,colnames(spe))
                spe <- spe[,speinds]   }
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
if (rlev_man==0) rlev <- min(rowSums(spe)) else rlev <- rlev_man
if (pa_flag==F & drare) spedr <- drarefy(spe,rlev) else spedr <- spe
# end chunk 6 further processing with spedr

# chunk 19
if (envexists)
{ # subchunk 1: read env data
  env <- fread(envname,dec=decpoint,data.table = F)
  rownames(env) <- env[,1]
  env <- env[,-1] # end subchunk 1  
  # subchunk3: adjust env data (only sites occurring in spe)
  source("adjustenv.R")
  env <- adjustenv(spe,env)
  # end subchunk 3
  # subchunk 4 sort spe according to constraint_var
  constr_ind <- grep(constraint_var,colnames(env))
  constr <- env[,constr_ind]
  spedr <- spedr[order(constr,decreasing=F),]
  spe <- spe[order(constr,decreasing=F),]
}  
# end chunk 19 read and process env data




# chunk7: generate dissimilarity matrix based on hellinger.bray_curtis.ochiai.sörensen distances
if (hel) 
{spe.h <- decostand(spedr,"hel")
if (pa_flag) dist_text <- "Ochiai" else dist_text <- "Hel" # Hellinger => Ochiai for pa
dist_method <- "euc"
}else 
{if (pa_flag) spe.h <- spedr else if (drare) spe.h <- decostand(spedr,"total") else spe.h <- decostand(spedr,margin=2,"max")
if (pa_flag) dist_text <- "Sörensen" else dist_text <- "Bray" # Bray-Curtis => Sörensen for pa
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
# end chunk 16



# chunk17: dendrogramm
source("hcoplot.R")
if (winsys) windows(title="Dendrogramm - Ward",16,12)
hcoplot(spe.dist.ward, spe.dist, k=k.best)
# end dendrogramm




# chunk18: Silhouette plot
cutg <- cutree(spe.dist.ward, k=k.best)
sil <- silhouette(cutg, spe.dist)
rownames(sil) <- row.names(spe)
maxlen <- max(length(rownames(sil)))
if (winsys) windows(title="Silhouette plot - Ward",12,12)
if (winsys) par(pin=c(5,5))
maintit <- paste("Silhouette Plot,",dist_text,spefileind)
lettersize <-  50 / nrow(spe)/2
plot(sil, max.strlen = maxlen, main=maintit,
     cex.names=lettersize, col=colorvec[1:k.best], nmax=250)
# end chunk18 Silhouette plot




# Multivariate regression trees
# *****************************
if (envexists)
{ if (winsys) windows(title="RIOJA chclust with sequential constraint", 14, 9)
  diss <- spe.dist
  clust <- chclust(diss)
  subtitle <- "CHCLUST coniss"
  maintitle <- paste("CHCLUST ",spefileind," (",dist_text,") ordered by ",constraint_var," rlev=",rlev, sep="")
  plot(clust, hang=-1,main=maintitle, sub=NULL, x.rev=x_orient)
  cutg <- cutree(clust,k=k.best)
  par(lwd=3)
  rect.hclust(clust, k=k.best, border=colorvec, cluster=cutg)
  par(lwd=1)
  legend("topright", paste("Cluster", 1:k.best), pch=pchvec[1], 
       pt.bg=colorvec[1:k.best], pt.cex=1.2,col="black",cex=0.85, bty="n")
 # if (winsys) {windows(title="bstick",12,8) 
 #                bstick(clust,12)}
}



# chunk19: indicator val
           
# *******************************************************************************
# indval characteristic OTUs for clusters
# *******************************************************************************
  #if (pval_lev==0.001) permnum=999 else permnum=9999
  #permnum <- 9999
  spe.h <- as.data.frame(spe.h)
  #indval <- multipatt(spe.h,cutg,max.order=maxorder,control=how(nperm=permnum))
  cat("calculate ",permnum," permutations (please wait..)\n")
  indvalori <- multipatt(spe.h,cutg,max.order=k.best,control=how(nperm=permnum))
  spefilenamesummary <- paste(spefileind,site,".txt",sep="")
  cat("writing summary\n")
  if (txt_file) sink(spefilenamesummary)
  summary(indvalori,alpha=pval_lev,indvalcomp=T)
  #if (txt_file) sink()
# *******************************************************************************
# end chunk19
  
  
# chunk 20: coverage
  sc <- coverage(spe.h,indvalori,At=A_t)
  names(sc) <- paste("cluster",names(sc),sep="")
  cat("Coverage [0,1] of clusters by indicative species A-threshold: ",A_t,"\n")
  print(sc)
  if (txt_file) sink()
# end chunk 20
  
  cat("runtime: ")
  print(Sys.time() - TI)  

  
# chunk21: plot coverage of cluster groups
if (cov_plot)
{  
  maintitle <- "Coverage of Site Clusters by Indicator values"
  if (winsys) windows(title=maintitle,12,12)
  plotcoverage(spe.h, indvalori, group="1", lty=1,lwd=3,col=colorvec[1],
               main=maintitle,xlab="Sensitivity Indicator A threshold (predictive)")
  for (i in 2:k.best)
    plotcoverage(spe.h, indvalori, group=as.character(i), lty=i,lwd=3, col=colorvec[i], add=TRUE)
  legend(x=0.02,y=80, paste("Cluster", c(1:k.best)), pch=pchvec[1:k.best], 
         pt.bg=colorvec[1:k.best], pt.cex=1.2,col="black",cex=0.85)
}  
# end chunk 21
  

#*********************************************************************************************************
# end of script