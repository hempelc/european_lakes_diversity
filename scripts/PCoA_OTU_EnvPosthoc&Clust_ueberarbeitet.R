# PCoA and ENVIRONMENT on rarefied data distance Hellinger or Bray ...
# V2.6 Author M.Jensen, June 8, 2016
# Biodiversity Group University of Duisburg-Essen, Germany
# based on Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
#**************************************************************************************


# chunk0: default filenames etc
#library(vegan3d) 3d plot possible with ordiplot3d, see line 244, but not very beautiful
library(gclus)
library(vegan)
library(dplyr)
library(data.table)
library(ellipse)
workdir <- "C:/Users/Chrisi/Desktop/Bachelorarbeit/Christopher-20Seen"
spe1fileprn <- "kingdom_Euk_without_Metazoa_and_Embryophyta_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2.csv"
taxfile <- "kingdom_Euk_without_Metazoa_and_Embryophyta_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2.tax"
envname <- "lakes_abiotic_data_just_numbers_to_work_with.csv" # name of environmental file
# end chunk 0 (files must be in the working directory)


# chunk1: some variables
# flags
rrare <- T  # no classical rarefaction
hel   <- F  # default: Hellinger/euclidic distance, if FALSE: Bray-Curtis distance
pa_flag <- F # presence-absence flag: all data are transformed to 0,1 and no rarefaction
minClust <- F # if False, the large cluster no is chosen according to both optimum criteria
winsys <- T # in case of Linux please set to FALSE
dendro <- T # if FALSE, no integrated dendrogram within NMDS plot
print_divs <- F
rlev_man <- 0
add_spe1c <- T
# variables for plotting
legend_later <- F
clust_ellip <- F
ellip_rad <- 0.95
legend_pos <- "bottomleft"
x_orient <- T
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
# end variables



# chunk2: data entry: modification of default values
# first: file-names and work-directory
data.entry(spe1fileprn,workdir,winsys)
data.entry(envexists,envname,decpoint)
# second: variables
data.entry(hel,pa_flag,add_spe1c,print_divs,rlev_man,dim_num,axes,minClust)
data.entry(legend_pos,legend_later,dendro,clust_ellip,ellip_rad,x_orient,y_orient,winsys)
# end data entry



# chunk3: read abundance data
setwd(workdir)
TI <- Sys.time()
spe1 <- fread(spe1fileprn,dec=decpoint,data.table = F)
tax <- read.csv(taxfile)
# end read data

rn <- paste("N",1:nrow(spe1),sep="") # =row.names(spe1)
rownames(spe1) <- rn
rownames(tax) <- rn

# chunk4: data structure: how many OTUs only in 1,2..10 lakes
spe1pa <- decostand(spe1,"pa")
spe1pasum <- rowSums(spe1pa)      # number of sites in which the OTUs (rows) occur...
sitecounts <- function(i,spe1pasum){sum(as.integer(spe1pasum==i))}
sitecounts10 <- sapply(1:10,sitecounts,spe1pasum=spe1pasum)
cat("community structure:\nhow many OTUs restricted to 1..10 sites: ",sitecounts10,"\n")
# end structure info


# chunk5: transpose data (for vegan: samples in rows and OTUs in columns) and pa transform if needed
spe1 <- t(spe1)
cat("number of sites: ",nrow(spe1),"\n")
if (pa_flag) spe1 <- decostand(spe1,"pa")
# now vegan functions can be used



# chunk6: rarefaction of sites with max possible level
source("rrare100.R")
if (rlev_man==0) rlev <- min(rowSums(spe1)) else rlev <- rlev_man
if (pa_flag==F & rrare) spe1rr <- spe1rrare <- rrarefy100(round(spe1, digits=0),rlev) else spe1rr <- spe1
# end chunk 6 further processing with spe1rr




# chunk7: generate dissimilarity matrix based on hellinger.bray_curtis.ochiai.sÃ¶rensen distances
if (hel) 
{spe1.h <- decostand(spe1rr,"hel")
if (pa_flag) dist_text <- "Ochiai" else dist_text <- "Hel" # Hellinger => Ochiai for pa
dist_method <- "euc"
}else 
{if (pa_flag) spe1.h <- spe1rr else if (rrare) spe1.h <- decostand(spe1rr,"total") else spe1.h <- decostand(spe1rr,margin=2,"max")
 if (pa_flag) dist_text <- "SÃ¶rensen" else dist_text <- "Bray" # Bray-Curtis => SÃ¶rensen for pa
 dist_method <- "bray"
}
spe1.dist <- vegdist(spe1.h,dist_method)
# end chunk7




# chunk 8 print diversity data
if (print_divs)
{ source("rarefy_putback (für Skript PCAOTU_EnvPosthoc&Clust).R")    # function rrarefy but with replace=TRUE. large rlevel possible
  spe1rr <- rarefy_putback(spe1,rlev)            # needed for diversity indexspe1rr <- spe1
  N0 <- rowSums(spe1rr > 0)         # spe1cies richness
  H <- diversity(spe1rr)            # Shannon entropy
  N1 <- exp(H)                     # Shannon diversity (number of abundant spe1cies)
  N2 <- diversity(spe1rr, "inv")    # Simpson diversity (number of dominant spe1cies)
  Pielou <- H/log(N0)                   # Pielou evenness
  Shan_eve <- N1/N0                     # Shannon evenness (Hill's ratio)
  Simp_eve <- N2/N0                     # Simpson evenness (Hill's ratio)
  div <- data.frame(N0, H, N1, N2, Shan_eve, Simp_eve, Pielou)
  colnames(div)[1:4] <- c("spe1c_R","Shan_Ent","Shan_div","Simp_div")
  View(div)
}
# end chunk 8


# chunk11: Ward clustering of dissimilarity matrix and extraction of groups
spe1.dist.ward <- hclust(spe1.dist, "ward.D2")
# end


# chunk12: optimal number of clusters according to silhouette widths
# (Rousseeuw quality index)
# First, create an empty vector in which the asw values will be written
asw <- numeric(nrow(spe1))
for (k in 2:(nrow(spe1)-1)) {
  sil <- silhouette(cutree(spe1.dist.ward, k=k), spe1.dist)
  asw[k] <- summary(sil)$avg.width
}
k.best <- which.max(asw)
# end optimal number




# chunk13: Plot average silhouette widths (using Ward clustering) for all partitions 
# except for the trivial partition in a single group (k=1)
if (winsys) windows(title="Numbers: Silhouettes - Ward - k = 2 to n-1",12,12)
plot(1:nrow(spe1), asw, type="h", 
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
kt <- data.frame(k=1:nrow(spe1), r=0)
for (i in 2:(nrow(spe1)-1)) {
  gr <- cutree(spe1.dist.ward, i)
  distgr <- grpdist(gr)
  mt <- cor(spe1.dist, distgr, method="pearson")
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
source("hcoplot (für Skript NMDSClust und PCAEnv&Clust).R")
if (winsys) windows(title="Dendrogramm - Ward",16,12)
par(mar=c(7,4,4,6)+0.1)
hcoplot(spe1.dist.ward, spe1.dist, k=k.best)
# end dendrogramm




# chunk18: Silhouette plot
cutg <- cutree(spe1.dist.ward, k=k.best)
sil <- silhouette(cutg, spe1.dist)
rownames(sil) <- row.names(spe1)
maxlen <- max(length(rownames(sil)))
if (winsys) windows(title="Silhouette plot - Ward",12,12)
if (winsys) par(pin=c(5,5))
maintit <- paste("Silhouette Plot,",dist_text,spe1fileprn)
lettersize <-  50 / nrow(spe1)/2
plot(sil, max.strlen = maxlen, main=maintit,
     cex.names=lettersize, col=colorvec[1:k.best], nmax=250)
# end chunk18 Silhouette plot











#*************************************************************************
# chunk 19 calculate PCoA
if (dist_method=="Bray") corr_meth <- "cailliez" else corr_meth <- "none"
dim_max <- nrow(spe1)-1
if (dim_num > dim_max) dim_num <- dim_max
spe1.b.pcoa <- wcmdscale(spe1.dist, k=dim_max, eig=T)
# end chunk 19
#*************************************************************************



# chunk19a: Shepard plot and goodness of fit
if (winsys) windows(title="PCoA - Shepard plot", 12, 12)
stress_title <- paste("Shepard plot of PCoA sites (dim=",dim_num,", ",dist_text," Dissim.)",sep="")
stressplot(spe1.b.pcoa, k=dim_num, main=stress_title)
gof <- spe1.b.pcoa$GOF
# end chunk19a


  
# chunk 20 plot eigenvalues and % of variance for each axis
source("evplot (für Skript PCAEnv&Clust und PCAOTU_EnvPosthoc&Clust).R")
ev <- spe1.b.pcoa$eig
if (winsys) windows(title="PCoA eigenvalues")
evplot(ev)
cat("\nPercentage/PCoA-Axis","\n")
print(100*ev/sum(ev))
# end chunk 20

  
  
  
# chunk 21 Plot of the sites and weighted average projection of spe1cies
if (winsys) windows(title="PCoA on OTUs",12,12)
sit.sc1 <- scores(spe1.b.pcoa)#[,c(1:dim_num)]
xrange <- range(sit.sc1[,axes[1]])
yrange <- range(sit.sc1[,axes[2]])
if (x_orient==F) xrange <- rev(xrange)
if (y_orient==F) yrange <- rev(yrange)
titletext <- 
      paste("PCoA + Ward-clusters ",spe1fileprn," (",dist_text," dist.) ", "rlev: ",rlev,sep="")
p <- ordiplot(sit.sc1, scaling=1,choices=axes,type="none",main=titletext,xlim=xrange,ylim=yrange)
# Add spe1cies
if (add_spe1c)
{ spe1.wa <- wascores(spe1.b.pcoa$points[,axes], spe1rr)
text(spe1.wa, rownames(spe1.wa), cex=0.7, col=colorvec[9])
}
abline(h=0, lty=3)
abline(v=0, lty=3)
spe1.bw.groups <- cutree(spe1.dist.ward, k=k.best)
grp.lev <- levels(factor(spe1.bw.groups))
for (i in 1:length(grp.lev))
{ points(sit.sc1[spe1.bw.groups==i,axes[1]],sit.sc1[spe1.bw.groups==i,axes[2]], 
         pch=pchvec[i], cex=1.2, col="black",bg=colorvec[i]) }
text(sit.sc1[ ,axes[1]], sit.sc1[ ,axes[2]],row.names(spe1), pos=4, cex=0.7, col="black")
# Add the dendrogram if dendro=TRUE
if (dendro & add_spe1c==F) ordicluster(p, spe1.dist.ward, col="dark grey")
if (legend_later==F)
legend(legend_pos, paste("Cluster", c(1:length(grp.lev))), pch=pchvec[1:length(grp.lev)], 
       pt.bg=colorvec[1:length(grp.lev)], pt.cex=1.2,col="black",cex=0.85)
# end chunk21 plot results


# chunk21a: draw ellipses around clusters
if (clust_ellip)
{sit_2d <- cbind(sit.sc1[,axes[1]],sit.sc1[,axes[2]])
for(i in 1:length(grp.lev))
{	sit_2di <- sit_2d[spe1.bw.groups==i,]
  cov <- cov(sit_2di)
  centre <- apply(sit_2di, 2, mean)
  lines(ellipse(cov, centre=centre, level=ellip_rad),col=colorvec[i]) }}
# end chunk21a




# chunk22: a posteriori projection of environmental variables in PCoA
# log nolog etc have to be modified depending on the env-file
if (envexists)
{ env <- fread(envname,dec=decpoint, data.table=F)
row.names(env) <- env[,1]
env <- select(env, -Seekürzel)
source("adjustenv (für Skript NMDSClust und PCAOTU_EnvPosthoc&Clust).R")
env <- adjustenv(spe1,env)
#if (pa_flag!=T) env$spe1c_rich <- rowSums(spe1rr)
env <- scale(env) # z-transformations are necessary here ...
spe1.nmds.env <- envfit(spe1.b.pcoa, env, choices=axes, na.rm=T, permutations=9999)
print(spe1.nmds.env)
plot(spe1.nmds.env,col="cyan")
# Plot significant variables with a different colour
plot(spe1.nmds.env, p.max=0.001, col="blue")
} #end chunk25 env posthoc


cat("runtime: ")
print(Sys.time() - TI)

if (legend_later)
  legend(locator(1), paste("Cluster", c(1:length(grp.lev))), pch=pchvec[1:length(grp.lev)], 
         pt.bg=colorvec[1:length(grp.lev)], pt.cex=1.2,col="black",cex=0.85)

detach("package:vegan")
detach("package:gclus")
detach("package:dplyr")
detach("package:data.table")
detach("package:ellipse")
#******************** end of script ******************************************************