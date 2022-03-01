# CPCoA & lin.discriminant.analysis on rarefied data distance Hellinger or Bray ...
# V2.9 Author M.Jensen, June 27, 2016
# Biodiversity Group University of Duisburg-Essen, Germany
# based on Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
#**************************************************************************************


# chunk0: default filenames etc
library(vegan)
#library(vegan3d) 3d plot possible with ordiplot3d, see line 244, but not very beautiful
library(MASS)
library(gclus)
library(dplyr)
library(data.table)
library(ellipse)
workdir <- "M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen"
spefileprn <- "christopher_lakes_data_OTUs_zusammengefasst.csv"
envname <- "lakes_abiotic_data_just_numbers_to_work_with.csv" # name of environmental file
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
add_TP.DN_DN.DOC <- F
txtprotocol <- F
lc_plot <- T
# variables for plotting
perm_steps <- 1000
legend_later <- F
clust_ellip <- T
ellip_rad <- 0.95
legend_pos <- "bottomleft"
x_orient <- T
y_orient <- T
x_orient_LDA <- T
y_orient_LDA <- T
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
data.entry(spefileprn,workdir,txtprotocol,winsys)
fileprotocol <- paste(spefileprn,".CPCoA.txt",sep="")
if (txtprotocol) sink(fileprotocol)
data.entry(envname,decpoint,add_TP.DN_DN.DOC,perm_steps)
# second: variables
if (envexists) dendro <- F
data.entry(hel,pa_flag,print_divs,rlev_man,dim_num,axes,minClust)
data.entry(legend_pos,legend_later,ellip_rad,x_orient,y_orient)
# end data entry



# chunk3: read abundance data
setwd(workdir)
TI <- Sys.time()
spe <- fread(spefileprn,dec=decpoint, data.table=F)
spe <- spe[,-1:-3]
#spe <- as.data.frame(apply(spe, 2, FUN=function (x) as.numeric(as.character(x))))
# end read data


# chunk4: data structure: how many OTUs only in 1,2..10 lakes
spepa <- decostand(spe,"pa")
spepasum <- rowSums(spepa)      # number of sites in which the OTUs (rows) occur...
sitecounts <- function(i,spepasum){sum(as.integer(spepasum==i))}
sitecounts10 <- sapply(1:10,sitecounts,spepasum=spepasum)
cat("community structure:\nhow many OTUs restricted to 1..10 sites: ",sitecounts10,"\n")
cat("No of OTUs: ",nrow(spe),"\n")
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
cat("\nChosen distance_method: ",dist_text,"\n\n")
# end chunk7




# chunk 8 print diversity data
if (print_divs)
{ source("rarefy_putback.R")    # function rrarefy but with replace=TRUE. large rlevel possible
  sperr <- rarefy_putback(spe,rlev)            # needed for diversity indexspedr <- spe
  N0 <- rowSums(sperr > 0)         # Species richness
  H <- diversity(sperr)            # Shannon entropy
  N1 <- exp(H)                     # Shannon diversity (number of abundant species)
  N2 <- diversity(sperr, "inv")    # Simpson diversity (number of dominant species)
  Pielou <- H/log(N0)                   # Pielou evenness
  Shan_eve <- N1/N0                     # Shannon evenness (Hill's ratio)
  Simp_eve <- N2/N0                     # Simpson evenness (Hill's ratio)
  div <- data.frame(N0, H, N1, N2, Shan_eve, Simp_eve, Pielou)
  colnames(div)[1:4] <- c("Spec_R","Shan_Ent","Shan_div","Simp_div")
  print(div)
}
# end chunk 8


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
maintit <- paste("Silhouette Plot,",dist_text,spefileprn)
lettersize <-  50 / nrow(spe)/2
plot(sil, max.strlen = maxlen, main=maintit,
     cex.names=lettersize, col=colorvec[1:k.best], nmax=250)
# end chunk18 Silhouette plot



# chunk 19a
if (envexists)
{ # subchunk 1: read env data
  env <- fread(envname,dec=decpoint,data.table = F)
  rownames(env) <- env[,1]
  env <- env[,-1]
  # end subchunk 1  
  
  # subchunk 2: log nolog etc have to be modified depending on the env-file
  # View(head(env))
  numbers <- c(1:ncol(env))
  nonlogs <- c(1, 2, 3, 20, 21, 23)
  cat("Select Variable Column Numbers of non-log variables\n")
  var_names <- colnames(env)
  data.entry(numbers, var_names, nonlogs)
  envnonlog <- env[, nonlogs]
  envlog <- env[, -nonlogs]
  envlog <- log1p(envlog)
  envscaled   <- cbind(envnonlog,envlog)
  # special additional env_vars possible
  if (add_TP.DN_DN.DOC) { envscaled$TP.DN <- envscaled$TP / envscaled$DN
                          envscaled$DN.DOC <- envscaled$DN / envscaled$DOC  }
  # end subchunk 2 log - nonlog - sort
  
  # subchunk3: adjust env data (only sites occurring in spe)
  source("adjustenv.R")
  envscaled <- adjustenv(spe,envscaled)
  envscaled$spec_rich <- rowSums(spedr)
  envscaled <- as.data.frame(scale(envscaled)) # z-transformations are necessary here ...
  # end subchunk 3
  
}  
# end chunk 19 read and process env data






#****************************************************************************
# chunk 19 calculate constrainedPCoA
if (dist_method=="Bray") corr_meth <- "cailliez" else corr_meth <- "none"
dim_max <- nrow(spe)-1
if (dim_num > dim_max) dim_num <- dim_max
spe.b.pcoa <- capscale(spe.dist ~ ., data = envscaled, dist = dist_method,add = T)
var_infl_fact <- vif.cca(spe.b.pcoa)
# subchunk 4 select Variable Column Numbers (vcn) = Variables of Interest
numbers <- c(1:ncol(envscaled))
varnums_of_interest <- c(1:26)
cat("Select Variable Column Numbers of Interest out of ",ncol(envscaled),"variables\n")
var_names <- colnames(envscaled)
data.entry(numbers, var_names, var_infl_fact,varnums_of_interest)
envsub <- envscaled[,varnums_of_interest]
spe.b.pcoa <- capscale(spe.dist ~ ., data = envsub, dist = dist_method,add = T)
# end chunk 19
#****************************************************************************



# chunk19a: Shepard plot and goodness of fit
if (winsys) windows(title="PCoA - Shepard plot", 12, 12)
stress_title <- paste("Shepard plot of PCoA sites (dim=",dim_num,", ",dist_text," Dissim.)",sep="")
stressplot(spe.b.pcoa, k=dim_num, main=stress_title)
gof <- spe.b.pcoa$GOF
# end chunk19a


  
# chunk 20 plot eigenvalues and % of variance for each axis
#source("evplot.R")
#ev <- spe.b.pcoa$CA$eig
#if (winsys) windows(title="PCoA eigenvalues",12,12)
#evplot(ev)
#cat("\nPercentage/PCoA-Axis","\n")
#print(100*ev/sum(ev))
# end chunk 20

#chunk 21
# function plotCPCoA  
plotCPCoA <- function(pcoa)
{# chunk 21: plot the sites with cluster symbols and colours (scaling 1)
 pcatitle <- paste("CPCoA (",dist_text,", lc) scal 1; ", spefileprn,
                   "; rlev: ",as.character(rlev),sep="")
 if (winsys) windows(title="only for scores",4,4)
 p <- plot(pcoa, choices=axes, display=c("lc","cn"), scaling=1, type="n",main="only for scores")
 sit.sc1 <- p$constraints
 xrange <- range(sit.sc1[,1])
 yrange <- range(sit.sc1[,2])
 if (x_orient==F) xrange <- rev(xrange)
 if (y_orient==F) yrange <- rev(yrange)
 if (winsys) 
 { windows(title="Ordination and clustering",32,20)
   par(mfrow=c(1,2))
 }
 p <- plot(pcoa, choices=axes, display=c("lc","cn"), scaling=1, type="text",
           main=pcatitle, xlim=xrange, ylim=yrange)
 text(sit.sc1, row.names(spe), cex=0.7, col="white",font=2)
 abline(v=0, lty="dotted")
 abline(h=0, lty="dotted")
 spe.bw.groups <- cutree(spe.dist.ward, k=k.best)
 grp.lev <- levels(factor(spe.bw.groups))
 for (i in 1:k.best)
 { points(sit.sc1[spe.bw.groups==i,1],sit.sc1[spe.bw.groups==i,2], 
          pch=pchvec[i], cex=1.2, col="black",bg=colorvec[i]) }
 text(sit.sc1,row.names(spe), pos=4, cex=0.7, col="black")
 # Add the dendrogram if dendro=TRUE
 if (dendro) ordicluster(p, spe.dist.ward, col="dark grey")
 if (legend_later==F) legend(legend_pos, paste("Cluster", c(1:k.best)), pch=pchvec[1:k.best], 
                             pt.bg=colorvec[1:k.best], pt.cex=1.2,col="black",cex=0.85)
 if (clust_ellip)
 {sit_2d <- cbind(sit.sc1[,1],sit.sc1[,2])
 for(i in 1:k.best)
 {	sit_2di <- sit_2d[spe.bw.groups==i,]
 cov <- cov(sit_2di)
 centre <- apply(sit_2di, 2, mean)
 lines(ellipse(cov, centre=centre, level=ellip_rad),col=colorvec[i]) }}
} #  end chunk 21 plot the sites with cluster symbolsend function plotCPCoA



# chunk 22 tests
print(anova(spe.b.pcoa))
print(anova(spe.b.pcoa,by="axis"))
cat("vif:\n")
print(vif.cca(spe.b.pcoa))
R2a.all <- RsquareAdj(spe.b.pcoa)
print(R2a.all)
# end chunk 22



# chunk 23 ordistep selected subset of env variables by step both directions
spe.b.1 <- capscale(spe.dist ~ 1, data = envsub, dist = dist_method,add = T)
spe.b.all <- capscale(spe.dist ~ ., data = envsub, dist = dist_method,add = T)
print(step.forward <- ordistep(spe.b.1,scope=formula(spe.b.all),direction="both", pstep=perm_steps))
print(anova(step.forward))
print(anova(step.forward,by="axis"))
vifs <- vif.cca(step.forward)
cat("vif:\n")
print(vifs)
# end chunk 23





# Linear discriminant analysis (LDA)
# **********************************

# chunk 24: test dispersion
gr <- cutg
#cutree(hclust(vegdist(spe.hel, "euc"), "ward.D2"), k.best)	# If not run above
env.pars3 <- envsub
env.pars3.d1 <- dist(env.pars3)
env.MHV2 <- betadisper(env.pars3.d1, gr)
#print(env.MHV2)
ptest <- permutest(env.MHV2)
print(ptest)
# end chunk 24


# chunk 25: Computation of LDA (discrimination)
env.pars3.df <- as.data.frame(env.pars3)
spe.lda <- lda(gr ~ ., data=env.pars3.df)
# The result object contains the information necessary to interpret the LDA
summary(spe.lda)
# Display the group means for the 3 variables
spe.lda$means
# Compute the normalized eigenvectors (matrix C, eq. 11.33)
# which are the standardized discriminant function coefficients
Cs <- spe.lda$scaling
#print(Cs)
# Compute the canonical eigenvalues
spe.lda$svd^2
# Position the objects in the space of the canonical variates
(Fp <- predict(spe.lda)$x)
# alternative way: Fp <- scale(env.pars3.df, center=TRUE, scale=FALSE) %*% C
# Classification of the objects
spe.class <- predict(spe.lda)$class
#print(spe.class)
# Posterior probabilities of the objects to belong to the groups
spe.post <- predict(spe.lda)$posterior
#print(spe.post)
# end chunk 25 LDA


# chunk 26: Contingency table of prior versus predicted classifications
cat("\nContingency table of prior versus predicted classifications\n")
spe.table <- table(gr, spe.class)
print(spe.table)
# Proportion of correct classification
cat("Proportion of correct classification\n")
print(diag(prop.table(spe.table, 1)))
# end chunk 26



# chunk 27: plot the sites with cluster symbols and colours (scaling 1)
plotLDA <- function(lda,envsub)
{ varnames <- colnames(envsub)[1]
  for (i in 2:length(envsub)) varnames <- paste(varnames,colnames(envsub)[i])
  lindisktitle <- paste("LDA; p_hom: ", ptest$tab$Pr[1],
                        "; Variables:\n", varnames,sep="")
  p <- plot(Fp[, 1], Fp[, 2], type="n",main=lindisktitle)
  abline(v=0, lty="dotted")
  abline(h=0, lty="dotted")
  spe.bw.groups <- gr
  grp.lev <- levels(factor(spe.bw.groups))
  for (i in 1:k.best)
  { points(Fp[spe.bw.groups==i,axes[1]],Fp[spe.bw.groups==i,axes[2]], 
           pch=pchvec[i], cex=1.2, col="black",bg=colorvec[i])  }
  text(Fp[ ,axes[1]], Fp[ ,axes[2]],row.names(spe), pos=4, cex=0.7, col="black")
  if (legend_later==F)
    legend(legend_pos, paste("Cluster", c(1:k.best)), 
         pch=pchvec[1:k.best], 
         pt.bg=colorvec[1:k.best], pt.cex=1.2,col="black",cex=0.85)
  # Draw 95% ellipses around the groups
  for(i in 1:k.best)
  {	cov <- cov(Fp[gr==i, ])
    centre <- apply(Fp[gr==i, ], 2, mean)
    lines(ellipse(cov, centre=centre, level=ellip_rad),col=colorvec[i]) }
}
# end chunk 27
              

# chunk 28: LDA with jackknife-based classification (i.e., leave-one-out cross-validation)
spe.lda.jac <- lda(gr ~ ., data=env.pars3.df, add=T, CV=TRUE)
#summary(spe.lda.jac)
# Numbers and proportions of correct classification
spe.jac.class <- spe.lda.jac$class
spe.jac.table <- table(gr, spe.jac.class)
#diag(prop.table(spe.jac.table, 1))
cat("\nLDA with jackknife-based classification (i.e., leave-one-out cross-validation)\n")
cat("Numbers and proportions of correct classification\n")
print(spe.jac.table)
# end chunk 28


# chunk 29 plot with envsub
plotCPCoA(spe.b.pcoa) # by function chunk21,21a,21b
plotLDA(spe.lda,envsub)
if (legend_later)
  legend(locator(1), paste("Cluster", c(1:k.best)), pch=pchvec[1:k.best], 
         pt.bg=colorvec[1:k.best], pt.cex=1.2,col="black",cex=0.85)
# end chunk 29


# chunk 30 plot with results of ordistep, create envsubsub
subsubvars <- names(vifs)
subsubinds <- match(subsubvars,colnames(envsub))
envsubsub <- envsub[, subsubinds]
spe.b.pcoa.subsub <- capscale(spe.dist ~ ., data = envsubsub, dist = dist_method,add = T)
plotCPCoA(spe.b.pcoa.subsub)
spe.lda.subsub <- lda(gr ~ ., data=envsubsub)
plotLDA(spe.lda.subsub,envsubsub)
if (legend_later)
  legend(locator(1), paste("Cluster", c(1:k.best)), pch=pchvec[1:k.best], 
         pt.bg=colorvec[1:k.best], pt.cex=1.2,col="black",cex=0.85)
spe.lda.jac <- lda(gr ~ ., data=envsubsub, add= T, CV=TRUE)
spe.jac.class <- spe.lda.jac$class
spe.jac.table <- table(gr, spe.jac.class)
#diag(prop.table(spe.jac.table, 1))
cat("\nreduced Var number according to ordistep\n")
cat("Variables: ",colnames(envsubsub))
cat("\nLDA with jackknife-based classification (i.e., leave-one-out cross-validation)\n")
cat("Numbers and proportions of correct classification\n")
print(spe.jac.table)
# end chunk 30


if (txtprotocol) sink()
#if (txtprotocol) sink()
cat("runtime: ")
print(Sys.time() - TI)


detach("package:vegan")
detach("package:gclus")
detach("package:dplyr")
detach("package:data.table")
detach("package:ellipse")
detach("package:MASS")
#******************** end of script ******************************************************