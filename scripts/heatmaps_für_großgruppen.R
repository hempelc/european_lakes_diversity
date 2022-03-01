library(vegan)
library(gclus)
library(gplots)
library(data.table)
library(dplyr)
library(RColorBrewer) #Für Farben (display.brewer.all())

setwd("M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen")
data <- "kingdom_Euk_without_Metazoa_and_Embryophyta_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2.agg"
spe <- fread(data, data.table=F)
row.names(spe) <- spe[,1]
spe <- spe[,-1]

ncolor <- 8 

#######
#for dendrogram of lakes

spe.h <- decostand(spe,"hel")

spe.dist <- vegdist(spe.h, "bray")

spe.dist.ward <- hclust(spe.dist, "ward.D2")

dend.lakes <-   as.dendrogram(reorder.hclust(spe.dist.ward, spe.dist))

###################################
#for dendrogram of groups

spe.h1 <- decostand(t(spe),"hel")

spe.dist1 <- vegdist(spe.h1, "bray")

spe.dist.ward1 <- hclust(spe.dist1, "ward.D2")

dend.groups <-   as.dendrogram(reorder.hclust(spe.dist.ward1, spe.dist1))
########
max <- max(spe.h)

breaks1 <- c(seq(0.01*max,0.04*max,length=1),seq(0.04*max,0.09*max,length=1),seq(0.09*max,0.16*max,length=1),seq(0.16*max,0.25*max,length=1),seq(0.25*max,0.36*max,length=1),seq(0.36*max,0.49*max,length=1),seq(0.49*max,0.64*max,length=1),seq(0.64*max,0.81*max,length=1),seq(0.81*max,1*max,length=1), seq(max,max))

x11(title="Heatmap dendrogram of lakes - colour breaks linear", 22, 20)
heatmap.2(t(spe.h), Rowv=NA, Colv=dend.lakes, dendrogram="col", density.info="none", trace="none", col=c("white", brewer.pal(ncolor, "YlOrRd")), scale="none", margin=c(8,10), keysize=0.5, lmat = rbind(c(0,3),c(2,1),c(0,4)), lwid=c(0.3,4),lhei=c(1.5,4,0.9), colCol=ifelse(rownames(spe.h) %in% c("A111AU_B", "A152WI_B", "S102LR_A", "S031BU_C", "O111BA_C", "S171MA_A", "S153TR_C", "Z042CP_C", "Z071SI_A", "Z122OU_C"), 2, 1), main="heatmap of rrarefied aggregated data - linear colourbreaks \n dendrogram of lakes")

x11(title="Heatmap dendrogram of lakes - colour breaks exponential", 22, 20)
heatmap.2(t(spe.h), Rowv=NA, Colv=dend.lakes, dendrogram="col", density.info="none", trace="none", breaks=breaks1, col=c("white", brewer.pal(ncolor, "YlOrRd")), scale="none", margin=c(8,10), keysize=0.5, lmat = rbind(c(0,3),c(2,1),c(0,4)), lwid=c(0.3,4),lhei=c(1.5,4,0.9), colCol=ifelse(rownames(spe.h) %in% c("A111AU_B", "A152WI_B", "S102LR_A", "S031BU_C", "O111BA_C", "S171MA_A", "S153TR_C", "Z042CP_C", "Z071SI_A", "Z122OU_C"), 2, 1), main="heatmap of rrarefied aggregated data - exponential colourbreaks \n dendrogram of lakes")

x11(title="Heatmap dendrogram of groups - colour breaks exponential", 22, 20)
heatmap.2(t(spe.h), Rowv=dend.groups, Colv=NA, dendrogram="row", density.info="none", trace="none", col=c("white", brewer.pal(ncolor, "YlOrRd")), scale="none", margin=c(8,10), keysize=0.5, lmat = rbind(c(0,3),c(2,1),c(0,4)), lwid=c(0.8,4),lhei=c(0.5,4,0.7), colCol=ifelse(rownames(spe.h) %in% c("A111AU_B", "A152WI_B", "S102LR_A", "S031BU_C", "O111BA_C", "S171MA_A", "S153TR_C", "Z042CP_C", "Z071SI_A", "Z122OU_C"), 2, 1), main="heatmap of rrarefied aggregated data - exponential colourbreaks \n dendrogram of groups")


jens.dendrogram <- spe.h[,c(2,3,1,4,7,8,9,6,5,10,11,12,24,22,23,15,16,17,14,13,18,19,20,21)]

x11(title="Heatmap dendrogram of groups - colour breaks exponential", 22, 20)
heatmap.2(t(data.matrix(jens.dendrogram)), Rowv=NA, Colv=NA, dendrogram="none", density.info="none", trace="none", col=c("white", brewer.pal(ncolor, "YlOrRd")), scale="none", margin=c(8,10), keysize=0.5, lmat = rbind(c(0,3),c(2,1),c(0,4)), lwid=c(0.8,4),lhei=c(0.5,4,0.7), colCol=ifelse(rownames(spe.h) %in% c("A111AU_B", "A152WI_B", "S102LR_A", "S031BU_C", "O111BA_C", "S171MA_A", "S153TR_C", "Z042CP_C", "Z071SI_A", "Z122OU_C"), 2, 1), main="heatmap of rrarefied aggregated data - exponential colourbreaks \n dendrogram of groups")

spe.h <- as.data.frame(t(spe.h))

high_lakes <- t(select(spe.h, A111AU_B, A152WI_B, S102LR_A, S031BU_C, O111BA_C, S171MA_A, S153TR_C, Z042CP_C, Z071SI_A, Z122OU_C))
low_lakes <- t(select(spe.h, -A111AU_B, -A152WI_B, -S102LR_A, -S031BU_C, -O111BA_C, -S171MA_A, -S153TR_C, -Z042CP_C, -Z071SI_A, -Z122OU_C))

high_lakes <- colSums(high_lakes)
low_lakes <- colSums(low_lakes)

lakes_high_low_sums <- cbind(high_lakes, low_lakes)

max2 <- max(lakes_high_low_sums)
breaks2 <- c(seq(0.01*max2,0.04*max2,length=1),seq(0.04*max2,0.09*max2,length=1),seq(0.09*max2,0.16*max2,length=1),seq(0.16*max2,0.25*max2,length=1),seq(0.25*max2,0.36*max2,length=1),seq(0.36*max2,0.49*max2,length=1),seq(0.49*max2,0.64*max2,length=1),seq(0.64*max2,0.81*max2,length=1),seq(0.81*max2,1*max2,length=1), seq(max2,max2))

x11(title="Heatmap summarized lakes", 10, 20)
heatmap.2(lakes_high_low_sums, Rowv=NA, Colv=NA, dendrogram="none", density.info="none", trace="none", breaks=breaks2, col=c("white", brewer.pal(ncolor, "YlOrRd")), scale="none", margin=c(8,10), keysize=0.5, lmat = rbind(c(0,3),c(2,1),c(0,4)), lwid=c(0.1,1),lhei=c(0.3,4,0.7), cexCol=2.5, srtCol=360, adjCol=c(NA,1.2))
