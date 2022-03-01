# IMPORTANT: maps have to be opened individually by "run", not by "source" (results in grey window)

library(ggplot2)
library(ggmap)
library(data.table)
library(dplyr)
library(rgdal)
library(akima)
library(RColorBrewer)
library(viridis)
library(ggthemes)
library(akima)
library(vegan)

setwd("C:/Users/Chrisi/Desktop/Bachelorarbeit/Christopher-20Seen")

data <- "kingdom_Euk_without_Metazoa_and_Embryophyta_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2.agg"
spe <- fread(data, data.table=F)
row.names(spe) <- spe[,1]
spe <- spe[,-1]
spe <- decostand(spe,"hel")


kml_data <- as.data.frame(readOGR("Kampagne2012.xlsx.kml", "Kampagne2012.xlsx"))
number <- 1:nrow(kml_data)
kml_data <- as.data.frame(cbind(number, kml_data))
kml_data <- kml_data[c(83,89,92,191,220,237,253,254,268,4,6,17,25,29,47,50,54,114,120,124,127),]


lon <- kml_data$coords.x1
lat <- kml_data$coords.x2
coordinates <- cbind(lon,lat)
spe <- cbind(spe,lon,lat) 

groupname <- "Viridi"   #     <= GROUP OF INTEREST (choose from spe, doesn't have to contain full name)
spe_group <- select(spe, lon, lat, contains(groupname))
group <- as.data.frame(spe_group[,3])

map <- get_map(location="europe", maptype="satellite", zoom = 4)

max <- max(spe_group[,3])

breakmode <- "exp"  #lin/exp for linear/exponential breaks of spe_group

if(breakmode=="exp")
{
spe_group0  <- subset(spe_group, spe_group[,3] >= 0        & spe_group[,3] <= 0.01*max)
spe_group1  <- subset(spe_group, spe_group[,3] >= 0.01*max & spe_group[,3] <= 0.04*max)
spe_group2  <- subset(spe_group, spe_group[,3] >= 0.04*max & spe_group[,3] <= 0.09*max)
spe_group3  <- subset(spe_group, spe_group[,3] >= 0.09*max & spe_group[,3] <= 0.16*max)
spe_group4  <- subset(spe_group, spe_group[,3] >= 0.16*max & spe_group[,3] <= 0.25*max)
spe_group5  <- subset(spe_group, spe_group[,3] >= 0.25*max & spe_group[,3] <= 0.36*max)
spe_group6  <- subset(spe_group, spe_group[,3] >= 0.36*max & spe_group[,3] <= 0.49*max)
spe_group7  <- subset(spe_group, spe_group[,3] >= 0.49*max & spe_group[,3] <= 0.64*max)
spe_group8  <- subset(spe_group, spe_group[,3] >= 0.64*max & spe_group[,3] <= 0.81*max)
spe_group9  <- subset(spe_group, spe_group[,3] >= 0.81*max & spe_group[,3] <=      max)
}

if (breakmode=="lin")
{
spe_group0  <- subset(spe_group, spe_group[,3] >= 0        & spe_group[,3] <=    max/10)
spe_group1  <- subset(spe_group, spe_group[,3] >= max/10   & spe_group[,3] <=  2*max/10)
spe_group2  <- subset(spe_group, spe_group[,3] >= 2*max/10 & spe_group[,3] <=  3*max/10)
spe_group3  <- subset(spe_group, spe_group[,3] >= 3*max/10 & spe_group[,3] <=  4*max/10)
spe_group4  <- subset(spe_group, spe_group[,3] >= 4*max/10 & spe_group[,3] <=  5*max/10)
spe_group5  <- subset(spe_group, spe_group[,3] >= 5*max/10 & spe_group[,3] <=  6*max/10)
spe_group6  <- subset(spe_group, spe_group[,3] >= 6*max/10 & spe_group[,3] <=  7*max/10)
spe_group7  <- subset(spe_group, spe_group[,3] >= 7*max/10 & spe_group[,3] <=  8*max/10)
spe_group8  <- subset(spe_group, spe_group[,3] >= 8*max/10 & spe_group[,3] <=  9*max/10)
spe_group9  <- subset(spe_group, spe_group[,3] >= 9*max/10 & spe_group[,3] <= 10*max/10)
}

final_map <- ggmap(map, extent='device') #+

if(nrow(spe_group0) > 0)  {layer0   <- stat_density2d(data=spe_group0,  aes(x=lon, y=lat, fill = "0", alpha=..level..), bins = 6, geom = 'polygon')} else {layer0 <- NULL}
if(nrow(spe_group1) > 0)  {layer1   <- stat_density2d(data=spe_group1,  aes(x=lon, y=lat, fill = "1", alpha=..level..), bins = 6, geom = 'polygon')} else {layer1 <- NULL}
if(nrow(spe_group2) > 0)  {layer2   <- stat_density2d(data=spe_group2,  aes(x=lon, y=lat, fill = "2", alpha=..level..), bins = 6, geom = 'polygon')} else {layer2 <- NULL}
if(nrow(spe_group3) > 0)  {layer3   <- stat_density2d(data=spe_group3,  aes(x=lon, y=lat, fill = "3", alpha=..level..), bins = 6, geom = 'polygon')} else {layer3 <- NULL}
if(nrow(spe_group4) > 0)  {layer4   <- stat_density2d(data=spe_group4,  aes(x=lon, y=lat, fill = "4", alpha=..level..), bins = 6, geom = 'polygon')} else {layer4 <- NULL}
if(nrow(spe_group5) > 0)  {layer5   <- stat_density2d(data=spe_group5,  aes(x=lon, y=lat, fill = "5", alpha=..level..), bins = 6, geom = 'polygon')} else {layer5 <- NULL}
if(nrow(spe_group6) > 0)  {layer6   <- stat_density2d(data=spe_group6,  aes(x=lon, y=lat, fill = "6", alpha=..level..), bins = 6, geom = 'polygon')} else {layer6 <- NULL}
if(nrow(spe_group7) > 0)  {layer7   <- stat_density2d(data=spe_group7,  aes(x=lon, y=lat, fill = "7", alpha=..level..), bins = 6, geom = 'polygon')} else {layer7 <- NULL}
if(nrow(spe_group8) > 0)  {layer8   <- stat_density2d(data=spe_group8,  aes(x=lon, y=lat, fill = "8", alpha=..level..), bins = 6, geom = 'polygon')} else {layer8 <- NULL}
if(nrow(spe_group9) > 0)  {layer9   <- stat_density2d(data=spe_group9,  aes(x=lon, y=lat, fill = "9", alpha=..level..), bins = 6, geom = 'polygon')} else {layer9 <- NULL}

x11(50,30)
final_map +
  list(layer0, layer1, layer2, layer3, layer4, layer5, layer6, layer7, layer8, layer9) +
  scale_alpha(range = c(0.07, 0.7), guide = FALSE) +
  geom_point(aes(x=lon, y=lat), colour="black", alpha=1, size=2, data=spe_group) +
  scale_fill_manual(name= "Abundanzlevel", values=c("lightyellow", "yellow", "orange", "orangered", "red"), label=c(paste("0 -",round((0.01*max), 2)), paste(round((0.04*max), 2), "-", round((0.09*max), 2)), paste(round((0.16*max), 2), "-", round((0.25*max), 2)), paste(round((0.49*max), 2), "-", round((0.64*max), 2)), paste(round((0.81*max), 2), "-", round((max), 2)))) +   #values=brewer.pal(9, "YlOrRd")
  theme(legend.position=c(0.9,0.5)) +
  theme(legend.background = element_rect(fill="transparent")) +
  theme(legend.text = element_text(colour="white", size=12, face="bold")) +
  theme(legend.title = element_text(colour="white", size=12, face="bold")) +
  theme(legend.key = element_rect(fill="transparent", colour="transparent")) +
  annotate("text", x=15, y=60, label=groupname, color="white", size=10, fontface="bold")


##########################################################################################################

# interpolation
fld <- with(spe_group, interp(x = lon, y = lat, z = spe_group[,3]))

df <- melt(fld$z, na.rm = TRUE)
names(df) <- c("x", "y", "Abundanz")
df$lon <- fld$x[df$x]
df$lat <- fld$y[df$y]

x11(50,30)
ggmap(map, extent='device') +
  geom_tile(data = df, aes(x = lon, y = lat, fill = df[,3]), alpha=0.6) +
  scale_fill_gradient(name = "Abundance level", low = "blue", high = "red", na.value = "transparent") +
  geom_point(aes(x=lon, y=lat), colour="black", alpha=1, size=2, data=spe_group) +
  theme(legend.position=c(0.9,0.5)) +
  theme(legend.background = element_rect(fill="transparent")) +
  theme(legend.text = element_text(colour="white", size=12, face="bold")) +
  theme(legend.title = element_text(colour="white", size=12, face="bold")) +
  annotate("text", x=15, y=60, label="Viridiplantae without Embryophyta", color="white", size=10, fontface="bold")

  