setwd("M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen/")
library(dplyr)
library(vegan)
library(data.table)

#Daten einlesen
spe <- t(fread("kingdom_Euk_without_Metazoa_and_Embryophyta_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2.csv", data.table=F))
rlev_man <- 0
write_div=T
if (rlev_man==0) rlev <- min(rowSums(spe)) else rlev <- rlev_man
source("rrare100.R")
sperr <- rrarefy100(round(spe, digits=0),rlev)
##### Diversity indices  #######

Species_richness <- rowSums(sperr>0)
Shannon_entropy <- diversity(sperr)
Shannon_diversity_number <- exp(Shannon_entropy)
Simpson_diversity_number <- diversity(sperr, "inv")
Pielou_evenness <- Shannon_entropy/log(Species_richness)                     
Shannon_evenness <- Shannon_diversity_number/Species_richness
Simpson_evenness <- Simpson_diversity_number/Species_richness               
div <- data.frame(Species_richness, Shannon_entropy, Shannon_diversity_number, Simpson_diversity_number, Pielou_evenness, Shannon_evenness, Simpson_evenness) #Tabelle für alle Indices
Lakes_shortnames <- row.names(div)
div <- cbind(Lakes_shortnames, div)
if (write_div==T) {write.csv(div, file="diversity_indices.csv", row.names=F)}
View(select(div, -Lakes_shortnames))

