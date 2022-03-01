setwd("M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen")  # workfolder for your data
library(vegan)
library(data.table)

filename_0 <- "final_lakes_TaxPrnInd_ncbi_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_0.prn"
csvmatrix_0 <- fread(filename_0, data.table=F)
csvmatrix_0 <- csvmatrix_0[,-1]
csvmatrix_0 <- apply(csvmatrix_0, FUN= function (x) as.numeric(as.character(x)), MARGIN=2)

filename_1 <- "final_lakes_TaxPrnInd_ncbi_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_1.prn"
csvmatrix_1 <- fread(filename_1, data.table=F)
csvmatrix_1 <- csvmatrix_1[,-1]
csvmatrix_1 <- apply(csvmatrix_1, FUN= function (x) as.numeric(as.character(x)), MARGIN=2)

filename_2 <- "final_lakes_TaxPrnInd_ncbi_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2.prn"
csvmatrix_2 <- fread(filename_2, data.table=F)
csvmatrix_2 <- csvmatrix_2[,-1]
csvmatrix_2 <- apply(csvmatrix_2, FUN= function (x) as.numeric(as.character(x)), MARGIN=2)

filename_3 <- "final_lakes_TaxPrnInd_ncbi_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_3.prn"
csvmatrix_3 <- fread(filename_3, data.table=F)
csvmatrix_3 <- csvmatrix_3[,-1]
csvmatrix_3 <- apply(csvmatrix_3, FUN= function (x) as.numeric(as.character(x)), MARGIN=2)

filename_4 <- "final_lakes_TaxPrnInd_ncbi_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_4.prn"
csvmatrix_4 <- fread(filename_4, data.table=F)
csvmatrix_4 <- csvmatrix_4[,-1]
csvmatrix_4 <- apply(csvmatrix_4, FUN= function (x) as.numeric(as.character(x)), MARGIN=2)

filename_5 <- "final_lakes_TaxPrnInd_ncbi_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_5.prn"
csvmatrix_5 <- fread(filename_5, data.table=F)
csvmatrix_5 <- csvmatrix_5[,-1]
csvmatrix_5 <- apply(csvmatrix_5, FUN= function (x) as.numeric(as.character(x)), MARGIN=2)

#######   Rarefaction curves of lakes       ##############################################################


A111AU_0 <- subset(csvmatrix_0, select=A111AU_B)
A111AU_1 <- subset(csvmatrix_1, select=A111AU_B)
A111AU_2 <- subset(csvmatrix_2, select=A111AU_B)
A111AU_3 <- subset(csvmatrix_3, select=A111AU_B)
A111AU_4 <- subset(csvmatrix_4, select=A111AU_B)
A111AU_5 <- subset(csvmatrix_5, select=A111AU_B)

A132OS_0 <- subset(csvmatrix_0, select=A132OS_C)
A132OS_1 <- subset(csvmatrix_1, select=A132OS_C)
A132OS_2 <- subset(csvmatrix_2, select=A132OS_C)
A132OS_3 <- subset(csvmatrix_3, select=A132OS_C)
A132OS_4 <- subset(csvmatrix_4, select=A132OS_C)
A132OS_5 <- subset(csvmatrix_5, select=A132OS_C)

A152WI_0 <- subset(csvmatrix_0, select=A152WI_B)
A152WI_1 <- subset(csvmatrix_1, select=A152WI_B)
A152WI_2 <- subset(csvmatrix_2, select=A152WI_B)
A152WI_3 <- subset(csvmatrix_3, select=A152WI_B)
A152WI_4 <- subset(csvmatrix_4, select=A152WI_B)
A152WI_5 <- subset(csvmatrix_5, select=A152WI_B)


x11(title="rarefaction curve Lake_A111AU_0")
rarecurve(t(A111AU_0))

x11(title="rarefaction curve Lake_A111AU_1")
rarecurve(t(A111AU_1))

x11(title="rarefaction curve Lake_A111AU_2")
rarecurve(t(A111AU_2))

x11(title="rarefaction curve Lake_A111AU_3")
rarecurve(t(A111AU_3))

x11(title="rarefaction curve Lake_A111AU_4")
rarecurve(t(A111AU_4))

x11(title="rarefaction curve Lake_A111AU_5")
rarecurve(t(A111AU_5))