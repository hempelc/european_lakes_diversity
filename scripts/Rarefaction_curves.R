setwd("M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen")  # workfolder for your data
library(vegan)
library(data.table)
library(dplyr)
filename <- "best151102-103935-final_swoncbi85_tax.csv"
csvmatrix <- fread(filename, data.table=F)
csvmatrix <- select(csvmatrix, -(seqid:seqs))
#######   Rarefaction curves of lakes       ##############################################################


A111AU <- subset(csvmatrix, select=A111AU_A:A111AU_C)
A132OS <- subset(csvmatrix, select=A132OS_A:A132OS_C)
A152WI <- subset(csvmatrix, select=A152WI_A:A152WI_C)
N121SA <- subset(csvmatrix, select=N121SA_A:N121SA_C)
N261LU <- subset(csvmatrix, select=N261LU_A:N261LU_C)
O022TU <- subset(csvmatrix, select=O022TU_A:O022TU_C)
O111BA <- subset(csvmatrix, select=O111BA_A:O111BA_C)
O121RA <- subset(csvmatrix, select=O121RA_A:O121RA_C)
O241PL <- subset(csvmatrix, select=O241PL_A:O241PL_C)
S022PA <- subset(csvmatrix, select=S022PA_A:S022PA_C)
S031BU <- subset(csvmatrix, select=S031BU_A:S031BU_C)
S102LR <- subset(csvmatrix, select=S102LR_A:S102LR_C)
S153TR <- subset(csvmatrix, select=S153TR_A:S153TR_C)
S171MA <- subset(csvmatrix, select=S171MA_A:S171MA_C)
S271TO <- subset(csvmatrix, select=S271TO_A:S271TO_C)
S281PM <- subset(csvmatrix, select=S281PM_A:S281PM_C)
S301MM <- subset(csvmatrix, select=S301MM_A:S301MM_C)
Z042CP <- subset(csvmatrix, select=Z042CP_A:Z042CP_C)
Z071SI <- subset(csvmatrix, select=Z071SI_A:Z071SI_C)
Z111VV <- subset(csvmatrix, select=Z111VV_A:Z111VV_C)
Z122OU <- subset(csvmatrix, select=Z122OU_A:Z122OU_C)

Lake_1_2_3 <- cbind(A111AU, A132OS, A152WI)

Lake_4_5_6 <- cbind(N121SA, N261LU, O022TU)

Lake_7_8_9 <- cbind(O111BA, O121RA, O241PL)

Lake_10_11_12 <- cbind (S022PA, S031BU, S102LR)

Lake_13_14_15 <- cbind(S153TR, S171MA, S271TO)

Lake_16_17_18 <- cbind(S281PM, S301MM, Z042CP)

Lake_19_20_21 <- cbind (Z071SI, Z111VV, Z122OU)

x11(title="rarefaction curve Lake_1_2_3")
rarecurve(t(Lake_1_2_3))

x11(title="rarefaction curve Lake_4_5_6")
rarecurve(t(Lake_4_5_6))

x11(title="rarefaction curve Lake_7_8_9")
rarecurve(t(Lake_7_8_9))

x11(title="rarefaction curve Lake_10_11_12")
rarecurve(t(Lake_10_11_12))

x11(title="rarefaction curve Lake_13_14_15")
rarecurve(t(Lake_13_14_15))

x11(title="rarefaction curve Lake_16_17_18")
rarecurve(t(Lake_16_17_18))

x11(title="rarefaction curve Lake_19_20_21")
rarecurve(t(Lake_19_20_21))
