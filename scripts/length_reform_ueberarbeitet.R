library(data.table)
library(dplyr)

filenamein <- "christopher_lakes_data.csv"
filenameout <- "christopher_lakes_data_OTUs_zusammengefasst.csv"
workdir <- "M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen"

data.entry(filenamein,filenameout,workdir)
T1 <- Sys.time()
setwd(workdir)

spe <- fread(filenamein, data.table=F)
#index <- order(spe$staxids,decreasing=F)
#spe <- spe[index, ]

spe_staxids <- spe$staxids
spe_taxonomy <- spe$taxonomy
spe_strs <- select(spe,seqid,seqs,qlen,length,pident,mismatch,evalue, taxlevel, staxids, taxonomy)
spe <- select(spe,-seqid,-seqs,-qlen,-length,-pident,-mismatch,-evalue,-staxids, -taxlevel, -taxonomy)

spenew  <-  aggregate(spe,by=list(spe_staxids),sum)
spenew_taxonomy <- aggregate(spe_staxids, by=list(spe_taxonomy),min)
order <- order(spenew_taxonomy$x)
spenew_taxonomy <- spenew_taxonomy[order,]
spenew_taxonomy <- spenew_taxonomy[-972,]
taxonomy <- spenew_taxonomy$Group.1

spenew_group1 <- select(spenew,Group.1)
spenew <- select(spenew,-Group.1)
readsums <- rowSums(spenew)
rs_order <- order(readsums,decreasing=T)
spenew <- cbind(readsums,spenew_group1, taxonomy, spenew)

spenew <- spenew[rs_order, ]
names(spenew)[1] <- "sums"
names(spenew)[2] <- "staxids"
write.csv2(spenew,file=filenameout,row.names=F)
