require(data.table)
require(dplyr)
require(vegan)


# chunk0: filenames etc
workdir <- "M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen/"
filename <- "christopher_lakes_data_OTUs_zusammengefasst"
# end filenames

# chunk1: some variables
filt_kingd <- "Euk"
del_group2 <- "Embryophyta"
del_group1 <- "Metazoa"
excl1_flag <- F
exclusive_group1 <- "Chrysophyceae"
excl2_flag <- F
exclusive_group2 <- "Synurophyceae"
excl3_flag <- F
exclusive_group3 <- "Oikomonadaceae"
write_exclusive_file <- F
write_IND <- F
tax_cuts <- T
hel <- T

otusumlevel <- 2 # =1 means singletons vanish etc. 
lakesumlevel <- 100 # =1000 means samples with reads < 1001 vanish
# end vars


# chunk2: data entry
# first: file-names and work-directory
data.entry(filename,workdir)
# second: variables, here mainly for filtering
data.entry(hel, filt_kingd,del_group1,del_group2,excl1_flag,exclusive_group1,
           excl2_flag,exclusive_group2,excl3_flag,exclusive_group3,
           write_exclusive_file, write_IND, tax_cuts)

# end data entry

TI <- Sys.time() # to find out runtime of script

# chunk3: read data
setwd(workdir)
otufilename <- paste(filename,".csv",sep="")
spe <- fread(otufilename, data.table=F)
# end read data


# chunk4: examine table and replace ", " in taxonomy by ";"
spe$taxonomy <- gsub("; ",";",spe$taxonomy)
# end inspection


# chunk5: change data format of 2 variables from character to numeric

#spe$pident <- gsub(",", ".",spe$pident)
#spe$evalue <- gsub(",", ".",spe$evalue)
#spe$pident <- as.numeric(spe$pident)
#spe$evalue <- as.numeric(spe$evalue)

# end change data format

# chunk6 filter OTUs with "subset": only Eukaryota (beginning with "Euk")
#                         ********
len <- nchar(filt_kingd)
spe <- subset(spe,substring(taxonomy,1,len)==filt_kingd)
# end filter groupname


# chunk7 filter OTUS with grep, for which indices is delfiltergroup substring of taxonomy ?

del1_indices <- grep(del_group1,spe$taxonomy)
if (length(del1_indices)>0)  spe <- spe[-del1_indices,]
del2_indices <- grep(del_group2,spe$taxonomy)
if (length(del2_indices)>0)  spe <- spe[-del2_indices,]
# end filter with grep


# chunk8: select and deselect columns (variables)

spe_seqid_seqs <- select(spe,sums:taxonomy)
spe <- select(spe,-c(sums:taxonomy))
rownames(spe) <- 1:nrow(spe)


# chunk9: select specific sites (if necessary)

spe <- select(spe, A111AU_B, A132OS_C, A152WI_B, N121SA_C, N261LU_C,
              O022TU_A, O111BA_C, O121RA_C, O241PL_C, S022PA_A, 
              S031BU_C, S102LR_A, S153TR_C, S171MA_A,S271TO_C, 
              S281PM_C, S301MM_C, Z042CP_C, Z071SI_A, Z111VV_C, Z122OU_C)



# chunk10: dataentry filter OTUsum etc

data.entry(otusumlevel,lakesumlevel)
cat("\nsamps not filtered \n")
print(colSums(spe))

# end chunk10


# chunk11: filter according to read numbers; lakes first

collast <- ncol(spe)
rlast <- nrow(spe)

spe <- t(spe)
filter_samp <- rowSums(spe)>lakesumlevel
spe <- subset(spe,filter_samp)
spe <- t(spe)
filter_otus <- rowSums(spe)>otusumlevel
spe <- subset(spe,filter_otus)
spe_seqid_seqs <- subset(spe_seqid_seqs,filter_otus)

taxlast <- nrow(spe)
sitelast <- ncol(spe)

tax <- select(spe_seqid_seqs,taxonomy)

cat("\nsamp_reads filtered\n")
print(colSums(spe))

cat("\nNumber of OTUs  left after  OTUSUM filtering: ",taxlast," of original",rlast,"\n")
cat("Number of sites left after LAKESUM filtering: ",sitelast, "of original",collast,"\n")

# end chunk11


# ********** data structure **************

cat("\nmost frequent OTU?s (No of reads):\n")
print(head(rowSums(spe)))

cat("\nmin reads: ",min(colSums(spe)),"\n")

spepa <- decostand(spe,"pa")
spepasum <- rowSums(spepa)      # number of sites in which the OTUs (rows) occur...
sitecounts <- function(i,spepasum){sum(as.integer(spepasum==i))}
sitecounts10 <- sapply(1:10,sitecounts,spepasum=spepasum)
cat("\ncommunity structure:\nhow many OTU?s restricted to 1..10 sites: ",sitecounts10,"\n")

rlevel <- min(colSums(spe))
source("rrare100.R")
spe <- rrarefy100(t(spe),rlevel)
if (hel==T) spe <- decostand(spe,"hel")

# chunk12 filter OTUs of specific taxonomic groups if respective excl_flags are set TRUE

spe_excl <- cbind(tax, t(spe))

if (excl1_flag)
{  exclusivetaxindex <- grep(exclusive_group1,spe_excl$taxonomy)
spe_exclusive_group1 <- spe_excl[exclusivetaxindex,]  }

if (excl2_flag)
{  exclusivetaxindex <- grep(exclusive_group2,spe_excl$taxonomy)
spe_exclusive_group2 <- spe_excl[exclusivetaxindex,]  }

if (excl3_flag)
{  exclusivetaxindex <- grep(exclusive_group3,spe_excl$taxonomy)
spe_exclusive_group3 <- spe_excl[exclusivetaxindex,]  }

if (excl1_flag==T | excl2_flag==T | excl3_flag==T) {exclusive_groups <- rbind(if (excl1_flag==T) {spe_exclusive_group1}, if (excl2_flag==T) {spe_exclusive_group2}, if (excl3_flag==T) {spe_exclusive_group3})}
if (excl1_flag==T | excl2_flag==T | excl3_flag==T) {exclusive_groups_tax <- exclusive_groups$taxonomy}
#exclusive_groups <- select(exclusive_groups, -taxonomy)

if (write_exclusive_file) 
{
  write_name <- paste("just_", if (excl1_flag==T) {paste0(exclusive_group1, "_")}, if (excl2_flag==T) {paste0(exclusive_group2, "_")}, if (excl3_flag==T) {paste0(exclusive_group3, "_")},  "of_dataset_", filename, "_lakesumlevel_bigger_than_", lakesumlevel, "_otusumlevel_bigger_than_", otusumlevel, sep="")
  write.csv2(exclusive_groups, file=paste(write_name, "_rrarefied", ".csv", sep=""), row.names=F)
  write.csv2(exclusive_groups_tax, file=paste(write_name, "_rrarefied", ".tax", sep=""), row.names=F)
  
}
# end write excl file
