# chunk0: filenames etc
library(vegan)
library(dplyr)
library(data.table)
workdir <- "M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen/"
filename <- "final_lakes"
# end filenames

# chunk1: some variables
otusumlevel <- 2
lakesumlevel <- 100
write_IND <- T
taxsplit_string <- "; "
tax_cuts <- 1  # only values 1 or 2 are reasonable
# end default variables


# chunk2: data entry
# first: file-names and work-directory
data.entry(filename,workdir)
spefilename <- paste(filename,".prn",sep="")
spefiletax <- paste(filename,"_TaxPrnInd_ncbi_lakesumlevel_bigger_than_", lakesumlevel, "_otusumlevel_bigger_than_", otusumlevel, ".tax",sep="") #for TAX information ...
spefileprn <- paste(filename,"_TaxPrnInd_ncbi_lakesumlevel_bigger_than_", lakesumlevel, "_otusumlevel_bigger_than_", otusumlevel, ".prn",sep="") #for clean abundance file
indfilename <- paste(filename,"_TaxPrnInd_ncbi_lakesumlevel_bigger_than_", lakesumlevel, "_otusumlevel_bigger_than_", otusumlevel, ".ind",sep="") # for indicspec.R
# second: variables, here mainly for filtering
data.entry(otusumlevel,lakesumlevel,write_IND,tax_cuts,taxsplit_string)
# end data entry


# chunk3: read data
setwd(workdir)
TI <- Sys.time()
spe <- read.csv(spefilename, nrows=12000, sep="\t")
# end read data


# chunk4: clean taxonomy information
tax <- read.csv("kingdom_Euk_without_Metazoa_and_Embryophyta_of_dataset_best151102-103935-final_swoncbi85_tax_(AUSGANGSDATEN!)_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2.tax", nrows=12000, sep="\t")
tax <- tax$taxonomy
tax <- gsub("environmental samples","envi_samp",tax)
tax <- gsub("unclassified sequences","unclassified",tax)
tax <- gsub("uncultured organims","uncultured",tax)
# end clean tax info


orig_cols <- ncol(spe) # original number of samples
orig_otus <- nrow(spe) # original number of otus


# chunk6: filter according to number of OTU_reads AND number of sample_reads
# first filter samples
spe <- t(spe) # spe: only integer matrix with reads left according to chunk5
spe <- subset(spe,rowSums(spe)>lakesumlevel) # samples with low reads vanish!
# second filter otus
spe <- t(spe)
filterotus <- rowSums(spe)>otusumlevel
spe <- subset(spe,filterotus)
tax <- subset(tax,filterotus)
# end filter samples and otus


# chunk7: create NAs in tax dataframe, if no BLAST result was found 
tax <- as.data.frame(sapply(tax,function(x) if(x!="") return(x) else return(NA)))
#colnames(tax)[1] <- "taxonomy"
# end create NAs


# chunk8: print controls
otulast <- nrow(spe)
sitelast <- ncol(spe)
cat("Number of OTUs  left after  OTUSUM filtering: ",otulast," of original",orig_otus,"\n")
cat("Number of sites left after LAKESUM filtering: ",sitelast, "of original",orig_cols,"\n")
# end print controls



# chunk10: if not already ordered, sort otus according to read numbers
rs <- rowSums(spe)
indrs <- order(rs,decreasing=T)
spe <- spe[indrs, ]
tax <- as.data.frame(tax[indrs, ])
colnames(tax)[1] <- "taxonomy"
rs <- rs[indrs]
# end sorting



# chunk11: write data to .prn and to .tax files: sites = colnames, taxa = rownames...
n<-nrow(spe)
rn <- paste("N",1:n,sep="") # =row.names(spe)
csv <- cbind(rn,spe)  # this write technique is compatible with data.frames AND data.tables
row.names(tax)<- rn     # different write technique, whether it is faster remains to be eludicated
write.table(csv,file=spefileprn,row.names=F,sep="\t")
write.table(tax,file=spefiletax,col.names=NA,sep="\t")
# end write data to files for further processing (taxonomy and statistics)




# chunk12: tshort = short tax info, only last 2 items; only needed for indicative species script
if (write_IND)
{ lastitems <- function(x) 
                 { xlen <- length(x) 
                   if (xlen==1 | tax_cuts==1) 
                      {x<-x[xlen]
                       return(x) } else
                   x <- paste(x[xlen-1],x[xlen],sep="")
                   return(x) 
                 }
  taxshort <- strsplit(as.character(tax$taxonomy),taxsplit_string)
  taxshort <- sapply(taxshort,lastitems)
}
# end shorten taxonomy


# chunk13: write abundance file with shortened tax-info in first column
if (write_IND)
{ taxshort <- paste(rn,"_",rs,"_",taxshort,sep="")
  spe_ind <- cbind(taxshort,spe)
  write.table(spe_ind,file=indfilename,row.names=F, col.names=T, dec=".",sep="\t", quote=FALSE)
} 
# end write indicative species file


################################################################################
# now some infos about data structure (no writing of files)
################################################################################


# chunk14: a posteriori sort according to colSums, i.e. reads in samples; print controls
cs <- colSums(spe)
indcs <- order(cs,decreasing=T)
spe <- spe[,indcs]
cat("No of reads in each site:\n")
print(colSums(spe))
cat("most frequent OTUs (No of reads):\n")
print(head(rs))  # = print(head(rowSums(spe)))
# end sorting 


# chunk15: structure: how many OTUs only in 1,2..10 lakes
spepa <- decostand(spe,"pa")
spepasum <- rowSums(spepa)      # number of sites in which the OTUs (rows) occur...
sitecounts <- function(i,spepasum){sum(as.integer(spepasum==i))}
sitecounts10 <- sapply(1:10,sitecounts,spepasum=spepasum)
cat("community structure:\nhow many OTUs restricted to 1..10 sites: ",sitecounts10,"\n")
# end structure info


# chunk16: info about elapsed time (should be some seconds at most)
cat("elapsed runtime: ")
print(Sys.time()-TI) # show runtime
# end time info
