setwd("M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen")  # workfolder for your data
csvsepchar <- "\t"   # try "," for english version of EXcEL ...
library(vegan)
library(compiler)
library(dplyr)
library(data.table)


filename <- "151102-103935-final" 

#....filtering of rows: elimnate rows with too low  number of reads of the corresponding OTU
#############################################################################################
#############################################################################################
# otusumlevel: filtering of OTU?s with too low number of reads;original OTU numbers survive #
#*******************************************************************************************#
otusumlevel <- 0
lakesumlevel <- 0
write_IND = 0
data.entry(filename,otusumlevel,lakesumlevel)
#############################################################################################
# lakesumlevel: filtering of sites with too low number of reads
#********************************************************************************************
TI <- Sys.time()
#############################################################################################
# suffix of site - names, should be eliminated, different for bacteria !!!!
#********************************************************************************************
suffixsitename <-""


# START of CONVERSION CODE
#********************************************************************************************

csvfiletax <- paste(filename,".tax",sep="") #for TAX information ...
indfilename <- paste(filename,".ind",sep="") # for indicspec.R

csvfile <- paste(filename,".csv",sep="")    #input-file-name automatically


# filtering of lakes and OTUS according to lakesumlevel and otusumlevel

######################################################################################################

csvmatrix <- fread(csvfile)
csvmatrix <- as.data.frame(csvmatrix)


MVfilename <- paste(filename,".prn",sep="") #for filtered numeric count data

tax <- csvmatrix$complete.taxonomy
filenamelen <- nchar(filename)
filenameend <- substring(filename,filenamelen-2,filenamelen)

csvmatrix <- select(csvmatrix, -sequence.length,-sequences,-c(sum:complete.taxonomy))
  
collast <- ncol(csvmatrix)
rlast <- nrow(csvmatrix)
csvmatrix <- t(csvmatrix)
csvmatrix <- subset(csvmatrix,rowSums(csvmatrix)>lakesumlevel)
csvmatrix <- t(csvmatrix)
filterotus <- rowSums(csvmatrix)>otusumlevel
csvmatrix <- subset(csvmatrix,filterotus)
tax <- subset(tax,filterotus)


# filtering finished already !! ######################################################################

# print controls  ***********************************************************

taxlast <- nrow(csvmatrix)
sitelast <- ncol(csvmatrix)

cat("Number of OTUs  left after  OTUSUM filtering: ",taxlast," of original",rlast,"\n")
cat("Number of sites left after LAKESUM filtering: ",sitelast, "of original",collast,"\n")


# ***************************************************************************

# write data to .prn and to .tax files: sites = colnames, taxa = rownames...

#****************************************************************************************************

# use transponed lakedatamatrix; write/read much faster !
n<-nrow(csvmatrix)
rn <- paste("N",1:n,sep="")#row.names(csvmatrix)

csv <- cbind(rn,csvmatrix)
tax <- cbind(rn,tax)
write.table(csv,file=MVfilename,row.names=F,sep="\t")
write.table(tax,file=csvfiletax,row.names=F,sep="\t")


#**************************************************************************************************

# ********** sort according to rowSums, i.e. frequency of OTUs **************

rs <- rowSums(csvmatrix)
ind <- order(rs,decreasing=T)
csvmatrix <- csvmatrix[ind,]

cs <- colSums(csvmatrix)
ind <- order(cs,decreasing=T)
csvmatrix <- csvmatrix[,ind]
#cat("No of reads in each site:\n")
#print(colSums(csvmatrix))
cat("most frequent OTU?s (No of reads):\n")
print(head(rowSums(csvmatrix)))


# ********** end sorting ****************************************************

cat("min reads: ",min(colSums(csvmatrix)),"\n")


# ********** structure: how many OTUs only in 1,2..10 lakes ****************

spepa <- decostand(csvmatrix,"pa")
spepasum <- rowSums(spepa)      # number of sites in which the OTUs (rows) occur...
sitecounts <- function(i,spepasum){sum(as.integer(spepasum==i))}
sitecounts10 <- sapply(1:10,sitecounts,spepasum=spepasum)
cat("community structure:\nhow many OTU?s restricted to 1..10 sites: ",sitecounts10,"\n")


# ********** end structure: how many OTU?s only in 1,2..10 lakes ****************

cat("runtime: ")
print(Sys.time()-TI) # show runtime
