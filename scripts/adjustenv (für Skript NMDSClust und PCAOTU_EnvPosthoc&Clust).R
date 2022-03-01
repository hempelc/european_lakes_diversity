#********************************************************************************************
# ADJUSTMENT OF SITE NAMES of a long environmental file 
# (with more sites than the species file, may even be wrongly ordered) 
# shorten it and put the sites into the same order as the species file
# envnames to be ordered must BEGIN with the spenames (identical character succession)
# Author: Manfred Jensen, biodiversity group, University of Duisburg-Essen, Germany
# V1.1 December 9, 2013   email: manfred.jensen@uni-due.de
#********************************************************************************************

adjustenv <- function(x,y)    # x <- spe , y <- env
{ 
nspe <- nrow(x)
spenames <- row.names(x)
envnames <- row.names(y)
envindices <- pmatch(spenames,envnames, nomatch=0, duplicates.ok=F)
if (0 %in% envindices) 
{print("error; site names not adjustable (incompatible)")
 stop
}
env2 <- y[0,]
for (i in 1:nspe) {env2 <- rbind(env2,y[envindices[i],])}
y  <- env2
return(y)
}

#test:
#csvsepchar <-";" # english spreadsheet: use ","
#spefilename <- "Austria2006Euk2.prn"
#envfilename <- "Austria2006ENV.csv"
#spe <- read.csv(spefilename, row.names=1,nrows=500,sep="\t") 
#env <- read.csv(envfilename, row.names=1,nrows=500,sep=csvsepchar) #longer file, wrong order
#env <- adjustenv(spe,env) # now spe and env match ...
