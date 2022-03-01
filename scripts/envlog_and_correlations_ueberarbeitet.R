# Modify environmental data and plot correlations
# V1.4 Author M.Jensen, May 25, 2016
# Biodiversity Group University of Duisburg-Essen, Germany
# based on Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
# select variables of interest in line 90,91
#**************************************************************************************

# chunk0: default filenames etc
library(gclus)
library(data.table)
library(dplyr)
workdir <- "M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen"
envname <- "lakes_abiotic_data_just_numbers_to_work_with.csv" # name of environmental file
# end chunk 0 (files must be in the working directory)



# chunk1: some variables
decpoint <- ","
splitchar <- "_"
#short_flag <- T
corr_flag <- T
winsys <- T
# end variables



# chunk2: data entry: modification of default values
# first: file-names and work-directory
data.entry(workdir,envname,decpoint,corr_flag,splitchar,winsys)
setwd(workdir)
# end data.entry



# chunk 3: read env file
env <- fread(envname,dec=decpoint, data.table=F)
rownames(env) <- env[,1]
env <- env[,-1]
# end chunk 3



# chunk 3a: shorten object = row names
#if (short_flag) {
#shorts <- strsplit(as.character(rownames(env)),splitchar)
#rownames(env) <- sapply(shorts,function(x) {xlen <- length(x) 
#                                if (xlen==1) x<-x[xlen] else
#                                    x <- paste(x[xlen-1],x[xlen],sep="")
#                                return(x) } ) }
# end chunk 3a


# chunk3b NAs are not allowed, eliminate the samples with NAs 
#env <- t(env)
#env <- select(env,-AH33,-HC33,-HE33)
#env <- as.data.frame(t(env))
# end chunk3b



# chunk 4: log nolog etc have to be modified depending on the env-file
View(head(env))
nologs <- c("Höhe [m]","pH Durchschnitt", "T Durchschnitt")
cat("entry of env no_log_variables: no logarithm\n")
data.entry(nologs)
envnolog <- env[,nologs]
inds <- match(nologs,colnames(env),nomatch=0)
envlog <- env[, -inds]
envlog <- log1p(envlog)
env   <- cbind(envnolog,envlog)
# special additional env_vars possible
#$TP.DN <- env$TP / env$DN
#env$DN.DOC <- env$DN / env$DOC
# end chunk 4 log - nolog - sort



# chunk 5: write envlog file
envlogname <- paste(envname, "_log")
cat("writing",envlogname,"\n\n")
write.csv2(env,file=envlogname)
# end chunk 5







# R-mode correlation matrices
# ***************************


# chunk 6: select variables of interest
if (corr_flag)
{ View(head(env))
  voi <- c(1:5)
  cat("Select Variables of Interest")
  data.entry(voi)
  env <- env[,voi]         # if area with NA = missing values, select only interesting vars
# end chunk6





# chunk 7: Pearson r linear correlation among environmental variables
# function pairs etc
source("panelutils_pvalue (für Skript envlog).R")
env.pearson <- cor(env)	# default method = "pearson"
# Reorder the variables prior to plotting
env.o <- order.single(env.pearson)
if (winsys) windows(title="Linear correlation matrix", 10, 10)
op <- par(mfrow=c(1,1), pty="s")
pairs(env[,env.o], lower.panel=panel.smooth, upper.panel=panel.cor,
      diag.panel=panel.hist, main="Pearson Correlation Matrix")
par(op)
# end chunk7


# chunk 8: Kendall tau rank correlation among environmental variables
env.ken <- cor(env, method="kendall")
env.o <- order.single(env.ken)
if (winsys) windows(title="Rank correlation matrix Kendall", 10, 10)
op <- par(mfrow=c(1,1), pty="s")
pairs(env[,env.o], lower.panel=panel.smoothb, 
      upper.panel=panel.cor, no.col=FALSE,
      method="kendall", diag.panel=panel.hist, 
      main="Kendall Correlation Matrix")
par(op)
# end chunk 8


# chunk 9: Spearman rank correlation among environmental variables
env.spear <- cor(env, method="spearman")
env.o <- order.single(env.spear)
if (winsys) windows(title="Rank correlation matrix Spearman", 10, 10)
op <- par(mfrow=c(1,1), pty="s")
pairs(env[,env.o], lower.panel=panel.smoothb, 
      upper.panel=panel.cor, no.col=FALSE,
      method="spearman", diag.panel=panel.hist, 
      main="Spearman Correlation Matrix")
par(op)
# end Spearman

# chunk 10: reorder env for partial correlation; x in col1, y in col2, z in col3:ncol
corvarx_colnum <- 1
corvary_colnum <- 2
data.entry(corvarx_colnum,corvary_colnum)
if (corvarx_colnum!=1 | corvary_colnum!=2) 
                        {envy <- env[,corvary_colnum]
                         envx <- env[,corvarx_colnum]
                         envxname <- colnames(env)[corvarx_colnum]
                         envyname <- colnames(env)[corvary_colnum]
                         env  <- env[,-c(corvarx_colnum,corvary_colnum)]
                         env  <- cbind(envx,envy,env)
                         colnames(env)[1] <- envxname
                         colnames(env)[2] <- envyname}
# end chunk 10


# chunk 11: example for partial correlation
library(ppcor)
cat("Partial Correlation: corvarx, corvarx in col1 col2\n")
cat("Correlation Vars = ",colnames(env)[1]," versus ",colnames(env)[2],"\n")
cat("Elimination of other Vars\n")
for (i in 3:ncol(env))
{ partial_alt_Temp_minus_Var <- pcor.test(env[,1],env[,2],env[,i],
                                          method="pearson")
  cat("Var = ",colnames(env)[i],"\n")
  print(partial_alt_Temp_minus_Var)
}
# end chunk 11


detach("package:ppcor")
detach("package:gclus")
detach("package:data.table")
detach("package:dplyr")
}

# end script envlog.R