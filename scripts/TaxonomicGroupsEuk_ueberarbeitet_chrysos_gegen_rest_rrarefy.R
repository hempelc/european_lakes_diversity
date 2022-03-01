#####################################################################################################
# Script: TaxonomicGroups_Boenigk&Jensen.R
# Author Manfred Jensen,Universität-Duisburg-Essen,Germany Version 4.2 Mar 5, 2015: means, sd    ###
# V4.3 compiler for necessary crucial double loops 233-240 Nov 11, 2015
# new source in wd: extraspe.R & adjustenv.R;  To be modified for your personal examination:
# LINES 10,11,16,21-28;32-40("#" to comment out if abiotic envi or sort not neccessary)  
# line 160,161:groups assigned; G0,GRest: possibly don`t appear within the plots
# CHARTS start: line 280; V2.9:preserve colors; V3.6 NEW sort on subgroups:line 28;32,33:external sort
# .prn and .tax files = output of Script "XLSTaxPrn.R", .agg files (tax groups) for further use
##*****************************************************************************************************
library(compiler)
enableJIT(3)

workfolder <- "M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen/"  # set working directory
#cat("OPEN Data Editor","\n","Click on task bar symbol","\n")
#data.entry(workfolder)
setwd(workfolder)

library(vegan)
library(reshape)
library(dplyr)
source("extraspe (für Skript TaxonomicGroupsEuk).R")
filename <- "kingdom_Euk_without_Metazoa_and_Embryophyta_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2"
outfilename <- "just_Chrysophyceae_Synurophyceae_Oikomonadaceae_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2"

rlevel <- F # if rlevel=F, rlevel is the minimum of readsums of lakes
revord <- F # flag: stacking of bars inverted
Inclunknown <- T      # GREST = up to now unknown sequence "unknown affiliation"
Inclothers <- T       # G0 = others: blast with a result beyond our groups
Sortsubs <- T
sortglob <- T
#************************************************************************************** #
# modify defaults (a posteriori); closing data editor = data entry finished             #
#************************************************************************************** #
cat("OPEN Data Editor","\n","(modify default values in first row, if necessary)","\n",
    "revord: 0 NO reverse sorting order of sites, 1 reverse [vertical bars: contrary]","\n",
    "Inclunknown:  0 for excluding OTUs without reasonable blast result, 1 otherwise ","\n",
    "Inclothers: 0 leave out blasts beyond our specified groups, 1 otherwise","\n",
    "Sortsubs: 0 sort according to species richness, 1 additional subgroup sorting","\n",
    "sortglob: 0 no sorting if no subgroups, 1 sorting globally if no subgroups","\n",
    "accept values by *** CLOSING the Data Editor *** press X at upper right corner of window","\n\n")
data.entry(filename,revord,Inclunknown,Inclothers,Sortsubs,sortglob, rlevel)
#****************** end of modify map options ********************************************************#
if (Inclunknown==T) {GREST <- "unknown affiliation"} else {GREST <- "N"}
if (Inclothers==T) {G0 <- "others"} else {G0 <- "N"}
#gi <- c(0,81,95,111,122,133,141,162,167,202,232) #indices for subclusters "handish"
#ginames <- c("Austria","Lake1","Colne","Whale","Biofilm","Borehole","Lake2","winterpico","Soil","HNF")
#gi <- c(0,16,52,61,138,147,176,192)
#ginames <- c("AuSeas","AuSpat","Baikal","ExpBac","LakeD","Soil","WhaleF")

spefilename <- paste(filename,".csv",sep="")
taxfilename <- paste(filename,".tax",sep="")
aggfilename <- paste(outfilename,".agg",sep="")


spe <- fread(spefilename, data.table=F)
tax <- read.csv(taxfilename, sep="\t")


#Umsortieren Eukarya
spe.A111AU <- select(spe,A111AU_B)
spe.A132OS <- select(spe,A132OS_C)
spe.A152WI <- select(spe,A152WI_B)
spe.N121SA <- select(spe,N121SA_C)
spe.N261LU <- select(spe,N261LU_C)
spe.O022TU <- select(spe,O022TU_A)
spe.O111BA <- select(spe,O111BA_C)
spe.O121RA <- select(spe,O121RA_C)
spe.O241PL <- select(spe,O241PL_C)
spe.S022PA <- select(spe,S022PA_A)
spe.S031BU <- select(spe,S031BU_C)
spe.S102LR <- select(spe,S102LR_A)
spe.S153TR <- select(spe,S153TR_C)
spe.S171MA <- select(spe,S171MA_A)
spe.S271TO <- select(spe,S271TO_C)
spe.S281PM <- select(spe,S281PM_C)
spe.S301MM <- select(spe,S301MM_C)
spe.Z042CP <- select(spe,Z042CP_C)
spe.Z071SI <- select(spe,Z071SI_A)
spe.Z111VV <- select(spe,Z111VV_C)
spe.Z122OU <- select(spe,Z122OU_C)
ginames <- c("A111AU","A132OS","A152WI","N121SA","N261LU","O022TU",
             "O111BA","O121RA","O241PL","S022PA","S031BU","S102LR",
             "S153TR","S171MA","S271TO","S281PM","S301MM","Z042CP","Z071SI","Z111VV","Z122OU")
gi <- c(0:21)
spe <- t(spe)
if (rlevel==F) rlevel <- min(rowSums(spe)) else rlevel=rlevel
if (Sortsubs==T) 
{ cat("OPEN Data Editor","\n","modify col 1+2: Subclust-Ranges,1:indices [first=0] 2:names","\n",
      "3: mean/var files: 1 English  2 German Excel","\n\n")
  csvsepchar <- "\t"    # english file format be aware of env file format ...
  English_German <- 1
  data.entry(gi,ginames,English_German)
  if (English_German==2) csvsepchar <- ";"
}

#speAU <- spe[grep("AU",rownames(spe)),]
#speWA <- spe[grep("WA",rownames(spe)),]
#speFu <- spe[grep("FU",rownames(spe)),]
#speFU <- spe[grep("Fu",rownames(spe)),]

#spe <- rbind(speAU,speWA,speFu,speFU)

#source("subgroups7Large.R")        # subclusters from external script; also used by ClusterNumSil.R
#source("taxnewclustgroups.R")      # subclusters from external script
if (revord) spe <- spe[rev(1:nrow(spe)),]    # reverse order of sites
if (nrow(spe)>60) {labsize <- 0.4} else labsize <- 0.8   #???
#env <- read.csv("AlpineSpatialEnv.env", row.names=1, nrows=500, sep="\t")
#env <- adjustenv(spe,env) ???
# sort on an environmental factor
#spe <- spe[order(env[,'Temp'],decreasing=F),] # sorting of sites on environ factor
# end sort


source("rrare100.R")
sperrare <- rrarefy100(spe,rlevel)
#sperrare <- extrasperrare(spe,sperrare,rlevel,chao=F)

spe_rich      <- rowSums(sperrare)  # ALMOST identical to rarefy(spe,rlevel), but extraspe


colorvec <- c("red","darkgreen","darkblue","cyan","darkred",
              "magenta","yellow","lightgreen","orange","lightblue","grey","black",
              "darkorange","lightyellow","brown","pink","green",
              "lightcyan","darkcyan","blue","lightblue","lightgrey",
              "darkmagenta","lightpink","white","darkgrey","beige","darkgreen",
              "turquoise","gray","white","purple","bisque")



taxkey <- c("Ciliophora","Dinophyceae","Apicomplexa","Alveolata",
            
            "Bicosoecida","Oomycetes","Chrysophyceae","Synurales","Synurophyceae",
            "Oikomonadaceae","Bacillariophyta","Stramenopiles",
            
            "Cercozoa","Rhizaria",
            
            "Heterolobosea","Euglenida","Kinetoplastida","Euglenozoa",
            "Diplomonadida","Dysnectes","Jakobida","Kipferlia","Oxymonadida",
            "Parabasalia","Retortomonadidae","Trimastix",
            
            "Choanoflagellida","Chytridiomycota",
            "Microsporidia",
            "Ascomycota","Basidiomycota","Glomeromycota","Fungi",
            
            "Amoebozoa",
            
            "Apusozoa",
            
            "Glaucocystophyceae","Rhodophyta","Embryophyta","Viridiplantae",
            
            "Cryptophyta","Katablepharidophyta",
            
            "Haptophyceae","Eccrinales","Breviata","Centroheliozoa","Ichthyosporea",
            "Corallochytrium","Nucleariidae","Telonemida","Metazoa",
            
            "Eukaryota","envi_samp","unclassified",NA)


G1 <- "Chrysophyceae"  
G2 <- "Others"

G3 <- "A3-Apicomplexa"
G4 <- "A4-Alveolata-rest"

G5 <- "B1-Bicosoecida"
G6 <- "B2-Oomycetes"
G7 <- "B3-Chrysophyceae"
G8 <- "B4-Synurophyceae"
G9 <- "BE-Oikomonadaceae"
G10<- "B5-Bacillariophyta"
G11<- "B6-Stramenopiles-rest"

G12<- "C1-Cercozoa"
G13<- "C2-Rhizaria-rest"

G14<- "D1-Heterolobosea"
G15<- "D2-Euglenida"
G16<- "D3-Kinetoplastida"
G17<- "D4-Euglenozoa-rest"
G18<- "D5-Excavata-rest"

G19<- "E1-Choanoflagellida"
G20<- "E2-Chytridiomycota"
G21<- "E3-Microsporidia"
G22<- "E4-Ascomycota"
G23<- "E5-Basidiomycota"
G24<- "E6-Glomeromycota"

G25<- "E7-Amoebozoa"
G26<- "E8-Apusozoa"

G27<- "F1-Glaucocystoph."
G28<- "F2-Rhodophyta"
G29<- "F3-Viridipl. no Embr"

G30<- "G1-Cryptophyta"
G31<- "G2-Katablepharidoph"
G32<- "H1-Haptophyceae"



groupvec <- c(G1,G2,G3,G4,G5,G6,G7,G8,G10,G11,G12,G13,G14,G15,G16,G17,G18,G19,G20,G21,
              G22,G23,G24,G25,G26,G27,G28,G29,G30,G31,G32,G0,GREST)      # no G9 ... 33 colors needed

# line 245-248: final determination of colors
# here: max. 33 groups and 33 colors defined


taxkey2 <- c(G2,G2,G2,G2,
             G2,G2,G1,G2,G1,G1,G2,G2,
             G2,G2,
             G2,G2,G2,G2,G2,G2,G2,G2,G2,G2,G2,G2,
             G2,G2,G2,G2,G2,G2,G2,
             G2,
             G2,
             G2,G2,G2,G2,
             G2,G2,
             G2,G2,G2,G2,G2,G2,G2,G2,G2,
             G2,G2,G2,G2)
colorvec <- colorvec[1:length(groupvec)]




taxkeylength <- length(taxkey)
taxlastrow <- nrow(tax)
taxlist <- list("ZZ")

taxsplit <- apply(tax,2,strsplit,";")

for (i in 1:taxlastrow)
{ taxlist[[c(1,i)]] <- "ZZ"
  for (j in 1:taxkeylength)  
    if (taxkey[j] %in% taxsplit[[c(1,i)]]) 
        { taxlist[[c(1,i)]]<-taxkey2[j]
          break
        }
}
totalreads <- sum(rowSums(spe))
averagereads <- totalreads/nrow(spe)
cat("total number of taxa : ",length(taxlist[[c(1)]])," OTU´s","\n",sep="\t")
cat("total number of reads: ",totalreads,"\t#of sites:",nrow(spe), 
    "   average reads/site:",averagereads,sep="\t","\n\n")

sperrare.trans <- t(sperrare)
sperrare.trans <- subset(sperrare.trans,taxlist[[1]] != "N")
spe.trans <- t(spe)
spe.trans <- subset(spe.trans,taxlist[[1]] != "N")

spe.norm <- decostand(t(spe.trans), "hellinger",na.rm=T)
spe.norm.trans <- t(spe.norm)

sperrare <- t(sperrare.trans)
spe <- t(spe.trans)
spe_rich <- rowSums(sperrare)
sri <- spe_rich  # temporal sorting vector on species richness


# calculate spe_otu_rich, i.e. bar length = OTU_richness, but subdivisions according to OTU numbers
spe_otu_pa    <- decostand(spe,method="pa")             # presence /absence data  : ok
spe_num       <- rowSums(spe_otu_pa)                    # specnumbers == specnumber(spe) : ok
spe_otu_rich  <- sweep(spe_otu_pa,1,spe_num,"/")        # -> rowSums = 1 in every row: ok 
spe_otu_rich  <- sweep(spe_otu_rich,1,spe_rich,"*")     # -> rowSums = spe_rich in every row: ok
# end calculate spe_otu_rich



# new matrix with sums = richness, but numbers according to OTU numbers
spe_otu_rich.trans <- t(spe_otu_rich)
#spe_otu_rich.trans <- subset(spe_otu_rich.trans,taxlist[[1]] != "N")
# end new matrix


# lines for sorting on subgroup species richness: producing sorting vector
#*************************************************************************
if (Sortsubs==T)
{
  if (revord) {gi <- gi[length(gi)] - rev(gi)
               ginames <- rev(ginames)}
  gilen <- length(gi)-1
  sri <- spe_rich[0]
  for (n in 1:gilen)
  { i <- gi[n]+1
    j <- gi[n+1]
    subg <- spe_rich[c(i:j)]  # subvector (group)
    subg <- subg[order(subg,decreasing=T)] # sort subvector according to species richness
    sri  <- c(sri,subg)   # combine subvectors to a large sort vector  
  } # end producing subgroup sorting vector
} else {if (sortglob) sri<-sort(sri,decreasing=T)}  #force sorting on species richness


# sort on [subgroup] species richness
spe_rich      <- sri    # rebuild of species richness vector in the new order
sri           <- as.data.frame(sri)
index         <- match(row.names(spe),row.names(sri))
sperrare      <- sperrare[order(index),]
spe           <- spe[order(index),]
spe_otu_rich.trans <- spe_otu_rich.trans[,order(index)]
sperrare.trans <- sperrare.trans[,order(index)]
spe.trans <- spe.trans[,order(index)]
spe.norm.trans <- spe.norm.trans[,order(index)]
# end of sort


#print(spe_rich)

taxlist[[1]] <- subset(taxlist[[1]],taxlist[[1]] != "N")


if (Inclunknown==F & Inclothers==F)
            write.table(spe.trans,paste(outfilename,".pNO",sep=""),row.names=T, col.names=T, 
            dec=".",sep="\t", quote=FALSE)

t.sperrare.agg <- aggregate(sperrare.trans,taxlist,sum)
t.spe.agg <- aggregate(spe.trans,taxlist,sum)
t.spe.otu <- aggregate(spe.trans>0,taxlist,sum)
t.spe.norm.agg <- aggregate(spe.norm.trans,taxlist,sum)

# new matrix with sums = richness, but numbers according to OTU numbers
t.spe.otu_rich <- aggregate(spe_otu_rich.trans,taxlist,sum)
spe_otu_rich.agg <- t(t.spe.otu_rich)
aggofilename <- paste(aggfilename,"o",sep="")
write.table(spe_otu_rich.agg,aggofilename,row.names=T, col.names=F, 
            dec=".",sep="\t", quote=FALSE)
spe_otu_rich.agg <- read.csv(aggofilename,row.names=1, nrows=500, sep="\t")
# end new matrix

speotu.agg <- t(t.spe.otu)
write.table(speotu.agg,file=aggfilename,row.names=T, col.names=F, 
            dec=".",sep="\t", quote=FALSE)
speotu.agg <- read.csv(aggfilename,row.names=1, nrows=500, sep="\t")

spe.agg <- t(t.spe.agg)
write.table(spe.agg,file=aggfilename,row.names=T, col.names=F, 
            dec=".",sep="\t", quote=FALSE)
spe.agg <- read.csv(aggfilename,row.names=1, nrows=500, sep="\t") 

spe.norm.agg <- t(t.spe.norm.agg)
write.table(spe.norm.agg,file=aggfilename,row.names=T, col.names=F, 
            dec=".",sep="\t", quote=FALSE)
spe.norm.agg <- read.csv(aggfilename,row.names=1, nrows=500, sep="\t") 

sperrare.agg <- t(t.sperrare.agg)
write.table(sperrare.agg,file=aggfilename,row.names=T, col.names=F, 
            dec=".",sep="\t", quote=FALSE)
sperrare.agg <- read.csv(aggfilename,row.names=1, nrows=500, sep="\t") 


# ************************************************* calculate subgroup means or SD ********
prop <- function(speagg,subgi,subginames,statflag)
{
  speprop <- prop.table(as.matrix(speagg),1)
  rownumbers <- length(subginames)
  readsagg <- speprop[0,]
  n_subg <- rep(0,rownumbers)
  for (i in 1:rownumbers)
  { rangeA <- subgi[i] + 1
    rangeB <- subgi[i+1]
    if (statflag==1){subgroupv <- t(colMeans(speprop[c(rangeA : rangeB),]))}
      else 
      { colSD <- apply(speprop[c(rangeA:rangeB),],2,sd)   
        subgroupv <- t(colSD)
      }
    readsagg <- rbind(readsagg,subgroupv)
    n_subg[i] <- rangeB - rangeA + 1
  }
  readsagg <- cbind(readsagg,n_subg)
  row.names(readsagg) <- subginames
  readsagg <- t(readsagg)              # transponed excel file only this line ...
  sites <- row.names(readsagg)
  readsagg <- cbind(sites,readsagg)
  names(readsagg)[1] <- "Sites"
  return(readsagg)
}
# ********************************* end function prop: calculate means or SD **************


# ********************************* function write mean and SD to file **************************
writeprop <- function(aggmatrix,fileext)
{   agg         <- prop(aggmatrix,gi,ginames,1)
    aggfilename <- paste(filename,fileext,"_mean.csv",sep="")
    write.table(agg,file=aggfilename,row.names=F, col.names=T,sep=csvsepchar,quote=F)  
    agg         <- prop(aggmatrix,gi,ginames,2)
    aggfilename <- paste(filename,fileext,"_SD.csv",sep="")
    write.table(agg,file=aggfilename,row.names=F, col.names=T,sep=csvsepchar,quote=F)
}
# ********************************* end write mean and SD to file *******************************

# *********************** write subgroup means and SD to FILES ****************************
#if (Sortsubs==T)
#{ 
#  writeagg <- writeprop(spe.agg,"_reads")
#  writeagg <- writeprop(spe.norm.agg,"_hel")
#  writeagg <- writeprop(speotu.agg,"_OTU")
#  writeagg <- writeprop(sperrare.agg,"_RARE")
#}  
#************************ end write means and SD *****************************************





# **************************************************************
# preserve colors before start of plot  ************************
colorvector <- subset(colorvec,groupvec %in% t.spe.agg[,1]) #colors preserved for every tax group
groupindices <- match(groupvec,t.spe.agg[,1],nomatch=0)

# groups which are absent **********************************************************
cat("Absence of groups:", sep="\n")
R_ENABLE_JIT <- 0 # disable compiling
for (i in 1:length(groupvec)) {if (groupindices[i]==0) {cat(groupvec[i],sep="\n")}}
R_ENABLE_JIT <- 3 # enable compiling
# **********************************************************************************














#****************************************************************************************
#                         BARCHART 1 (species richness)
#****************************************************************************************

windows(title="Species Richness: megasystematic groups ",16,11)

#spe.table <- as.table(as.matrix(spe_otu_rich.agg))
spe.table <- as.table(as.matrix(sperrare.agg))
xmax <- max(rowSums(sperrare.agg))+5
xlimits <- c(0,xmax)

#xlabel <- paste("species richness and OTU-composition of TAX groups [rarefaction at reads= ", 
#                as.character(rlevel),"]", sep="")
xlabel <- paste("species richness and composition of TAX groups [by rarefaction at reads= ", 
                as.character(rlevel),"]", sep="")
ylabel <- "rarefacted No. of OTU´s"
textlegends <- as.vector(dimnames(spe.table)[[2]])
titlelegend <- "Colour legend"



specrich <- barchart(spe.table,xlab=ylabel,ylab=xlabel, horizontal=T,xlim=xlimits, col=colorvector,
                     scales=list(y=list(cex=labsize,rot=0, col=ifelse(row.names(spe.table) %in% c("A111AU_B", "A152WI_B", "S102LR_A", "S031BU_C", "O111BA_C", "S171MA_A", "S153TR_C", "Z042CP_C", "Z071SI_A", "Z122OU_C"), 2, 1)), 
                                 x=list(cex=0.8,rot=0)), key=list(col=colorvector,space="right",title=titlelegend,
                                                                  text=list(textlegends,cex=1,col="black"),rect=F))

plot(specrich)
#************************ END OF BARCHART ******* horizontal=F: rot=90, horizontal=T: rot=0 

#****************************************************************************************
#                         BARCHART 1a (species richness but relative)
#****************************************************************************************
windows(title="Proportions of groups ", 16, 11)

spe.table <- as.table(as.matrix(sperrare.agg))


xlimits <- c(0,1)
xlabel <- paste("taxonomic groups, fraction of species richness [rarefaction at reads= ",
                as.character(rlevel),"]", sep="")
ylabel <- "relative probabilities of groups"
textlegends <- as.vector(dimnames(spe.table)[[2]])

specprop <- barchart(prop.table(spe.table,margin=1),xlab=ylabel,ylab=xlabel, horizontal=T,xlim=xlimits, col=colorvector,
                     scales=list(y=list(cex=labsize,rot=0, col=ifelse(row.names(spe.table) %in% c("A111AU_B", "A152WI_B", "S102LR_A", "S031BU_C", "O111BA_C", "S171MA_A", "S153TR_C", "Z042CP_C", "Z071SI_A", "Z122OU_C"), 2, 1)),x=list(cex=0.8,rot=0)),
                     key=list(col=colorvector,space="right",title=titlelegend,
                              text=list(textlegends,cex=1,col="black"),rect=F))

plot(specprop)
#************************ END OF BARCHART ************************************************





#****************************************************************************************
#                         BARCHART2 (proportions reads)
#****************************************************************************************
windows(title="Proportions of groups ", 16, 11)

spe.table <- as.table(as.matrix(spe.agg))

xlimits <- c(0,1)
xlabel <- "proportion of taxonomic groups (based on number of Read´s)"
ylabel <- "relative No. of Read´s"
textlegends <- as.vector(dimnames(spe.table)[[2]])

plotreads <- barchart(prop.table(spe.table,margin=1),xlab=ylabel,ylab=xlabel, horizontal=T,xlim=xlimits, col=colorvector,
                      scales=list(y=list(cex=labsize,rot=0, col=ifelse(row.names(spe.table) %in% c("A111AU_B", "A152WI_B", "S102LR_A", "S031BU_C", "O111BA_C", "S171MA_A", "S153TR_C", "Z042CP_C", "Z071SI_A", "Z122OU_C"), 2, 1)),x=list(cex=0.8,rot=0)),
                      key=list(col=colorvector,space="right",title=titlelegend,
                               text=list(textlegends,cex=1,col="black"),rect=F))
plot(plotreads)
#************************ END OF BARCHART ************************************************


#****************************************************************************************
#                         BARCHART2a (proportions drare/hellinger)
#****************************************************************************************
windows(title="Proportions of groups ", 16, 11)

spe.table <- as.table(as.matrix(spe.norm.agg))

#xlimits <- c(0,max(rowSums(spe.table)))
xlimits <- c(0,1)
xlabel <- "proportion of taxonomic groups (based on number of Read´s / Hellinger)"
ylabel <- "relative No. of Read´s - Hellinger transformed"
textlegends <- as.vector(dimnames(spe.table)[[2]])

plotreads <- barchart(prop.table(spe.table,margin=1),xlab=ylabel,ylab=xlabel, horizontal=T,xlim=xlimits, col=colorvector,
                      scales=list(y=list(cex=labsize,rot=0, col=ifelse(row.names(spe.table) %in% c("A111AU_B", "A152WI_B", "S102LR_A", "S031BU_C", "O111BA_C", "S171MA_A", "S153TR_C", "Z042CP_C", "Z071SI_A", "Z122OU_C"), 2, 1)),x=list(cex=0.8,rot=0)),
                      key=list(col=colorvector,space="right",title=titlelegend,
                               text=list(textlegends,cex=1,col="black"),rect=F))
plot(plotreads)
#************************ END OF BARCHART ************************************************







#****************************************************************************************
#                         BARCHART3 (proportions OTU´s)
#****************************************************************************************
windows(title="Proportions of groups ", 16, 11)

spe.table <- as.table(as.matrix(speotu.agg))

xlimits <- c(0,1)
xlabel <- "proportion of taxonomic groups (based on number of OTU´s)"
ylabel <- "relative No. of OTU´s"
textlegends <- as.vector(dimnames(spe.table)[[2]])

plototus <- barchart(prop.table(spe.table,margin=1),xlab=ylabel,ylab=xlabel, horizontal=T,xlim=xlimits, col=colorvector,
                     scales=list(y=list(cex=labsize,rot=0, col=ifelse(row.names(spe.table) %in% c("A111AU_B", "A152WI_B", "S102LR_A", "S031BU_C", "O111BA_C", "S171MA_A", "S153TR_C", "Z042CP_C", "Z071SI_A", "Z122OU_C"), 2, 1)),x=list(cex=0.8,rot=0)),
                     key=list(col=colorvector,space="right",title=titlelegend,
                              text=list(textlegends,cex=1,col="black"),rect=F))
plot(plototus)
#************************ END OF BARCHART ************************************************

#*****************************************************************************************
#                       PIECHART
#*****************************************************************************************
windows(title="Average composition in sites",10,10)

spe.table <- as.table(as.matrix(sperrare.agg))

textlegends <- as.vector(dimnames(spe.table)[[2]])
maintitle <- paste(filename," (rarefaction at: ",as.character(rlevel),")")
pienums <- colSums(sperrare.agg)
textrotation <- 0
pie(pienums, labels=textlegends, main=maintitle,
    col=colorvector,clockwise=T,init.angle=90,cex=0.8,srt=textrotation)
#*********************  END OF PIECHART **************************************************