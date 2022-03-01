# extrapolation of exspected OTU numbers 
# input:  spe.row = site to be estimated, rlev = desired rarefactionlevel
# output: schao, sACE = exspected number of OTU´s for rlevel according to chao1 and ACE
# based on Shen T-J, Chao A, Lin C-F (2003):
# Predicting the number of new species in further taxonomic sampling.
# Ecology 84, 798-804


extraspe <- function(spe.row,rlev,ch)
{ f1        <- sum(spe.row==1)  # number of OTU´s with reads=1 (frequency f1)
  n         <- sum(spe.row)    # number of reads (too low!, therefore extrapolation)
  m         <- rlev - n        # number of extrareads (extending n, n+m=rlevel)
  sextra    <- estimateR(spe.row) # calculation of specnumber [1], chao1 [2] and ACE [4]
  fchao     <- sextra[2,] - sextra[1,] # chao1 number of extra-OTU´s estimated for rlevel
  fACE      <- sextra[4,] - sextra[1,] # ACE   number of extra-OTU´s estimated for rlevel
  schao <- sextra[1,] + fchao * (1- (1-f1/(n*fchao))** m) # formula by Shen et al. 2003
  sACE  <- sextra[1,] + fACE  * (1- (1-f1/(n*fACE)) ** m) # formula by Shen et al. 2003
  if (ch==T) {return(schao)}
      else   {return(sACE)}
}

# (partial) extrapolation of spedrare to extrapolated species richness schao or sACE
extraspedrare <- function(sp,sper,rlevel,chao)
{reads <- rowSums(sp)                            #number of reads in spe
 ilast <- nrow(sp)                               #number of sites
 for (i in 1:ilast)
   if (reads[i]< rlevel)                          #extrapolation only for low reads
      { sorig <- sp[i,]                           #row i
        specs <- specnumber(sp[i,])               #number of OTU´s in row i
        skorr <- extraspe(sorig,rlevel,chao)      #skorr = extrapolated OTU number für row i
        korrfactor <- skorr/specs                 #scalar: correction factor > 1
        sper[i,]   <- sper[i,] * korrfactor       #correct spedrare, row i
        cat(row.names(sorig),"\t")
        cat("extrapol.factor: ",korrfactor,"\t\t")  #control of correction: site and korrfactor
        cat("#OTU´s after extrapol.= ",skorr,"\n")  #print extrapolated number
      }
 return(sper)                                     #function result: corrected spedrare
}
