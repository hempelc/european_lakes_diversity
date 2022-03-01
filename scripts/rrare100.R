require(vegan)
rrarefy100 <- function(spe,rlev)
{ rs <- rowSums(spe)
  num <- round(max(rs)/rlev)
  if (num<=100) num <- 100 else num <- 1000
  cat("rrarefy number: ",num,"\n")  
  spe100 <- rrarefy(spe,rlev)
  for (i in 2:num)
  { spe_i <- rrarefy(spe,rlev)
    spe100 <- spe100 + spe_i   }
  spe100 <- spe100/num  
}