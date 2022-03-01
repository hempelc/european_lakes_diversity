# Function hcoplot()
# Reorder and plot dendrogram with colors for groups and legend
#
# Usage:
# hcoplot(tree = hclust.object, diss = dissimilarity.matrix, k = nb.clusters, 
#	title = paste("Reordered dendrogram from",deparse(tree$call),sep="\n"))
#
# License: GPL-2 
# Author: Francois Gillet, 23 August 2012

"hcoplot" <- function(tree, diss, k, 
	title=paste("Reordered dendrogram from", deparse(tree$call), sep="\n"))
{
  colorvec <- c("blue","red","green","brown","cyan",
                "darkmagenta","orange","forestgreen","gold","steelblue1","gray55","hotpink",
                "darkorange1","darkgreen","brown3","pink","darkolivegreen1",
                "lightcyan","aquamarine3","skyblue","lightblue","lightgrey",
                "magenta","pink2","lightyellow",
                "darkgrey","beige","white","turquoise","gray","purple","bisque")
	require(gclus)
  require(dendextend)
	gr <- cutree(tree, k=k)
  tor <- as.dendrogram(reorder.hclust(tree, diss))
 	tor <- color_labels(tor, col=ifelse(labels(tor) %in% c("A111AU_B", "A152WI_B", "S102LR_A", "S031BU_C", "O111BA_C", "S171MA_A", "S153TR_C", "Z042CP_C", "Z071SI_A", "Z122OU_C"), 2, 1)) 
 	tor <- color_branches(tor, k=k, col=colorvec)
  #labels_cex(tor) <- 2

 	plot(tor,main=title)
	mtext(paste(length(gr),"sites"), side=1, line=4.2)
	mtext(paste(k,"clusters"), side=1, line=5.2)
	so <- gr[tor$order]
	gro <- numeric(k)
	for (i in 1:k)
	{
		gro[i] <- so[1]
		if (i<k) so <- so[so!=gro[i]]
	}
	legend("topright", paste("Cluster",1:k), fill=colorvec[1:k], bty="o", inset=c(-0.1,0.2), xpd=T)
}

