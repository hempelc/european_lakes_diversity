setwd("M:/Allgemeine Daten/Christopher/Bachelor/Christopher-20Seen/")
library(data.table)
library(dplyr)
library(vegan)

spe <- fread("just_Chrysophyceae_Synurophyceae_Oikomonadaceae_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2_rrarefied.csv")
tax <- fread("just_Chrysophyceae_Synurophyceae_Oikomonadaceae_of_dataset_christopher_lakes_data_OTUs_zusammengefasst_lakesumlevel_bigger_than_100_otusumlevel_bigger_than_2_rrarefied.tax")
tax <- as.data.frame(lapply(tax,  FUN=function(x) substr(x, nchar(x)-32, nchar(x))))
spe <- select(spe, -taxonomy)
spe <- apply(spe, 2, function (x) scan(text=x, dec=",", sep=".")) #gsub schneller, musste aber schnell gehen und hat nach mehreren versuchen nicht geklappt
spe <- cbind(tax, spe)
rownames(spe) <- spe[,1]
spe <- spe[,-1]
spe <- t(spe)
spe <- as.data.frame(spe)


Chromophyton <- select(spe, contains("Chromophyton"))
Chromophyton <- colSums(t(Chromophyton))

Chromulina <- select(spe, contains("Chromulina"))
Chromulina <- colSums(t(Chromulina))

Dinobryon <- select(spe, contains("Dinobryon"))
Dinobryon <- colSums(t(Dinobryon))

Hibberdia <- select(spe, contains("Hibberdia"))
Hibberdia <- colSums(t(Hibberdia))

Mallomonas <- select(spe, contains("Mallomonas"))
Mallomonas <- colSums(t(Mallomonas))

Ochromonas <- select(spe, contains("Ochromonas"))
Ochromonas <- colSums(t(Ochromonas))

Paraphysomonas <- select(spe, contains("Paraphysomonas"))
Paraphysomonas <- colSums(t(Paraphysomonas))

Pedospumella <- select(spe, contains("Pedospumella"))
Pedospumella <- colSums(t(Pedospumella))

Poterioochromonas <- select(spe, contains("Poterioochromonas"))
Poterioochromonas <- colSums(t(Poterioochromonas))

Spumella <- select(spe, contains("Spumella"))
Spumella <- colSums(t(Spumella))

Synura <- select(spe, contains("Synura"))
Synura <- colSums(t(Synura))

Uroglena <- select(spe, contains("Uroglena"))
Uroglena <- colSums(t(Uroglena))


Mixotrophic <- data.matrix(rbind(Uroglena, Dinobryon, Chromulina, Ochromonas, Poterioochromonas))
#Mixotrophic <- decostand(Mixotrophic,"hel")
Mixotrophic_row.names_axis_labels <- row.names(Mixotrophic)
Mixotrophic_row.names_axis_labels [5] <- paste("Poterio-", "\n", "ochromonas")
Heterotrophic <- data.matrix(rbind(Pedospumella, Paraphysomonas, Spumella))
#Heterotrophic <- decostand(Heterotrophic,"hel")
Phototrophic <- data.matrix(rbind(Chromophyton, Hibberdia, Mallomonas, Synura))
#Phototrophic <- decostand(Phototrophic,"hel")


colorvec <- c("red","blue","darkgreen","orange",
              "purple","turquoise","lightgreen","brown","magenta","grey","black",
              "yellow3","lightyellow","lightblue","pink","green",
              "lightcyan","darkcyan","darkblue", "lightgrey",
              "darkmagenta","lightpink","white","darkgrey","beige",
              "yellow","gray","purple","bisque")

#########################################################################################################################################
x11(40,25,title="Trophic groups lakes sums")
par(mfrow=c(2,2))

heterotrophic_plot1 <- barplot(colSums(Heterotrophic), main="Heterotrophic", ylab="reads", col=colorvec, las=2, ylim=c(0,1.1*max(colSums(Heterotrophic))))
text(x = heterotrophic_plot1, y = colSums(Heterotrophic), label = colSums(Heterotrophic), pos = 3, cex = 0.9)

mixotrophic_plot1 <- barplot(colSums(Mixotrophic), main="Mixotrophic", ylab="reads", col=colorvec, las=2, ylim=c(0,1.1*max(colSums(Mixotrophic))))
text(x = mixotrophic_plot1, y = colSums(Mixotrophic), label = colSums(Mixotrophic), pos = 3, cex = 0.9)

phototrophic_plot1 <- barplot(colSums(Phototrophic), main="Phototrophic", ylab="reads", col=colorvec, las=2, ylim=c(0,1.3*max(colSums(Phototrophic))))
text(x = phototrophic_plot1, y = colSums(Phototrophic), label = colSums(Phototrophic), pos = 3, cex = 0.9)
############################################################################################################################
x11(40,25,title="Trophic groups species in lakes")
par(mfrow=c(2,2))

barplot(Heterotrophic, main="Heterotrophic", ylab="reads", las=2, col=colorvec, legend.text=row.names(Heterotrophic), args.legend=list(x="top"), ylim=c(0,1.1*max(colSums(Heterotrophic))))
text(x = heterotrophic_plot1, y = colSums(Heterotrophic), label = colSums(Heterotrophic), pos = 3, cex = 0.9)

barplot(Mixotrophic, main="Mixotrophic", ylab="reads", las=2, col=colorvec, legend.text=row.names(Mixotrophic), args.legend=list(x="top"), ylim=c(0,2*max(Mixotrophic)))
text(x = mixotrophic_plot1, y = colSums(Mixotrophic), label = colSums(Mixotrophic), pos = 3, cex = 0.9)

barplot(Phototrophic, main="Phototrophic", ylab="reads", las=2, col=colorvec, legend.text=row.names(Phototrophic), args.legend=list(x="top"), ylim=c(0,1.3*max(Phototrophic)))
text(x = phototrophic_plot1, y = colSums(Phototrophic), label = colSums(Phototrophic), pos = 3, cex = 0.9)


x11(40,25,title="Trophic groups species in lakes (relative)")
par(mfrow=c(2,2), mar=c(5,4,4,8)+0.1)

barplot(prop.table(Heterotrophic, margin=2), main="Heterotrophic", ylab="reads", las=2, col=colorvec)
legend("topright",legend=rev(row.names(Heterotrophic)), inset=c(-0.238,0.2), xpd=T, fill=rev(colorvec[1:3]), bty="n")

barplot(prop.table(Mixotrophic, margin=2), main="Mixotrophic", ylab="reads", las=2, col=colorvec)
legend("topright",legend=rev(row.names(Mixotrophic)), inset=c(-0.25,0.2), xpd=T, fill=rev(colorvec[1:5]), bty="n")

barplot(prop.table(Phototrophic, margin=2), main="Phototrophic", ylab="reads", las=2, col=colorvec)
legend("topright",legend=rev(row.names(Phototrophic)), inset=c(-0.2,0.2), xpd=T, fill=rev(colorvec[1:4]), bty="n")

##############################################################################################################
x11(40,25,title="Trophic groups species sums")
par(mfrow=c(2,2))

heterotrophic_plot2 <- barplot(rowSums(Heterotrophic), main="Heterotrophic", ylab="reads", col="red", las=1, ylim=c(0,1.3*max(rowSums(Heterotrophic))))
text(x = heterotrophic_plot2, y = rowSums(Heterotrophic), label = rowSums(Heterotrophic), pos = 3, cex = 0.9)

mixotrophic_plot2 <- barplot(rowSums(Mixotrophic), main="Mixotrophic", ylab="reads", col="blue", las=1, ylim=c(0,100000))
text(x = mixotrophic_plot2, y = rowSums(Mixotrophic), label = rowSums(Mixotrophic), pos = 3, cex = 0.9)

phototrophic_plot2 <- barplot(rowSums(Phototrophic), main="Phototrophic", ylab="reads", col="darkgreen", las=1, ylim=c(0,1.1*max(rowSums(Phototrophic))))
text(x = phototrophic_plot2, y = rowSums(Phototrophic), label = rowSums(Phototrophic), pos = 3, cex = 0.9)

###########################################################################################################################

x11(40,25,title="Trophic groups compared")
groups_total <- data.matrix(rbind(sum(colSums(Heterotrophic)), sum(colSums(Mixotrophic)), sum(colSums(Phototrophic))))
rownames(groups_total) <- c("Heterotrophic", "Mixotrophic", "Phototrophic")
colnames(groups_total) <- "Total"
groups_total_plot <- barplot(t(groups_total), main="Trophic groups compared", ylab="reads", las=1, col="aquamarine3")
text(x = groups_total_plot, y = groups_total, label = groups_total, pos = 3, cex = 0.9)
mtext(sum(Mixotrophic), cex=0.9, at=1.94, line=0.3)


################################################################################################################################

x11(40,25,title="Species compared")
groups_species <- rbind(Heterotrophic, Mixotrophic, Phototrophic)
groups_species_plot <- barplot(rowSums(groups_species), main="Species compared", ylab="reads", las=1, col=c("red","red", "red", "blue", "blue", "blue", "blue", "blue", "darkgreen", "darkgreen", "darkgreen", "darkgreen"), cex.names=0.75, names.arg=c(row.names(Heterotrophic), Mixotrophic_row.names_axis_labels, row.names(Phototrophic)))
#text(x = groups_species_plot, y = rowSums(groups_species), label = rowSums(groups_species), pos = 3, cex = 0.9)
#mtext(sum(Mixotrophic[2,]), cex=0.9, at=5.5, line=0.3)
legend("topleft", legend=c("Heterotrophic", "Mixotrophic", "Phototrophic"), fill=c("red", "blue", "darkgreen"), inset=c(0.06,0))


############################################################################################################
groups_sums <- rbind(colSums(Heterotrophic), colSums(Mixotrophic), colSums(Phototrophic))
row.names(groups_sums) <- c("Heterotrophic", "Mixotrophic", "Phototrophic")

x11(40, 25, title="Trophic groups in lakes (absolute)")
groups_sums_plot <- barplot(groups_sums, main="Trophic groups in lakes (absolute)", ylab="reads", las=2, col=colorvec, legend.text=row.names(groups_sums), args.legend=list(x="topright", inset=0.3), ylim=c(0,1.1*max(groups_sums)), cex.names=0.9)
text(x = groups_sums_plot, y =colSums(groups_sums), label=colSums(groups_sums), pos = 3, cex = 0.9)


x11(40, 25, title="Trophic groups in lakes (relative)")
par(mar=c(5,4,4,8)+0.1)
barplot(prop.table(groups_sums, margin=2), main="Trophic groups in lakes (relative)", ylab="reads (relative)", las=2, col=colorvec, cex.names=0.9)
legend("topright",legend=rev(row.names(groups_sums)), inset=c(-0.1,0.2), xpd=T, fill=rev(colorvec[1:3]))
