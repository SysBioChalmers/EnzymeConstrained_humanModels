library(pheatmap)
library(readxl)
library(RColorBrewer)

#raw data from matlab
HepG2 <- as.matrix(read_excel("~/Dropbox/Academia/PhD/Cancer/Enzyme-Constrained HMR/GEM/EnzymeConstrained-HMR-GEM/Results/KO_fluxDist_HepG2.xlsx"))

HepG2 <- data.frame(HepG2)
rownames(HepG2) <- HepG2[,1]
HepG2 <- HepG2[,-1]

HepG2_filtered <- as.matrix(HepG2)
class(HepG2_filtered) <- "numeric"

#2. replace zeros by minimum non-zero flux obtained in the WT solution
minval <- 3.8114e-13
for(i in 1:length(HepG2_filtered)){
  if(HepG2_filtered[i] < 0){
    HepG2_filtered[i] <- HepG2_filtered[i]-minval
  }
  else
    HepG2_filtered[i] <-  HepG2_filtered[i]+minval
}

#3. Apply ratio over wild-type
HepG2_filtered <- data.frame(HepG2_filtered)
HepG2_FC <- HepG2_filtered[,2:length(HepG2_filtered)]
for(i in 1:length(HepG2_FC[1,])){
  HepG2_FC[,i] <- HepG2_FC[,i]/HepG2_filtered[,1]
}
#4. Apply log and filter zero fluxes
HepG2_FC <- log(abs(HepG2_FC))*sign(HepG2_FC)
HepG2_FC <- HepG2_FC[which(rowSums(HepG2_FC) != 0),]
HepG2_FC <- HepG2_FC[,-8]
#5. heatmap representation
colr <- rev(brewer.pal(9,"RdBu"));colr[5] <- "#FFFFFF";colr <- c("#103458",colr,"#6A0D19")

pheatmap(HepG2_FC[1101:1128,],
         fontsize = 4,
         fontsize_col = 5,
         #scale = 'row',
         color = colr)

#6. Interesting ratios: 80% of the columns same direction or opposite direction
ratios <- as.matrix(read.csv("~/Dropbox/Academia/PhD/Cancer/Enzyme-Constrained HMR/GEM/EnzymeConstrained-HMR-GEM/Results/InterestingRatios_Normal.csv", sep=";",header = F))
pheatmap(HepG2_FC[ratios,],
         fontsize = 6,
         fontsize_col = 7,
         #scale = 'row',
         color = colr,
         border_color=NA)
         



