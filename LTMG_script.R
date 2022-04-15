

library(LTMGSCA)
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Rtsne)
library(RColorBrewer)
library(KRLS)
library(RSpectra)
library(SwarmSVM)
library(reshape2)


#devtools::install_github("clwan/LTMG",force=TRUE)

args = commandArgs(trailingOnly=TRUE)


File_data <- as.matrix(read.csv(args[1], sep = ',', header = TRUE, row.names = 1))
File_data <- File_data[rowSums(File_data) > 0, colSums(File_data) > 0]
Zcut_G<-log(Global_Zcut(File_data))


Gene<- args[2]
VEC<-File_data[,Gene]
VEC_LTMG<-LTMG(VEC,Zcut_G,k=5)

### Visualize the LTMG gene expression distribution
plot_gene(VEC = VEC,Data_LTMG = VEC_LTMG,Zcut = Zcut_G)




Gene_use<-rownames(File_data)[order(apply(File_data, 1, var),decreasing = T)[1:2000]]
Gene_use<-intersect(Gene_use,rownames(File_data))

File_LTMG<-LTMG_MAT(MAT = File_data,Zcut_G = Zcut_G,Gene_use = Gene_use)


write.table(t(as.matrix(File_LTMG$State)), args[3] , sep = "\t", col.names = NA, row.names = TRUE)




