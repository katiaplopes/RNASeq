# Katia Lopes
# August 02, 2019
# Basic plots for RNASeq data

# The data used here was downloaded from Embl Expression Atlas: https://www.ebi.ac.uk/gxa/home 
# We are using TPM, but you should consider another normalizations and correction by effects. 
# For demonstration purpose only! 

#This command clean all variables. BE CAREFULL!!! 
rm(list = setdiff(ls(), lsf.str()))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("RColorBrewer")) install.packages("RColorBrewer"); library("RColorBrewer")
if(!require("gplots")) install.packages("gplots"); library("gplots")
if(!require("Hmisc")) install.packages("Hmisc"); library("Hmisc")
if(!require("amap")) BiocManager::install("amap"); library("amap")
if(!require("ggfortify")) BiocManager::install("ggfortify"); library("ggfortify")
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("factoextra")) install.packages("factoextra"); library("factoextra")

work_plots = "/Users/katia/Desktop/kt/rnaseq/"

# Input data 
tpm_data = as.matrix(read.table("/Users/katia/Desktop/kt/rnaseq/fantom5_2016_tpm.txt", header = TRUE, row.names = 1, check.names = FALSE))
dim(tpm_data)

#####################################
# Trees to check outliers samples
#####################################
tpm_data_t = t(tpm_data)

#Euclidean distance with the transposed matrix
sampleTree = hclust(dist(tpm_data_t), method = "average")
#png(paste0(work_plots, "tree_euclidean.png"), width = 20, height = 10, res = 300, units = "in")
plot(sampleTree, main = "Sample clustering to detecting outliers", 
     sub = "Function: hclust, Method: average from Euclidean distance", xlab = "samples", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)
#dev.off()

#Spearman correlation
#png(paste0(work_plots, "tree_Spearman.png"), width = 20, height = 10, res = 300, units = "in")
plot(hcluster(tpm_data_t, method= "spearman", link="average" ), 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, main= "Hierarchical Clustering",
     sub = "Function: hcluster, Method: Spearman", xlab = "samples")
#dev.off()

#####################################
# PCA
#####################################
res.pca = prcomp(tpm_data_t)
#fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F     # Avoid text overlapping
)

#####################################
# MDS
#####################################
d <- dist(tpm_data_t)
mds <- cmdscale(d)
plot(mds[,1],mds[,2],pch=21,xlab="First dimension",ylab="Second dimension")
title("MDS - Euclidean Distance")

#####################################
# Boxplots
#####################################
boxplot(log2(tpm_data[, 1:56] +1), col=rainbow(56), main="TPM", ylab="Log2(TPM+1)",las=2,cex.axis=0.8)

#####################################
# Heatmaps
#####################################
#Spearman
cormatrix = rcorr(as.matrix(tpm_data), type='spearman')
corrdata = as.matrix(cormatrix$r)
heatmap.2(corrdata, main="TPM - Spearman",trace="none", col = c(sort(brewer.pal(9,"Blues")),brewer.pal(9,"Reds")), margins = c(8,8))

#Pearson
cormatrix = rcorr(as.matrix(tpm_data), type='pearson')
corrdata = as.matrix(cormatrix$r)
heatmap.2(corrdata, main="TPM - Pearson", trace="none", col = c(sort(brewer.pal(9,"Blues")),brewer.pal(9,"Reds")), margins = c(8,8))

#####################################
# Density plot
#####################################
ltpm = log2(tpm_data)
nsamples <- ncol(tpm_data)
samplenames <- colnames(tpm_data)

colfunc <- colorRampPalette(c("#4DBBD5FF", "#3C5488FF")) # You can specify other colors 
col = alpha(colfunc(nsamples), alpha = 0.1) # To put all samples in light blue

plot(density(ltpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="TPM ", xlab="Log2(tpm)")
for (i in 2:nsamples){
  den <- density(ltpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("right", samplenames, text.col=col, bty="n") # It's not a good idea when you have 56 samples! 


