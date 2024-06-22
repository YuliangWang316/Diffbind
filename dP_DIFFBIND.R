library(DiffBind)
tamoxifen <- dba(sampleSheet = "E:/CJJ202404166BCD19DPDNK27CUTTAG/Bam/Diffbind_CJJ202404166bcd19dpdnl27cuttag.csv")
tamoxifen <- dba.blacklist(tamoxifen, blacklist=DBA_BLACKLIST_MM10,cores = 20,greylist = FALSE)
tamoxifen <- dba.count(tamoxifen,bParallel = TRUE)
tamoxifenpeak<-cbind(as.data.frame(tamoxifen$peaks[1]),as.data.frame(tamoxifen$peaks[2]),as.data.frame(tamoxifen$peaks[3]),as.data.frame(tamoxifen$peaks[4]))
tamoxifenpeak_new<-tamoxifenpeak[,c(6,14,22,30)]
colnames(tamoxifenpeak_new)<-c("104DP","107DP","106DP","110DP")
library(tidyverse)
library(DESeq2)
#import data
setwd("e:/")
mycounts_lz<-tamoxifenpeak_new
condition_lz<-factor(c(rep("KO",2),rep("WT",2)),levels = c("WT","KO"))
colData_lz<-data.frame(row.names = colnames(mycounts_lz),condition_lz)

dds_lz <- DESeqDataSetFromMatrix(mycounts_lz, colData_lz, design= ~ condition_lz)
dds_lz<-estimateSizeFactors(dds_lz)
sizeFactors(dds_lz)<-c(1.0641347,1.0127698,1.19111486,0.8828803)
dds_lz <- DESeq(dds_lz)

res_lz= results(dds_lz)
res_lz = res_lz[order(res_lz$pvalue),]
head(res_lz)
summary(res_lz)
res_lz<-as.data.frame(res_lz)
res_lz<-na.omit(res_lz)
tamoxifenpeak_new_new<-tamoxifenpeak[rownames(res_lz),]
res_lz<-cbind(tamoxifenpeak_new_new[,c(1:4)],res_lz)
write.csv(res_lz,file="All_results_DP.csv")
