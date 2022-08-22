
setwd("~/Documents/H3ABIONET_GAW20_SNPs/Analysis")
rm(list=ls())

##=========================================================================================
#Loading packages and functions
##=========================================================================================
##(1) Packages 
lod.pckgs<-c( "HardyWeinberg", "data.table",   "pROC","foreign", "reshape2", "lattice", "limma",
              "readxl", "readr", "multtest","rocc","biomaRt","mice", "Hmisc","randomForest",
              "annotate","AnnotationDbi", "hgu219.db", "hgu133plus2.db",'illuminaHumanv2.db',
              "hgu133a.db", "hgu133b.db","GEOquery", "impute","corrplot",  "sva", "caret", "Rtsne",
              "WGCNA","gplots","ggplot2", "glmnet","metaMA", "metap", "epiR", "dplyr", "arsenal")
lod.pckgs
for(x in 1:length(lod.pckgs))
  library(package = lod.pckgs[x], character.only = TRUE)

##(2) Functions
source('~/Documents/H3ABIONET_GAW20_SNPs/Analysis/GAW20_Functions2022June06.R')

##=========================================================================================
##Phenotype data
##=========================================================================================
pdat<- read_csv("COVAR.csv")
names(pdat)
names(pdat)[names(pdat)%in%c( "_smoking","_atp_msdx","_idf_msdx")]<-c("Smoking","atp_msdx","idf_msdx")
names(pdat)


pdat<- pdat %>%
      mutate(hdl.pre=rowMeans(cbind(hdl1,hdl2), na.rm = TRUE),
             hdl.post=rowMeans(cbind(hdl3,hdl4), na.rm = TRUE),
             hdl.change=hdl.post-hdl.pre,
             Center=factor(center, levels = 0:1, labels=c("Minnesota","Utah")),
             ATP_msdx=factor(x= atp_msdx, levels=0:1, labels=c("Negative", "Positive")),
             IDF_msdx=factor(x= idf_msdx, levels=0:1, labels=c("Negative", "Positive")))
names(pdat)
hist(pdat$hdl.change, main="HDL change", xlab="HDL change", col="deepskyblue4")


pheno.tab1<-tableby( ~ age + Center + Smoking +  hdl.change +
                       ATP_msdx +  IDF_msdx, data=pdat)

pLabVar<-c( age="Age in years",Center="Study site",Smoking="Smoking status",
            hdl.change="HDL change (post-pre)",
            ATP_msdx="Metabolic disease syndrom (ATP-III)",
            IDF_msdx="Metabolic disease syndrom (IDF)")
(pheno.tab2<-summary(pheno.tab1, text=TRUE, labelTranslations=pLabVar))
write.csv(pheno.tab2, file="PhenoTypeDescriptibeTable2022June13.csv")
##=======================================================================================
##Data filtering using the SNPsFilt.fn function
##=======================================================================================
##Eligible file
snps.path1<-"~/Documents/H3ABIONET_GAW20_SNPs/Analysis/snps"
(geno.files<-list.files(path = snps.path1))
(chroms1<-sapply(geno.files, function(y) gsub(pattern = "geno", replacement = "Chromo",
                                    gsub(pattern = ".csv", replacement = "", x =  y))))
(SNPs.files1<-setNames(geno.files,chroms1))
length(SNPs.files1)
names(pdat)
# 
# ##Expected time (minutes)
# (length(SNPs.files1)*8.5)/60


##=======================================================================================
##Statistically significant SNPs
##=======================================================================================
##Stats p-value function

StatPval.fn<-function(pheno111){
  print(rep(pheno111,5))
  df.list11<-lapply(SNPs.files1, function(x)
    SNPsFilt.fn(geno.file1 = x, pheno1 = pheno111, stat.pval1 = 0.01))
  length(df.list12<-df.list11[!is.na(df.list11)])
  filt.df11<-Reduce(rbind.data.frame,lapply(df.list12, function(x) x$filt.df))
  stats.pval11<-filt.df11$StatFilter
  rm(df.list11,  filt.df11)
  return(stats.pval11) 
}

##Application of the StatPval.fn function
# age.pvals<-StatPval.fn(pheno111="age")
# smoke.pvals<-StatPval.fn(pheno111="Smoking")
# loc.pvals<-StatPval.fn(pheno111="Center")
# hdl.pvals<-StatPval.fn(pheno111="hdl.change")
# atp.pvals<-StatPval.fn(pheno111="ATP_msdx")
# idf.pvals<-StatPval.fn(pheno111="IDF_msdx")

# pval.list<-list(Age=age.pvals, Smoking=smoke.pvals,
#                     Location=loc.pvals,HDL_change=hdl.pvals,
#                     ATP_msdx=atp.pvals, IDF_msdx= idf.pvals)
# save(pval.list, file = "AssociationPvalues2022June16.RData")
load(file = "AssociationPvalues2022June16.RData")
sapply(pval.list, function(x) table(x<0.05))
sapply(pval.list, function(x) table(x<0.01))

##======================================================================================= 
##International diabetes  federation (IDF)
##=======================================================================================
# df.list1<-lapply(SNPs.files1, function(x)
#   SNPsFilt.fn(geno.file1 = x, pheno1 = "IDF_msdx", stat.pval1 = 0.01))
# 
# ##Data organization
# sapply(df.list1, names); table(is.na(df.list1))
# length(df.list2<-df.list1[!is.na(df.list1)]) #Remove missing values
# 
# ##SNPs data
# snp.df.list<-lapply(df.list2, function(x) x$snps.df)
# sapply(snp.df.list, names)
# table(sapply(snp.df.list, ncol))
# snp.df<-Reduce(f =  function(x1,y1) merge(x1,y1, by=c("GAWSUBJ", "IDF_msdx")),x =snp.df.list)
# dim(snp.df)
# save(snp.df, file="SNPsDF2022June17.RData")
# load(file="SNPsDF2022June17.RData"); dim(snp.df)
# # 
# ##Filter data frame
# filt.df.list<-lapply(df.list2, function(x) x$filt.df)
# filt.df<-Reduce(rbind.data.frame,filt.df.list)
# save(filt.df, file="FilterDF2022June17.RData")
load(file="FilterDF2022June17.RData"); dim(filt.df)
##=========================================================================================
##Filtering summary table
##=========================================================================================
names(filt.df)
table(filt.df$Keep)
sapply(pval.list, function(x) table(x<0.05))
sapply(pval.list, function(x) table(x<0.01))


(filt.ob1<-c(Total=nrow(filt.df), Missing05=sum(filt.df$Missingness<5),
            MAF=sum(filt.df$MAF>0.05), 
            HWB01=sum(filt.df$HWB.pVAL>0.01),
            Pv05=sapply(pval.list, function(x) sum(x<0.05)),
            Pv01=sapply(pval.list, function(x) sum(x<0.01)),
            FinalFilt=sum(filt.df$Keep)))
            
(filt.tab1<-stack(filt.ob1)[, c(2,1)])
names(filt.tab1)<-c("Filter", "mPass")
filt.tab1$PassPercent=round(100*filt.tab1$mPass/nrow(filt.df),2)
filt.tab1$mFail=nrow(filt.df)-filt.tab1$mPass
filt.tab1$FailPercent=round(100-filt.tab1$PassPercent, 2)
filt.tab1$Pass<-paste(filt.tab1$mPass, "(",filt.tab1$PassPercent, "%", ")", sep="")
filt.tab1$Fail<-paste(filt.tab1$mFail, "(",filt.tab1$FailPercent, "%", ")", sep="")
write.csv(filt.tab1, file="FilterSummaryDF2022June17.csv")
##=======================================================================================
##Manhattan plot
##=======================================================================================

dim(filt.df2<-filt.df[filt.df$Keep,])
filt.df2$FDR<-p.adjust(p=filt.df2$StatFilter,method = "fdr")
hist(filt.df2$FDR, main="Association P-values", xlab = "FDR", col="deepskyblue4")

table(filt.df2$Chromosome)
table(filt.df2$Chrom<-as.numeric(gsub("Chromosome", "", filt.df2$Chromosome)))
filt.df2<-filt.df2[order(filt.df2$Chrom),]
(chroms22<-setNames(1:22, paste("Chromo", 1:22,sep="")))

##The plot
pdf(file = "ManhattanPlot2022June08.pdf", height = 10, width = 15)
par(mar=c(8,6,4,2))
(xlab1<-tapply(1:nrow(filt.df2), paste("Chromo", filt.df2$Chrom, sep=""), median))
(xlab2<-tapply(1:nrow(filt.df2), paste("Chromo", filt.df2$Chrom, sep=""), max))

plot(-log10(filt.df2$FDR), col=ifelse(filt.df2$Chrom%in%seq(1,23,2),"gold", "blue"),
     las=2,pch=ifelse(filt.df2$Chromosome%in%seq(1,23,2),16, 17), xaxt="n",
     xlab="", ylab="", main="Metabolic disease syndrome (IDF)", cex.main=2)
abline(h=seq(2.1,4,0.2) , v=xlab2, lwd=1, lty=2, col="gray")
axis(1,at=xlab1,labels=names(xlab1), cex=1.5, las=2)
mtext(text="Chromosomes",side = 1, line=6,cex=2)
mtext(text="-log10(FDR)",side = 2, line=3.8,cex=2)
par(mar=c(5,5,4,2))

dev.off()
##=======================================================================================
##Missing data imputation
##=======================================================================================
table(nmiss.snp<-apply(snp.df,2, function(x) sum(is.na(x))))

##Mode function 
Mode.fn<-function(yy1) {
tab1<-table(yy1)
mod1<-names(tab1[which.max(tab1)])
return(mod1)
}
Mode.fn(yy1 = filt.df$Chromosome)

##Imputation of the snp.df data frame using the Mode.fn function
for (x in names(snp.df)){
  y1<-snp.df[,x]
  y1[is.na(y1)]<-Mode.fn(yy1 = y1)
  snp.df[,x]<-y1
}

class(snp.df); table(is.na(snp.df)) 
##=======================================================================================
##ML feature selection
##=======================================================================================
#----------------------------------------------------------------------------------------
##(1)FCBF
#----------------------------------------------------------------------------------------
# library("FCBF")
# table(snp.df$IDF_msdx)
# length(sel.snps2<-names(snp.df)[-c(1:2)])
# su_plot(x = t(snp.df[, sel.snps2]), y = snp.df$IDF_msdx)
# #
# # ##Thresholds
# (fcbf.thresholds0<-seq(0.008, 0.02, by=0.001))
# (fcbf.thresholds<-setNames(fcbf.thresholds0, paste("FCBF", fcbf.thresholds0, sep="")))
# fcbf.list<- lapply(fcbf.thresholds, function(x)
#    fcbf(x =snp.df[, sel.snps2], y = snp.df$IDF_msdx, samples_in_rows=TRUE, thresh = x))
#  sapply(fcbf.list, nrow); length(fcbf.list)
# 
# duplicated(fcbf.list)

##----------------------------------------------------------------------------------------
##(2)Boruta algorithms 
# ##----------------------------------------------------------------------------------------
# library("Boruta")
# boruta.res01<-Boruta(x=snp.df[, sel.snps2], y = snp.df$IDF_msdx, pValue = 0.01)
# boruta.res05<-Boruta(x=snp.df[, sel.snps2], y = snp.df$IDF_msdx, pValue = 0.05)
# 
# boruta.snps.list0<-list(P05a=boruta.res05$finalDecision[!boruta.res05$finalDecision%in%"Rejected"],
#                        P05b=boruta.res05$finalDecision[boruta.res05$finalDecision%in%"Confirmed"],
#                        P01a=boruta.res01$finalDecision[!boruta.res01$finalDecision%in%"Rejected"],
#                        P01b=boruta.res01$finalDecision[boruta.res01$finalDecision%in%"Confirmed"])
# 
# boruta.snps.list<-lapply(boruta.snps.list0, names)
# sapply(boruta.snps.list,length)

##----------------------------------------------------------------------------------------
##(3)Organizing the final list
##----------------------------------------------------------------------------------------
# feat.list<-c(lapply(fcbf.list, row.names), boruta.snps.list)
# sapply(feat.list, length)
# save(feat.list, file = "CandFeatList2022June13.RData")
load(file = "CandFeatList2022June13.RData")
sapply(feat.list, length)
