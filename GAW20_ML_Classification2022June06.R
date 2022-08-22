
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
              "WGCNA","gplots","ggplot2", "glmnet","metaMA", "metap", "epiR", "dplyr")
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
         ATP_msdx=factor(x= atp_msdx, levels=0:1, labels=c("Negative", "Positive")),
         IDF_msdx=factor(x= idf_msdx, levels=0:1, labels=c("Negative", "Positive")))
names(pdat)
##=========================================================================================
##SNPs data 
##=========================================================================================
##Loading save data frames
load(file="FilterDF.RData"); dim(filt.df)
load(file="SNPsDF.RData"); dim(snp.df)

####Imputation of missing data using the Mode.fn function
table(nmiss.snp<-apply(snp.df,2, function(x) sum(is.na(x))))
Mode.fn(yy1 = snp.df$rs1541318XXX)

for (x in names(snp.df)){
  y1<-snp.df[,x]
  y1[is.na(y1)]<-Mode.fn(yy1 = y1)
  snp.df[,x]<-y1
}

class(snp.df); table(is.na(snp.df))

##=======================================================================================
##Candidate SNP features 
##=======================================================================================
##(1) Candidate feature list
load(file = "CandFeatList2022June06.RData")
sapply(feat.list, length)

##=======================================================================================
##ML classification using the ClassfyML.fn function 
##=======================================================================================
##Eligible methods 
BinMethods.list1<-list(KNN ='knn',NB = 'naive_bayes', ROCC= 'rocc',
                       LDA='lda2', GLMNET='glmnet', SVM = 'svmRadial',RF = 'rf')
##SNPS lists
sapply(feat.list, length); length(feat.list)

# ##Predictions
# pred.list<-mapply(function(x,y)
#   ClassfyML.fn(Outcome="IDF_msdx", snp.set1=x, snpset.name1=y,
#                           K1=5, R1=3, K2=3,R2=3, ml.methods1=BinMethods.list1, seed1=12345),
#   x=feat.list, y=names(feat.list), SIMPLIFY = FALSE)

##Saving the performance list
#save(pred.list, file="PredPerfList2022June06.RData")
load(file="PredPerfList2022June06.RData")
sapply(pred.list, dim)
sapply(pred.list, class)


