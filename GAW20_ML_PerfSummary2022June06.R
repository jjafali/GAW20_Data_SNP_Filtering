
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
##Loading data
##=========================================================================================
##(1)Phenotype data
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

##(2)Candidate SNP features  
load(file = "CandFeatList2022June13.RData")
sapply(feat.list, length)

##(3)Eligible ML algorithms
BinMethods.list1<-list(KNN ='knn',NB = 'naive_bayes', ROCC= 'rocc',
                       LDA='lda2', GLMNET='glmnet', SVM = 'svmRadial',RF = 'rf')
##(4)Performance 
load(file="PredPerfList2022June06.RData")
sapply(pred.list, dim)
sapply(pred.list, class)

##=======================================================================================
##Performance metric function
##=======================================================================================

##Test input:dim(pdf1<-pred.list[[1]])
Perf.fn<-function(pdf1){
   names(pdf1)
  ##Confusion matrices 
 table(pdf1$Outcome)
 (methods1<-setNames(names(BinMethods.list1),names(BinMethods.list1)))
  (plist1<-lapply(methods1, function(x) 
   unlist2(confusionMatrix(data=pdf1[,x], reference=pdf1$Outcome, positive="Positive"))))
 
  metrics1<-c(Accuracy="Accuracy", BAccuracy="Balanced Accuracy", Sensitivity= "Sensitivity",
              Specificity="Specificity",PPV="Pos Pred Value", NPV="Neg Pred Value") 
  (perf.mat1<-t(sapply(plist1, function(x) sapply(metrics1, function(y) as.numeric(x[y])))))
  
  (output.df<-data.frame(Outcome=unique(pdf1$OutcomeVariable), SNPset=unique(pdf1$SNPset), 
                    SNPs=unique(pdf1$SNPs), Algorithm=methods1, perf.mat1)) 
##End function 
 
  return(output.df)
}

##=======================================================================================
##Summary plot (all the SNPs sets)
##=======================================================================================
##(i)Performance metrics data 
names(pred.list)
perf.df1<-Reduce(rbind, lapply(pred.list, function(x)  Perf.fn(pdf1=x)))

##(ii)Plot items
#Plot colors and pch
unique(perf.df1$Algorithm)
bin.pch<-c(KNN=1,NB=2, ROCC=8, LDA=17, GLMNET=18, SVM=15,RF=9)
bin.cols<-c(KNN="cyan",NB="darkgreen", ROCC="deepskyblue4", LDA="blue",
                               GLMNET="red", SVM="violet",RF="gold")
(methods1<-names(bin.pch))

##Data cleaning
names(perf.df1)
(keep.vars1<-c("SNPset","SNPs","Algorithm","Accuracy"))
plot.df<-reshape(data=perf.df1[,keep.vars1], v.names = "Accuracy",
                 timevar = "Algorithm", idvar = "SNPset", direction = "wide")
(names(plot.df)<-gsub(pattern = "Accuracy.",replacement = "", x = names(plot.df)))
plot.df<- plot.df[order( substr(plot.df$SNPset,1,2), plot.df$SNPs, decreasing = FALSE),]
plot.df$Median<-apply(plot.df[, methods1],1, median)
  
##Bar plot
pdf(file="PerfSummaryPlot2022June16.pdf", height = 12, width = 12)
par(mar=c(11,6,4,6))
bpx1<-barplot(plot.df$SNPs, ylim=c(0, 150), las=2, 
                names.arg =paste(plot.df$SNPset,"(", plot.df$SNPs, ")", sep = ""),
              density = c(rep(100, 13), rep(30, 4)))
text(c(3, 18), 0, c("FCBF", "BORUTA"), col="blue", font=2, pos=3, cex=1.5)
  
##Line plots
  par(new = TRUE)
  plot(bpx1,plot.df$Median, type="b", ylim=c(0,0.9), pch=16, lwd=2, axes=FALSE,xlab="", ylab="")
 
##Looping over methods 
lapply(methods1, function(x) lines(bpx1,plot.df[,x], type = "b", 
                                       col=bin.cols[x], pch=bin.pch[x]))
##Labels 
axis(side=4, at=seq(0,0.9, length=10), las=1, col="blue4") 
mtext("SNPset(size)", side=1, cex=2, font=1, line=8)
mtext("Number of SNPs", side=2, cex=2, font=1, line=3.5)
mtext("Perfomance Accuracy", side=4, cex=2, font=1, line=4)
mtext("FCBF=Fast correlation based filter", side=1, cex=1, font=3, line=8.5, adj=0)
mtext("BORUTA=Random forest based filter", side=1, cex=1, font=3, line=9.5, adj=0)

abline(h=seq(0,1, 0.1), col="gray10", lty=3, lwd=0.5)
abline(v=bpx1[which.max(plot.df$Median)],col="gray10", lty=3, lwd=3)

##Legend
legend("bottom", legend = names(bin.pch), col=bin.cols, pch=bin.pch, 
       pt.cex=3, cex=1.1, ncol=2, y.intersp = 1.5,bg = "gray90")
dev.off()

##---------------------------------------------------------------------------------------
##Summary plot of best SNP set
##---------------------------------------------------------------------------------------
pdf(file="BestSNPsetPerfPlot2022June16.pdf", height = 10, width = 10)
par(mar=c(8,7,4,2))
 ##Bar Plot here  
(max.snpset<-plot.df$SNPset[which.max(plot.df$Median)]) 
plot2.df<-perf.df1[perf.df1$SNPset%in%max.snpset, ]
plot2.df<-plot2.df[order(plot2.df$BAccuracy),]
    
plot(plot2.df$Accuracy, ylim=c(0.5, 1), pch=16, cex=3, type="b", lty=3, 
     xaxt="n", las=1, xlab="", ylab="")
abline(h=seq(0,1,0.1), v=seq(1.5, 7,1), lty=3, lwd=0.5, xpd=FALSE)
errbar(x=1:7, y =plot2.df$Accuracy, yplus = plot2.df$Sensitivity,
       yminus = plot2.df$Specificity, add=TRUE) 
points(1:7, plot2.df$Sensitivity, col="red", cex=2,pch=17)
points(1:7, plot2.df$Specificity, col="blue", cex=3,pch=18)

##Labels
axis(side = 1, at=1:7, labels = plot2.df$Algorithm, las=2 )
mtext("ML algorithms", side=1, line=6, cex=2)
mtext("Perfomance accuracy", side=2, line=3.5, cex=2)
legend("bottomright", c("Accuracy", "Sensitivity", "Specifity"),cex=1.5, 
       y.intersp = 1.5, col=c("black", "red", "blue"), pch=c(16,17,18), pt.cex=3)
dev.off()



