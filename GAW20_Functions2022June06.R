##Date: 2022 June 06
##H3ABIONET (Machinelearning work package)
#This scrpt contains R functions for the main analyses of the GAW20 data 
##=========================================================================================
##Data filtering function
##=========================================================================================
##Test input:geno.file1<-"geno13_15.csv"; pheno1<-"IDF_msdx"; stat.pval1<-0.05

SNPsFilt.fn<-function(geno.file1, pheno1, stat.pval1){
  
  t00<-proc.time()
  print(rep(paste(pheno1, geno.file1, sep = "||"),3))
  ##Loading data
  dim(gdf<-read_csv(paste(snps.path1,"/",  geno.file1, sep="")))
  length(intersect(gdf$GAWSUBJ, pdat$GAWSUBJ))
  dim(mdf1<-merge(pdat, gdf, by="GAWSUBJ")); rm(gdf)
  length(snps1<-setdiff(names(mdf1),names(pdat)))
  
  ##Monomorphic SNPs
  table(monomorphics1<-apply(mdf1[,snps1], 2, function(x) length(table(x))<2))
  length(snps2<-snps1[!monomorphics1])
  
  ##Missing values (samples)
  pmiss.samps<-apply(mdf1[,snps2],1, function(x) 100*(sum(is.na(x))/length(snps2)))
  table(pmiss.samps>5)
  dim(mdf1<-mdf1[pmiss.samps<5,])
  
  ##Missing values (SNPs)
  (pmiss.snps<-apply(mdf1[,snps2],2, function(x) 100*(sum(is.na(x))/ncol(mdf1))))
  table(pmiss.samps>5)
  table(pmiss.snps>5)
  length(pmiss.snps); length(snps2)
  
  ##Minor allele frequency (MAF)>5%
  (MAFs <- apply(mdf1[,snps2], 2, function(x) 1- mean(x, na.rm=TRUE)/2))
  table(MAFs>0.05)
  
  ##Hardy-Weinberg filter (p-value>1%)
  head(hwd.df<-t(apply(mdf1[,snps2],2,  function(x)
    c(AA=sum(x==0, na.rm=TRUE), BB=sum(x==1,na.rm=TRUE), AB=sum(x==2,na.rm=TRUE)))))
  dim(hwd.df)
  
  HWB.Chisq.Pvals <- HWChisqStats(X =hwd.df, x.linked=FALSE,pvalues=TRUE)
  table(HWB.Chisq.Pvals>0.01)
  
  ##Outcome variable 
  table(yy<-mdf1[, pheno1])
  table(is.na(yy))
  dim(mdf1<-mdf1[!is.na(yy),])
  
  ##Statistical filtering
  if(length(unique(yy))>5){
    print("Continous outcome (Kruskal Wallis test)!!")
    
    (StatFilter<-sapply(snps2,function(y) 
      kruskal.test(x = yy,g = factor(mdf1[,y]))$p.value))
    
  } else {
    print("Categorical outcome (Chi-square test)!!")
    (StatFilter<-sapply(snps2,function(y1) chisq.test(x = yy,y=factor(mdf1[,y1]))$p.value))

  }
  #hist(StatFilter, main="Statistics filter", col="deepskyblue4", xlab="Association P-values")
  table(StatFilter<0.05)
  
  ##Output files
  (chromo<-gsub(pattern = "geno", replacement = "Chromosome", 
                x = unlist(strsplit(x = geno.file1, "_"))[1]))
  filt.df<-data.frame(Chromosome=chromo, SNPs=snps2, Missingness=pmiss.snps,
                      MAF=MAFs,HWB.pVAL=HWB.Chisq.Pvals, StatFilter=StatFilter,
                      stringsAsFactors =FALSE )
  table(filt.df$Keep<-(filt.df$Missingness < 5 & filt.df$MAF>0.05 &
                         filt.df$HWB.pVAL>0.01 & filt.df$StatFilter < stat.pval1))
  
  (sel.snps<-filt.df$SNPs[filt.df$Keep])
  
  if(length(sel.snps)>1) { 
    snps.df<-data.frame(mdf1[, c("GAWSUBJ", pheno1)], 
                        apply(mdf1[,sel.snps], 2, function(y) 
                          factor(x = y, levels = 0:2, labels = c("AA", "AB", "BB"))))
    
    output.list<-list(snps.df= snps.df, filt.df=filt.df)
  } else {
    output.list<-NA  
  }
  sapply(output.list, dim)
  
  ##End of Function
  print(proc.time()-t00)
  
  return(output.list)
}


##=======================================================================================
##Mode function
##=======================================================================================
Mode.fn <- function(yy1) {
  uniqx <- unique(na.omit(yy1))
  uniqx[which.max(tabulate(match(yy1, uniqx)))]
}


##=======================================================================================
##ML Classification function
##=======================================================================================
# ##Test inputs
# snp.set1<-feat.list$FCBF0.016; Outcome<-"IDF_msdx"
# snpset.name1<-"FCBF"; K1<-5; R1<-3; K2<-5; R2<-3; seed1<-1235
# ml.methods1<- list(KNN ='knn',NB = 'naive_bayes',ROCC= 'rocc')

ClassfyML.fn<-function(Outcome, snp.set1, snpset.name1, 
                       K1, R1, K2,R2, ml.methods1, seed1) {
  t00<-proc.time()
  print(rep(paste(Outcome,  "||", snpset.name1, "=", length(snp.set1),sep = ""),5))
  
  ##(1) Model data.frame
  dim(mdf<-snp.df[, c(Outcome, snp.set1 )])
  dim(mod.mat<- model.matrix(~.-1, data = mdf[,-1]))
  dim(mod.df<-data.frame(Outcome=snp.df[row.names(mod.mat), Outcome],mod.mat))
  ##Remove AA variable
  (keep.vars2<-setdiff(names(mod.df), grep("AA", names(mod.df), value = TRUE)))
  mod.df<-mod.df[, keep.vars2]
  
  ##(2)Splitting data into cross-validation folds
  set.seed(seed1)
  (DataFolds<-createMultiFolds(y = mod.df$Outcome,k=K1, times = R1))
  sapply(DataFolds,length)
  
  ##[3]Hyper-parameter tuning function
  TrainByEachFold.fn<-function(train.fold1) {
    ##Training and testing sets
    dim(train.df<-mod.df[train.fold1,])
    dim(test.df<-mod.df[-train.fold1,])
    
    ##Parameter tuning & model fitting
    train_control <- trainControl(method="repeatedcv", number=K2, repeats=R2,search = "grid")
    fit.list<-lapply(ml.methods1, function(x)
      train(Outcome ~., data = train.df, method = x, trControl = train_control,
            preProcess = c("center","scale"), tuneLength = 10))
    
    ##Model prediction
    (pred.mat<- sapply(fit.list, function(x) predict(object = x,newdata=test.df[,-1])))
    
    ##Output data
    (out.df<-data.frame(SampleID=row.names(test.df), Outcome=test.df$Outcome, pred.mat))
    
    return(out.df)
  }
  
  ##[4] Looping over the data folds (application of the TrainByEachFold.fn function)
  pred.df.list<-lapply(DataFolds, function(y)  TrainByEachFold.fn(train.fold1 = y))
  lapply(pred.df.list, dim)
  dim(pred.df<-Reduce(rbind.data.frame, pred.df.list))
  
  ##Consolidation of replicates(mode/mean value)
  pred.list2<-split(x = pred.df, f = pred.df$SampleID)
  lapply(pred.list2, dim)
  dim(output.df<-data.frame(Reduce(rbind, 
                                   lapply(pred.list2, function(D) 
                                     apply(D,2, function(k) Mode.fn(yy1 = k))))))
  
  ##Output
  output.df<-output.df%>% mutate(SNPset=snpset.name1, SNPs=length(snp.set1))
  output.df$OutcomeVariable<-Outcome
  print(proc.time()-t00)
  ##END of FUNCTION
  
  return(output.df)
}