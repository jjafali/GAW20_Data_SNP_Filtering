# SNP_Filtering
SNPs filtering  and machine learning classification  of the metabolic disease syndrome using the GAW20 data 

GAW20 Data Analysis
Introduction 
This analysis applied four R scripts (described below) to filter multiple SNPs files (snps folder) and their corresponding annotation files (annotation folder) for 1105 participants (described in COVAR.csv metadata file).
Data 
1)	COVAR.csv: a metadata of multiple phenotypes for 1105 participants in the GAW20 data sets
2)	snps: a folder containing multiple genotype (SNPs) data
3)	annotation: a folder for genotype information (annotations) for the SNPs in the snps folder
 R scripts 
GAW20_Functions2022June06.R
This R script compiled the functions for filtering and classification of SNPs data.

GAW20_FilteringCode2022June06.R
This R script was applied for all the SNPs filtering including machine learning-based approaches (calling the GAW20_Functions2022June06.R script)

GAW20_ML_Classification2022June06.R
This script was applied to assess classification performance of the selected SPNs candidate sets using the Metabolic Disease Syndrome (IDF) phenotype (calling the GAW20_Functions2022June06.R script).

GAW20_ML_PerfSummary2022June06.R
This script was applied to summarize the outputs of the GAW20_ML_Classification2022June06.R scripts 

![image](https://user-images.githubusercontent.com/53820359/185928462-e9eaad20-6bd0-425c-bf31-d0f019ecd254.png)
