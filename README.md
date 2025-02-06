# PGS_LDPred2_GenR
This repository was made based on Polygenic scores and inference using LDpred2 tutorial by Florian Privé. However we made some changes based on our need for calculating Poly Genetic risk score here in Generation R study.
To run this script you need : 
1.Summary statistics of the trait. This summary statistics should contain marginal effect sizes, their standard errors, and the corresponding sample size(s).
2.An LD (linkage disequilibrium) matrix derived from individuals sharing the same genetic ancestry as those included in the GWAS.Here we used the Hapmap+ LD refrence that can be downloaded from "https://ndownloader.figshare.com/files/25503788". Remember that partitioning LD matrices into independent LD blocks can enhance robustness and significantly improve computational efficiency based on this papar https://www.sciencedirect.com/science/article/pii/S2666247722000525?via%3Dihub.
This script concists og multiple steps 
Step1 : loading and matching 
# Load the required packages
library(bigsnpr)


 Obtain HapMap3 SNPs and LD correlation matrix downloadable at https://ndownloader.figshare.com/files/25503788
 Use https://figshare.com/ndownloader/files/37802721 for HapMap3+ variants
info <- readRDS("/path/to/map.rds") # retrieved from https://figshare.com/articles/dataset/European_LD_reference/13034123?file=25503788
#########################################
 LOAD AND TRANSFORM SUMMARY STATISTICS #
#########################################

 Read in the summary statistic file
sumstats_all <- bigreadr::fread2("/path/to/summary_statistics.txt") 

 Select the needed variables
 LDpred 2 requires the headers with the following exact naming
 Modify the name of sumstats_all$columns in order to match the required info
sumstats <- data.frame(chr= sumstats_all$CHROM,
                       pos= sumstats_all$POS,
                       rsid= sumstats_all$ID,
                       a0= sumstats_all$ALT,
                       a1= sumstats_all$REF, #effect allele for which BETA/OR are calculated
                       # n_eff= sumstats_all$NEFFDIV2,
                       n_case= sumstats_all$NCAS,
                       n_control= sumstats_all$NCON,
                       beta_se= sumstats_all$SE,
                       p= sumstats_all$PVAL,
                       beta= sumstats_all$BETA
                       # OR= sumstats_all$OR
                       # INFO= sumstats_all$INFO
)


# If OR is provided instead of beta, transform the OR into log(OR)
sumstats$beta <- log(sumstats$OR)

# In case of binary trait with separate N for cases and controls
sumstats$n_eff <- 4 / (1 / sumstats$n_case + 1 / sumstats$n_control)
sumstats$n_case <- sumstats$n_control <- NULL

# Filter out HapMap3 SNPs
sumstats <- sumstats[sumstats$rsid %in% info$rsid,] 




