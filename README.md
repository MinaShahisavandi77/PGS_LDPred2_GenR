 # PGS_LDPred2_GenR
 ---
This repository was made based on Polygenic scores and inference using LDpred2 tutorial by Florian Privé.
- Florian Privé, Julyan Arbel, Bjarni J Vilhjálmsson, LDpred2: better, faster, stronger, Bioinformatics, Volume 36, Issue 22-23, 1 December 2020, Pages 5424–5431, <https://doi.org/10.1093/bioinformatics/btaa1029>
---
However, we made some modifications based on our needs for calculating the Polygenic Risk Score (PRS) in the Generation R study.
To run this script you need : 
- Summary statistics of the trait. This summary statistics should contain marginal effect sizes, their standard errors, and the corresponding sample size(s).
- An LD (linkage disequilibrium) matrix derived from individuals sharing the same genetic ancestry as those included in the GWAS.Here we used the [HapMap3+](https://ndownloader.figshare.com/files/25503788) 
LD refrence.
- Remember that partitioning LD matrices into independent LD blocks can enhance robustness and significantly improve computational efficiency based on [this papar](https://www.sciencedirect.com/science/article/pii/S2666247722000525?via%3Dihub).
---
This script concists of multiple steps
-
Step1 : loading and matching

Load the required packages
```
library(bigsnpr)
```
```
Obtain HapMap3 SNPs and LD correlation matrix downloadable at https://ndownloader.figshare.com/files/25503788
Use https://figshare.com/ndownloader/files/37802721 for HapMap3+ variants
```
```
info <- readRDS("/path/to/map.rds") # retrieved from https://figshare.com/articles/dataset/European_LD_reference/13034123?file=25503788
```

Read in the summary statistic file
```
sumstats_all <- bigreadr::fread2("/path/to/summary_statistics.txt") 

#Select the needed variables
#LDpred2 requires the headers with the following exact naming
#Modify the name of sumstats_all$columns in order to match the required info
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

#If OR is provided instead of beta, transform the OR into log(OR)
sumstats$beta <- log(sumstats$OR)

#In case of binary trait with separate N for cases and controls
sumstats$n_eff <- 4 / (1 / sumstats$n_case + 1 / sumstats$n_control)
sumstats$n_case <- sumstats$n_control <- NULL

#Filter out HapMap3 SNPs
sumstats <- sumstats[sumstats$rsid %in% info$rsid,] 



```
duplicates need to be remove! it is important to always check ! 
```
 # Remove duplicates based on 'rsid' and create a clean dataset
table(duplicated(sumstats$rsid))
sumstats <- sumstats[!duplicated(sumstats$rsid), ]
# Filter out HapMap3 SNPs
sumstats <- sumstats[sumstats$rsid %in% info$rsid,]
#filter out the variants with low sample size 
max_sample_size <- max(sumstats$n_eff)

# Set a threshold for minimum sample size
threshold <- 0.7 * max_sample_size

# Filter variants that have a sample size below the threshold
filtered_sumstats <- sumstats[sumstats$n_eff >= threshold, ]

# Output the filtered variants
#write.table(filtered_variants, file="filtered_variants.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
##########################################################
# MATCH VARIANTS BETWEEN GENOTYPE AND SUMMARY STATISTICS #
##########################################################

# Extract the SNP information from the genotype
map <- data.frame(chr= info$chr,
                  rsid= info$rsid,
                  pos= info$pos,
                  a1= info$a1,
                  a0= info$a0)

# Perform SNP matching
filtered_sumstats$chr <- as.character(filtered_sumstats$chr)
map$chr <- as.character(map$chr)

df_beta <- snp_match(filtered_sumstats, map)
```
---
STEP 2:QC

1.Filter based on the differences between standard deviations of allele frequencies between UKBB genotypes and SUMSTATS
1.	Create [info_snp] = [info] restricted to [sumstats]
2.	Create a var called [sd_ldref] directly by extracting allele info from [info_snp] 
3.	Calculate [sd_ss] by extracting beta and beta_se from [df_beta] and calculating [sd_y] (for continuous traits by using a formula (the estimate should be a number around 1)

```

###################################################
# IF YOU WANT TO PERFORM QC OF SUMMARY STATISTICS #
###################################################

info_snp <- tidyr::drop_na(tibble::as_tibble(info))
info_snp<- info[info$rsid %in% df_beta$rsid,]

# Better to use af of GWAS and INFO scores as well (then can use 0.7 instead of 0.5 in L35)
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))

# IF BINARY TRAIT
sd_ss <- with(df_beta, 2 / sqrt(n_eff * beta_se^2 + beta^2))

# IF CONTINUOUS TRAIT
#sd_y = with(df_beta, sqrt(quantile(0.5 * (n_eff * beta_se^2 + beta^2), 0.01)))
#sd_ss = with(df_beta, sd_y / sqrt(n_eff * beta_se^2 + beta^2))

# Estimate SNPs to remove
sd_ldref <- sd_ldref[1:length(sd_ss)]
is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05
```
check the graph for removed SNPs
```
library(ggplot2)
qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "SD from allele frequencies of the LD reference",
       y = "SD from summary statistics",
       color = "Removed?")
ggsave("SCZ1.png") 
```
```
df_beta_no_qc <- df_beta
df_beta <- df_beta[!is_bad, ]
#str(df_beta)
```
STEP 3:
LDPRED correlations using blocks
```
######################################
# COMPUTE LDpred2 SCORES GENOME-WIDE #
######################################

# Create a correlation matrix

NCORES <- nb_cores()
tmp <- tempfile(tmpdir = "temp/")

for (chr in 1:22) {
  
  cat(chr, ".. ", sep = "")
  
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map$chr == chr))
  
  corr_chr <- readRDS(paste0("ldref/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3] # "LD_chr" folder includes LD blocks retrieved from: https://figshare.com/articles/dataset/European_LD_reference_with_blocks_/19213299/1?file=34133082
  
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

```
STEP4: 
Heritability estimate

```
# Estimate h2 from LD Score regression

# Perform matching
matched_indices <- match(df_beta$rsid, info$rsid)

# Filter out NAs from the matched indices
valid_indices <- !is.na(matched_indices)

# Update df_beta and ld to contain only valid (non-NA) SNPs
df_beta <- df_beta[valid_indices, ]
ld <- info[matched_indices[valid_indices], ]

ldsc <- with(df_beta, snp_ldsc(ld$ld, length(ld$ld), chi2 = (beta / beta_se)^2,
                               sample_size = n_eff, blocks = NULL))
ldsc_h2_est <- ldsc[["h2"]]

ldsc_h2_est

 ```
Step 5 
LDPRED auto
 ```

############################################################
# ESTIMATE ADJUSTED BETAS WITH THE LDpred2 AUTOMATIC MODEL #
############################################################

coef_shrink <- 0.4  # Reduce this up to 0.4 if you have some (large) mismatch with the LD ref

set.seed(42)  # to get the same result every time

# Can take minutes to hours
multi_auto <- snp_ldpred2_auto(corr,
                               df_beta,
                               h2_init = ldsc_h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES), # 0.2 as tutorial, 0.9 as in the paper
                               ncores = 42,
                               #use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
                               allow_jump_sign = FALSE, shrink_corr = coef_shrink)
str(multi_auto, max.level = 1)                               
str(multi_auto[[1]], max.level = 1)
```
Step 6: CHECK CONVERGENCE CHAINS + FILTER 
Plot chains from multi.auto (check if they look good)
```
# Verify whether the chains “converged” by looking at the path of the chains
library(ggplot2)
auto <- multi_auto[[1]]  # first chain
plot_grid(
  qplot(y = auto$path_p_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv")

ggsave("SCZ2.png")
```
STEP 7 : BETAS

 Filter the chains and check for valid 'beta_est'

 In the LDpred2 paper, we proposed an automatic way of filtering bad chains by comparing the scale of the resulting predictions.
 Here we recommend an equivalent and simpler alternative:
```
range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

# To get the final effects you should only use chains that pass this filtering
beta_auto_means <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
```
Step 8 
Calculate PGS using genotype data

```
# Save adjusted betas
adjusted_betas <- data.frame(CHR= df_beta$chr,
                             POS= df_beta$pos,
                             RSID = df_beta$rsid,
                             A0= df_beta$a0,
                             A1= df_beta$a1,
                             BETA_ADJ = beta_auto_means)
adjusted_betas_reverse <- data.frame(CHR= df_beta$chr,
                                     POS= df_beta$pos,
                                     RSID = df_beta$rsid,
                                     A1= df_beta$a1,
                                     A2= df_beta$a0,
                                     BETA_ADJ = beta_auto_means)

# For PLINK format, we use "adjusted beta reverse" because PLINK assumes A1 as effect allele and A2 as other allele
write.table(adjusted_betas_reverse, "SCZ_adjusted_weights.txt", sep="\t", dec =".", quote= FALSE, row.names=FALSE, col.names= FALSE)
```
To calculate PRSs for each given phenotype, you can use PLINK. To download PLINK:
	- Windows OS and MacOS: you can download and install it from https://www.cog-genomics.org/plink/
	- Linux OS: you can directly install it from the terminal by typing "sudo apt install plink1.9"

Please note that if you are using Windows OS and Mac OS, you can run PLINK from the terminals while having the executable file (e.g., "plink.exe" or "plink") in the same working directory of the weights files.



# CALCULATE POLYGENIC RISK SCORES #


Once you have installed PLINK, you can obtain PRSs for your data using the following command (example for Linux OS):
```
	plink1.9 --bfile my_genetic_data --score my_scores.txt 3 4 6 sum --out output_filename
```
Where:
	- "my_genetic_data" = file name (without extension) of your genetic data after your QC and imputation protocols.
	- "my_scores.txt" = file name (with extension) of the provided adjusted weights.
	- "output_filename" = file name (without extension) of the desired output files.
	- "3 4 6" = columns corresponding to SNPID, Effect allele, and Adjusted weigths, respectively.
	- "sum" = option to calculate PRSs as the sum of valid per-allele scores (instead of the mean, which is calculated as default).

Please note that the input file "my_genetic_data" must be in PLINK format (i.e., you should have the .bim, .bed, and .fam versions of the file).
If you have a VCF file, you can convert it in PLINK format with the following command:
	plink1.9 --vcf my_genetic_data.vcf.gz --keep-allele-order --make-bed --out my_genetic_data

For detailed information about PLINK commands and flags, please see: https://www.cog-genomics.org/plink/1.9/index



 

 
 


 


 



