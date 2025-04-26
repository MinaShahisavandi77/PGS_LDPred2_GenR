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
This script consists of multiple steps
-
Step1 : loading and matching

Load the required packages
```
library(bigsnpr)
```

Obtain HapMap3 SNPs and LD correlation matrix downloadable at https://ndownloader.figshare.com/files/25503788
Use https://figshare.com/ndownloader/files/37802721 for HapMap3+ variants

```
info <- readRDS("/path/to/map.rds") # retrieved from https://figshare.com/articles/dataset/European_LD_reference/13034123?file=25503788
```

Read in the summary statistics file
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
Duplicates SNPs need to be removed! it is important to always check ! 
```
# Remove duplicates based on 'rsid' and create a clean dataset
table(duplicated(sumstats$rsid))
sumstats <- sumstats[!duplicated(sumstats$rsid), ]
```
 Filter out HapMap3 SNPs
```
sumstats <- sumstats[sumstats$rsid %in% info$rsid,]
#filter out the variants with low sample size 
max_sample_size <- max(sumstats$n_eff)
```
 Set a threshold for minimum sample size
```
threshold <- 0.7 * max_sample_size
# Filter variants that have a sample size below the threshold
filtered_sumstats <- sumstats[sumstats$n_eff >= threshold, ]
#filtered_variants<- sumstats$n_eff >= threshold
# Output the filtered variants
#write.table(filtered_variants, file="filtered_variants.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
```
 Match the variables of the summary statistics with genotype 
```
# Extract the SNP information from the genotype
map <- data.frame(chr= info$chr,
                  rsid= info$rsid,
                  pos= info$pos,
                  a1= info$a1,
                  a0= info$a0)

# Perform SNP matching
#filtered_sumstats$chr <- as.numeric(filtered_sumstats$chr) only if you get the error in the next step
#map$chr <- as.numeric(map$chr)


df_beta <- snp_match(filtered_sumstats, map)
df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)  # use rsid instead of pos as a solution to previous error
# If the error persists, try:
map_mod<- snp_modifyBuild(info_snp=map,
                          liftOver='/path/to/liftOver', # retrieved from: https://privefl.github.io/bigsnpr/reference/snp_modifyBuild.html
                          from = "hg18",
                          to = "hg19",
                          check_reverse = TRUE)
df_beta <- snp_match(sumstats, map_mod)
# Alternatively, according to the documentation, try:
df_beta <- snp_match(sumstats, map, join_by_pos = FALSE, match.min.prop=0.05)
```
---
Step 2:QC

Filter based on the differences between standard deviations of allele frequencies between UKBB genotypes and SUMSTATS
1.	Create [info_snp] = [info] restricted to [sumstats]
2.	Create a var called [sd_ldref] directly by extracting allele info from [info_snp] 
3.	Calculate [sd_ss] by extracting beta and beta_se from [df_beta] and calculating [sd_y] (for continuous traits by using a formula (the estimate should be a number around 1)

```


# IF YOU WANT TO PERFORM QC OF SUMMARY STATISTICS #


#info_snp <- tidyr::drop_na(tibble::as_tibble(info)) # most of the time needs to be silent
info_snp<- info[info$rsid %in% df_beta$rsid,]

# Better to use af of GWAS and INFO scores as well (then can use 0.7 instead of 0.5 in L35)
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))

# IF BINARY TRAIT
sd_ss <- with(df_beta, 2 / sqrt(n_eff * beta_se^2 + beta^2))

# IF CONTINUOUS TRAIT
sd_y = with(df_beta, sqrt(quantile(0.5 * (n_eff * beta_se^2 + beta^2), 0.01)))
sd_ss = with(df_beta, sd_y / sqrt(n_eff * beta_se^2 + beta^2))

# Estimate SNPs to remove
sd_ldref <- sd_ldref[1:length(sd_ss)]
is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05
```
you can check via the graph for removed SNPs
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
ggsave("QC1.png") 
```
```
#df_beta_no_qc <- df_beta  # if you want to skip the QC steps
df_beta <- df_beta[!is_bad, ]
#str(df_beta)
```
Step 3:
LDPRED correlations using blocks
```

# COMPUTE LDpred2 SCORES GENOME-WIDE #


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
Step 4: 
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


# ESTIMATE ADJUSTED BETAS WITH THE LDpred2 AUTOMATIC MODEL #


coef_shrink <- 0.95  # Reduce this up to 0.4 if you have some (large) mismatch with the LD ref

set.seed(42)  # to get the same result every time

# Can take minutes to hours
multi_auto <- snp_ldpred2_auto(corr,
                               df_beta,
                               h2_init = ldsc_h2_est,
                               vec_p_init = seq_log(1e-4, 0.2, length.out = NCORES), # 0.2 as tutorial, 0.9 as in the paper
                               ncores = 42,
                               #use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
                               allow_jump_sign = FALSE, shrink_corr = coef_shrink)
str(multi_auto, max.level = 1)                               
str(multi_auto[[1]], max.level = 1)
```
Step 6: Check Convergence and Filter
Plot Chains from multi.auto: Verify if the chains converge properly and look good.
Apply Filtering: Ensure only well-converged chains are retained for subsequent analyses.

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

ggsave("chains.png")
```
Step 7 : BETAS

 Filter the chains and check for valid 'beta_est'

 In the LDpred2 paper, it proposed an automatic way of filtering bad chains by comparing the scale of the resulting predictions.
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
write.table(adjusted_betas_reverse, "PRS_adjusted_weights.txt", sep="\t", dec =".", quote= FALSE, row.names=FALSE, col.names= FALSE)
```
To calculate PRSs for each given phenotype, you can use PLINK.

How to Download and Install PLINK:
Windows OS and macOS: Download and install it from the official PLINK website.
Linux OS: Install it directly from the terminal by typing:
```
sudo apt install plink1.9
```
Important Notes:
If you are using Windows OS or macOS, you can run PLINK from the terminal. Ensure the PLINK executable file (e.g., plink.exe or plink) is in the same working directory as your weight files.

# Calculate the polygenic risk score 

Once you have installed PLINK, you can obtain PRSs for your data using the following command (example for Linux OS):
```
	plink1.9 --bfile my_genetic_data --score my_scores.txt 3 4 6 sum --out output_filename
```
Explanation of Parameters:

- my_genetic_data: File name (without extension) of your genetic data after applying QC and imputation protocols.
- my_scores.txt: File name (with extension) of the provided adjusted weights.
- output_filename: Desired file name (without extension) for the output.
- 3 4 6: Columns representing the SNP ID, effect allele, and adjusted weights, respectively.
- sum: Option to calculate PRSs as the sum of valid per-allele scores (instead of the default mean).

Important Note:
The input file my_genetic_data must be in PLINK format (.bim, .bed, and .fam).

If your genetic data is in VCF format, you can convert it to PLINK format using the following command:
```
plink1.9 --vcf my_genetic_data.vcf.gz --keep-allele-order --make-bed --out my_genetic_data
```
For a comprehensive guide on PLINK commands and flags, visit the PLINK documentation.
 

 
 


 


 



