# PGS_LDPred2_GenR
This repository was made based on Polygenic scores and inference using LDpred2 tutorial by Florian Privé. However we made some changes based on our need for calculating Poly Genetic risk score here in Generation R study.
To run this script you need : 
1.Summary statistics of the trait. This summary statistics should contain marginal effect sizes, their standard errors, and the corresponding sample size(s).
2.An LD (linkage disequilibrium) matrix derived from individuals sharing the same genetic ancestry as those included in the GWAS.Here we used the Hapmap+ LD refrence that can be downloaded from "https://ndownloader.figshare.com/files/25503788". Remember that partitioning LD matrices into independent LD blocks can enhance robustness and significantly improve computational efficiency based on this papar https://www.sciencedirect.com/science/article/pii/S2666247722000525?via%3Dihub.
