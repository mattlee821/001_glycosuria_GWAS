#script for producing tables

cd ./tables

R

#script for producing tables

library(data.table)
library(dplyr)
meta <- fread("./meta_analysis/meta_analysis_final_sig.txt", header = T)
GWAS_sig <- fread("./GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final_sig.txt", header = T)
GWAS <- fread("./GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final.txt", header = T)
LD_rs13337037 <- fread("./GWAS/ALSPAC/step9_LD/step2/output/data_chr16_sample_filter_subsetted.ld", header = T)
LD_rs10991823 <- fread("./GWAS/ALSPAC/step9_LD/step2/output/data_chr09_sample_filter_subsetted.ld", header = T)
NFBC <- fread("./GWAS/NFBC/NFBC_GWAS_final_parental.txt", header = T)

############ make table of SNPs in block of high LD (>= 0.8) with lead SNP
dim(LD_rs13337037)
head(LD_rs13337037)
LD_rs13337037 <- LD_rs13337037[,c(6,7)]
colnames(LD_rs13337037)[1] <- "SNP"
GWAS_r2 <- GWAS_sig
GWAS_r2 <- left_join(GWAS_sig, LD_rs13337037, by = "SNP")
dim(GWAS_r2)
GWAS_r2 <- subset(GWAS_r2, R2 >= 0.8)
dim(GWAS_r2)
head(GWAS_r2)
GWAS_r2 <- GWAS_r2[,c(1,2,3,4,5,6,7,10,13,14,15,16,17,18,19)]
write.table(GWAS_r2, "Supplementary Table 3 block of SNPs in high LD_rs13337037.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

dim(LD_rs10991823)
head(LD_rs10991823)
LD_rs10991823 <- LD_rs10991823[,c(6,7)]
colnames(LD_rs10991823)[1] <- "SNP"
GWAS_r2 <- GWAS
GWAS_r2 <- left_join(GWAS, LD_rs10991823, by = "SNP")
dim(GWAS_r2)
GWAS_r2 <- subset(GWAS_r2, R2 >= 0.8)
dim(GWAS_r2)
head(GWAS_r2)
GWAS_r2 <- GWAS_r2[,c(1,2,3,4,5,6,7,10,13,14,15,16,17,18,19)]
write.table(GWAS_r2, "Supplementary Table 4 block of SNPs in high LD_rs10991823.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

############ make table of SNPs reaching genome wide sig in ALSPAC
GWAS_sig <- GWAS_sig[,c(1,2,3,4,5,6,7,10,13,14,15,16,17,18)]
colnames(GWAS_sig) <- c("rsID", "chr", "position", "EA", "NEA", "info", "N", "EAF", "all_OR", "all_OR_lower", "all_OR_upper", "P", "B", "SE")
GWAS_sig <- GWAS_sig[,c("rsID", "chr", "position", "EA", "NEA", "EAF", "P", "B", "SE", "all_OR", "all_OR_lower", "all_OR_upper", "info", "N")]
write.table(GWAS_sig, "Supplementary Table 2 SNPs reaching genome wide significance in ALSPAC.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

############ make table of SNPs reaching genome wide sig in NFBC
NFBC_sig <- subset(NFBC, P <= 5e-8)
dim(NFBC_sig)
head(NFBC_sig)
NFBC_sig <- NFBC_sig[order(P),] 
head(NFBC_sig)
names(NFBC_sig)
NFBC_sig <- NFBC_sig[,c(1,2,3,4,5,7,8,11,12,13,14,15,9,6)]
head(NFBC_sig)
colnames(NFBC_sig) <- c("rsID", "chr", "position", "EA", "NEA", "EAF", "P", "B", "SE", "all_OR", "all_OR_lower", "all_OR_upper", "info", "N")
head(NFBC_sig)
write.table(NFBC_sig, "Supplementary Table 5 SNPs reaching genome wide significance in NFBC.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

############ make table of SNPs reaching genome wide sig in meta-analysis
head(meta)
names(meta)
meta_sig <- meta[,c(1:6,10:13,18:20)]
colnames(meta_sig) <- c("rsID", "chr", "position", "EA", "NEA", "EAF", "effect", "SE", "P", "direction", "OR", "ci_lower", "ci_upper")
meta_sig <- meta_sig[,c("rsID", "chr", "position", "EA", "NEA", "EAF", "P", "effect", "SE", "OR", "ci_lower", "ci_upper", "direction")]
write.table(meta_sig, "Supplementary Table 6 SNPs reaching genome wide significance from meta-analysis.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


############ check meta SNPs are all in GWAS and no additions
GWAS_sig <- GWAS_sig[order(P),] 
meta_sig <- meta_sig[order(P),]
dim(GWAS_sig)
dim(meta_sig)
head(GWAS_sig)
head(meta_sig)
tail(GWAS_sig)
tail(meta_sig)

# what SNP is only present in ALSPAC GWAS
GWAS_only <- anti_join(GWAS_sig, meta_sig, by = "rsID", all = T)
dim(GWAS_only)
head(GWAS_only)

# what SNP is only present in meta analysis
meta_only <- anti_join(meta_sig, GWAS_sig, by = "rsID", all = T)
dim(meta_only)
head(meta_only)

# what is the SNP overlap in ALSPAC and meta-analysis
overlap <- left_join(GWAS_sig, meta_sig, by = "rsID", all = T)
dim(overlap)
head(overlap)
tail(overlap)

write.table(overlap, "overlap_ALSPAC_meta.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
