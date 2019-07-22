cd ./other_analysis/genetic_overlap/

module add languages/R-3.5-ATLAS-gcc-7.1.0
R

library(data.table)

####### load ALSPAC GWAS significant hits
gwas_sig <- read.table("./GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final_sig.txt", header = T)
dim(gwas_sig)
colnames(gwas_sig)
head(gwas_sig)
gwas_sig <- gwas_sig[,c(1,2,3,4,5,10,17,18,16)]
colnames(gwas_sig)[1] <- "snp"
gwas_sig$SNPID <- paste0(gwas_sig$chr, ":", gwas_sig$position)
head(gwas_sig)

####### MAGIC
magic_hba1c <- fread("./other_analysis/genetic_overlap/MAGIC/HbA1c_METAL_European.txt", header = T, sep = "\t")
dim(magic_hba1c)
colnames(magic_hba1c)
head(magic_hba1c)
magic_hba1c <- magic_hba1c[,c(1,4,5,6,7,8,9)]
colnames(magic_hba1c) <- c("snp", "magic_hba1c_EA", "magic_hba1c_OA", "magic_hba1c_eaf-hapmap-CEU", "magic_hba1c_B", "magic_hba1c_SE", "magic_hba1c_P")

magic_fg <- fread("./other_analysis/genetic_overlap/MAGIC/MAGIC_FastingGlucose.txt", header = T, sep = "\t")
dim(magic_fg)
colnames(magic_fg)
head(magic_fg)
colnames(magic_fg) <- c("snp", "magic_fg_EA", "magic_fg_OA", "magic_fg_MAF", "magic_fg_B", "magic_fg_SE", "magic_fg_P")

magic_fi <- fread("./other_analysis/genetic_overlap/MAGIC/MAGIC_ln_FastingInsulin.txt", header = T, sep = "\t")
dim(magic_fi)
colnames(magic_fi)
head(magic_fi)
colnames(magic_fi) <- c("snp", "magic_fi_EA", "magic_fi_OA", "magic_fi_MAF", "magic_fi_B", "magic_fi_SE", "magic_fi_P")

####### DIAGRAM
diagram <- fread("./other_analysis/genetic_overlap/DIAGRAM/METAANALYSIS_DIAGRAM_SE1.txt", header = T)
dim(diagram)
colnames(diagram)
diagram <- diagram[,c(1,2,3,4,5,6)]
colnames(diagram) <- c("SNPID", "diagram_EA", "diagram_OA", "diagram_B", "diagram_SE", "diagram_P")
head(diagram)

####### GIANT
giant <- fread("./other_analysis/genetic_overlap/GIANT/SNP_gwas_mc_merge_nogc.tbl.uniq", header = T)
dim(giant)
colnames(giant)
giant <- giant[,c(1,2,3,4,5,6,7)]
colnames(giant) <- c("snp", "giant_EA", "giant_OA", "giant_Freq1.Hapmap", "giant_B", "giant_SE", "giant_P")
head(giant)

####### estimated glomerular filtration rate
eGFR <- fread("./other_analysis/genetic_overlap/eGFR/CKDGen_1000Genomes_DiscoveryMeta_eGFRcrea_overall.csv")
dim(eGFR)
colnames(eGFR)
eGFR <- eGFR[,c(1,2,3,4,5,6,7)]
colnames(eGFR) <- c("snp", "eGFR_EA", "eGFR_OA", "eGFR_Freq1", "eGFR_B", "eGFR_SE", "eGFR_P")
head(eGFR)


####### merge data
library(dplyr)

# merge MAGIC with ALSPAC
data <- left_join(gwas_sig, magic_hba1c, by = "snp")
dim(data)
head(data)

data <- left_join(data, magic_fg, by = "snp")
dim(data)
head(data)

data <- left_join(data, magic_fi, by = "snp")
dim(data)
head(data)

# merge DIAGRAM with ALSPAC
data <- left_join(data, diagram, by = "SNPID")
dim(data)
head(data)

# merge GIANT with ALSPAC
data <- left_join(data, giant, by = "snp")
dim(data)
head(data)

# merge eGFR with ALSPAC
data <- left_join(data, eGFR, by = "snp")
dim(data)
head(data)

####### save merged data frame
data <- data[,c(1:9,11:45)]
write.table(data, "genetic_overlap_ALSPAC_sig_SNPs.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
###################


