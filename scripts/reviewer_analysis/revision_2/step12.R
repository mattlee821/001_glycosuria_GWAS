# Reviewer comment: Did the authors check the association between maternal glycosuria in ALSPAC and the fetal genotype at the confirmed chr16 SNP? This is recommended for 2 reasons: (i) showing that the fetal genotype association in ALSPAC is approx. half that of the maternal genotype would add confidence that the NFBC1986 analysis (using fetal genotype as proxy) is valid; (ii) since the glycosuria is measured in pregnancy, there is a chance that the fetal genotype would have a role in influencing it. Therefore the association should be checked, and in addition, an analysis of maternal genotype, conditional on fetal genotype should be carried out. # ====
### glycosuria correlation ====

# data
## extract rs13337037 dosage
module add apps/qctool-2.0
qctool -g ./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/dosage_bgen/data_chr16.bgen \
-s ./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/data.sample \
-incl-range 31384553-31522650  \
-assume-chromosome 16 \
-ofiletype dosage \
-og ./001_glycosuria/reviewer_analysis/revision_2/GWAS_revision2/step12/all_chr16_sig.dosage


# open R and load packages
R
library(dplyr)

sig_snps <- read.table("./001_glycosuria/GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final_sig.txt", header = T, sep = "\t")
head(sig_snps)[1:10]
sig_snps <- subset(sig_snps, chr == 16)
SNPs <- sig_snps[[1]]
SNPs

# transpose dosage files and ensure each aln_qlet has correct dosage 
dosage <- read.table("./001_glycosuria/reviewer_analysis/revision_2/GWAS_revision2/step12/all_chr16_sig.dosage", header = T, fill = T, sep = "", quote = "")
dosage2 <- dosage
dim(dosage2)
head(dosage2)[1:7]
dosage2 <- dosage2[dosage2$rsid %in% SNPs, ]
dim(dosage2)
head(dosage2)[1:7]

dosage3 <- t(dosage2)
dim(dosage3)
head(dosage3, 7)
dosage3 <- dosage3[-1:-2,]
dosage3 <- dosage3[-2:-4,]
rownames(dosage3) <- gsub("X([0-9]+)", "\\1", rownames(dosage3))

write.table(dosage3, "./001_glycosuria/reviewer_analysis/revision_2/GWAS_revision2/step12/all_chr16_final.dosage", 
            row.names = T, col.names = T, quote = FALSE, sep = "\t", )

# edit row 1 header manually
dosage <- read.table("./001_glycosuria/reviewer_analysis/revision_2/GWAS_revision2/step12/all_chr16_final-1.dosage", header = T, sep = "\t")
head(dosage)

## offspring genetic data 
genetic_data <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/data.sample", skip=2)
colnames(genetic_data) <- c("ID_1","ID_2","missing","father","mother","sex","plink_pheno")
colnames(genetic_data)[1] <- "aln_qlet"

## split genetic_data aln and qlet column into two
genetic_data$qlet <- gsub("[0-9]", "", genetic_data$aln_qlet)
genetic_data$aln <- gsub("[a-zA-Z]", "", genetic_data$aln_qlet)
genetic_data$qlet <- as.factor(genetic_data$qlet)
table(genetic_data$qlet)
head(genetic_data)

## remove all mums
genetic_data <- subset(genetic_data, qlet != "M")
table(genetic_data$qlet)
head(genetic_data)
dim(genetic_data)

## save offspring sample file
genetic_data <- genetic_data[,c(1:7)]
write.table(genetic_data, "./001_glycosuria/reviewer_analysis/revision_2/GWAS_revision2/step12/offspring.sample", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# merge offspring sample file with dosage file
head(genetic_data)
str(dosage)
head(dosage)
library(data.table)

colnames(dosage)[1] <- "join_col"
colnames(genetic_data)[5] <- "join_col"
dosage <- as.data.frame(dosage)
data1 <- left_join(genetic_data, dosage, by = "join_col")
head(data1)[1:10]

data1$ID_1 <- data1$aln_qlet
data1$ID_1 <- gsub("([0-9]+)[A-Z a-z]", "\\1", data1$ID_1)
data1$ID_1 <- paste(data1$ID_1, "M", sep = "")
head(data1)
offspring_dosage <- data1[,c(1,8:ncol(data1))]
head(offspring_dosage)[1:6]
dim(offspring_dosage)

## load mothers pheno_file
pheno_file <- read.table("./001_glycosuria/other_analysis/regressions/pheno_file.txt", header = T)
head(pheno_file)
dim(pheno_file)

## join mothers_phenofile with offspring_dosage
pheno_file <- inner_join(pheno_file, offspring_dosage, by = "ID_1")
dim(pheno_file)
head(pheno_file)
str(pheno_file)

# phenotype
phenotype <- "Self-reported glycosuria (yes/no)"
group <- "offspring genotype"
note <- "raw"

# logistic regression of motehrs glycosuria and offspring genotype
data <- pheno_file
df_list <- list()

for(i in 25:ncol(data)) {
  SNP <- colnames(data)[i]
  #print(metab)
  eval(parse(text = paste("temporary <- glm(",SNP,"~ e113, data = data)", sep =""))) 

  #print(temporary)
  
  rsid <- as.character(temporary[["terms"]][[2]])
  df <- temporary[["df.residual"]]
  n <- temporary[["df.residual"]] + temporary[["rank"]]
  b <- summary(temporary)[["coefficients"]][2,1]
  se <- summary(temporary)[["coefficients"]][2,2]
  p <- summary(temporary)[["coefficients"]][2,4]
  
  temp_df <- data.frame(SNP = rsid, n = n, b = b, se = se, p = p)
  j=i-24
  df_list[[j]] <- temp_df
}  

table <- do.call(rbind, df_list)
table
write.table(table, "./001_glycosuria/reviewer_analysis/revision_2/GWAS_revision2/step12/glm_logit_glycosuria_child-dosage.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
#####

ls()
mothers <- read.table("./001_glycosuria/GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final_sig.txt", header = T, sep = "\t")
mothers <- subset(mothers, chr == 16)
offspring <- read.table("./001_glycosuria/reviewer_analysis/revision_2/GWAS_revision2/step12/glm_logit_glycosuria_child-dosage.txt", header = T, sep = "\t")
data <- left_join(mothers, offspring, by = "SNP")

names(mothers)
names(offspring)
mean(data$b / data$BETA)
png("./001_glycosuria/reviewer_analysis/revision_2/GWAS_revision2/step12/histogram_glycosuria_ratio.png")
hist(data$b / data$BETA, main = "Glycosuria: chromsome 16, 54 SNPs")
dev.off()



