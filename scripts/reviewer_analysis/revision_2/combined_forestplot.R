# step 10

# packages
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggforestplot)

# glycosuria GWAS
## glycosuria GWAS ALSPAC
glycosuria_ALSPAC <- fread("./001_glycosuria/GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final_sig.txt", header = T)
glycosuria_ALSPAC <- subset(glycosuria_ALSPAC, P <= 5e-8)
glycosuria_ALSPAC <- subset(glycosuria_ALSPAC, P == min(glycosuria_ALSPAC$P))
glycosuria_ALSPAC$data <- "Mother glycosuria ALSPAC"
glycosuria_ALSPAC <- glycosuria_ALSPAC[, c(1:7,10,17,18,16,19)]
colnames(glycosuria_ALSPAC) <- c("SNP", "CHR", "BP", "EA", "NEA", "info", "N", "EAF", "BETA", "SE", "P", "data")

## glycosuria GWAS NFBC
glycosuria_NFBC <- fread("./001_glycosuria/GWAS/NFBC/NFBC_GWAS_final.txt", head = T, sep = "\t")
glycosuria_NFBC <- subset(glycosuria_NFBC, SNP == glycosuria_ALSPAC$SNP[[1]])
mum_NFBC <- glycosuria_NFBC[,c(1:5,11,6,7:10)]
mum_NFBC$data <- "Offspring raw glycosuria NFBC"
colnames(mum_NFBC) <- c("SNP", "CHR", "BP", "EA", "NEA", "info", "N", "EAF", "BETA", "SE", "P", "data")

offspring_NFBC <- glycosuria_NFBC[,c(1:5,11,6,7:10)]
offspring_NFBC$data <- "Offspring doubled glycosuria NFBC"
offspring_NFBC$BETA <- offspring_NFBC$BETA * 2
offspring_NFBC$SE <- offspring_NFBC$SE * 2
colnames(offspring_NFBC) <- c("SNP", "CHR", "BP", "EA", "NEA", "info", "N", "EAF", "BETA", "SE", "P", "data")

## combine
glycosuria <- rbind(glycosuria_ALSPAC, mum_NFBC, offspring_NFBC)
glycosuria$GWAS <- "glycosuria (yes/no)"
dim(glycosuria)

# test GWAS
### data
mother <- fread("./001_glycosuria/reviewer_analysis/GWAS_revision2/step5/output/GWAS_output_final.txt", header = T)
mother_SNP <- subset(mother, P == min(mother$P))
mother_SNP$data <- "Mother hair-colour raw"
mother_SNP <- mother_SNP[, c(1:7,10,17,18,16,19)]
colnames(mother_SNP) <- c("SNP", "CHR", "BP", "EA", "NEA", "info", "N", "EAF", "BETA", "SE", "P", "data")

offspring <- fread("./001_glycosuria/reviewer_analysis/GWAS_revision2/step9/output/GWAS_output_final.txt", header = T)
offspring_SNP <- subset(offspring, SNP == mother_SNP$SNP[[1]])

offspring_SNP <- offspring_SNP[, c(1:7,10,17,18,16)]
offspring_SNP$data <- "Offspring hair-colour raw"
colnames(offspring_SNP) <- c("SNP", "CHR", "BP", "EA", "NEA", "info", "N", "EAF", "BETA", "SE", "P", "data")

offspring_doubled_SNP <- offspring_SNP
offspring_doubled_SNP$BETA <- offspring_doubled_SNP$BETA * 2
offspring_doubled_SNP$SE <- offspring_doubled_SNP$SE * 2
offspring_doubled_SNP$data <- "Offspring hair-colour doubled"

## combine
hair_colour <- rbind(mother_SNP, offspring_SNP, offspring_doubled_SNP)
hair_colour$GWAS <- "blond hair (yes/no)"
dim(hair_colour)

# combine both GWAS
data <- rbind(glycosuria, hair_colour)
write.table(data, "GWAS_comparisons.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# make plot
data <- read.table("reviewer_analysis/GWAS_revision2/mother_GWAS/GWAS_comparisons.txt", header = T, sep = "\t")

levels(data$data)
data$data = factor(data$data, levels(data$data)[c(6,3,1,5,4,2)])

data$OR <- exp(data$BETA)

##### Plot
ci <- 0.95
psignif <- 0.05
pdf("reviewer_analysis/plots/GWAS_forestplot.pdf", width = 14)
forestplot(df = data, 
           name = SNP, 
           estimate = OR, 
           se = SE, 
           pvalue = P, 
           colour = data, 
           shape = GWAS, 
           logodds = FALSE, 
           psignif = psignif, 
           ci = ci, 
           xlab = "OR (95% CI) per effect allele") +
  theme(
    legend.position = "bottom",
    legend.title.align = 0,
    legend.text.align = 0,
    legend.title = element_blank(),
    axis.title.x = element_text(hjust = 0.5)
  )
dev.off()




