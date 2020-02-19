# step 10

# packages
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggforestplot)

### data
mother <- fread("./001_glycosuria/reviewer_analysis/GWAS_revision2/step5/output/GWAS_output_final.txt", header = T)
mother_sig <- subset(mother, P <= 5e-8)
mother_SNP <- subset(mother, P == min(mother$P))
mother_SNP$data <- "Mother raw"
offspring <- fread("./001_glycosuria/reviewer_analysis/GWAS_revision2/step9/output/GWAS_output_final.txt", header = T)
offspring_sig <- subset(offspring, P <= 5e-8)
offspring_SNP <- subset(offspring, SNP == mother_SNP$SNP[[1]])
offspring_SNP$data <- "Offspring raw"

# make significant SNP list from mums and extract offspring betas for those SNPs
data_sig <- mother_sig
data_sig$group <- "mum"
sig_SNPs <- data.frame(data_sig[,1])
mother_sig_offspring <- left_join(sig_SNPs, offspring, by = "SNP")
mother_sig_offspring$group <- "offspring"
data_sig <- rbind(data_sig, mother_sig_offspring)
write.table(data_sig, "sig_snps_mum-offspring.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_sig <- read.table("reviewer_analysis/GWAS_revision2/mother_GWAS/sig_snps_mum-offspring.txt", header = T, sep = "\t")
mum <- subset(data_sig, group =="mum")
offspring <- subset(data_sig, group == "offspring")
offspring <- offspring[,c("SNP","BETA")]
data <- left_join(mum, offspring, by = "SNP")

# double offspring beya
data$offspring_doubled <- data$BETA.y * 2

# divide mum beta by offspring doubled beta
data$mum_beta_divided_by_offspring_doubled_beta <- data$BETA.x / data$offspring_doubled

# correlation of mum new beta by mum EAF
cor.test(data$mum_beta_divided_by_offspring_doubled_beta, data$EAF)             
pdf("mum-beta_offspring-beta_mum-eaf.pdf")
ggplot(data = data, 
       aes(x = mum_beta_divided_by_offspring_doubled_beta, y = EAF)) +
  geom_point(aes(x = mum_beta_divided_by_offspring_doubled_beta, y = EAF)) +
  ylab("Mum EAF")
dev.off()



# test difference between betas (mum, offpsirng, offspring doubled)
data <- rbind(mother_SNP, offspring_SNP)
data <- rbind(data, offspring_SNP)
data <- data[,c(1,2,3,4,5,7,8,10,16,17,18,19)]
data[3]$data <- "Offspring doubled"
data[3]$BETA <- data[3]$BETA*2
data[3]$SE <- data[3]$SE*2
data$OR <- exp(data$BETA)
data$ci_lower <- exp(data$BETA - 1.96 * data$SE) 
data$ci_upper <- exp(data$BETA + 1.96 * data$SE) 

##### Plot
ci <- 0.95
psignif <- 0.05
pdf("GWAS_beta_forestplot.pdf")
forestplot(df = data, 
           name = SNP, 
           estimate = BETA, 
           se = SE, 
           pvalue = P, 
           colour = data, 
           shape = NULL, 
           logodds = FALSE, 
           psignif = psignif, 
           ci = ci, 
           xlab = "OR (95% CI) per G allele of blond hair") +
  theme(
    legend.position = "bottom",
    legend.title.align = 0,
    legend.text.align = 0,
    legend.title = element_blank(),
    axis.title.x = element_text(hjust = 0.5)
  )
 dev.off()