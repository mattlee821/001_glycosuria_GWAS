# Reviewer comment: Did the authors check the association between maternal glycosuria in ALSPAC and the fetal genotype at the confirmed chr16 SNP? This is recommended for 2 reasons: (i) showing that the fetal genotype association in ALSPAC is approx. half that of the maternal genotype would add confidence that the NFBC1986 analysis (using fetal genotype as proxy) is valid; (ii) since the glycosuria is measured in pregnancy, there is a chance that the fetal genotype would have a role in influencing it. Therefore the association should be checked, and in addition, an analysis of maternal genotype, conditional on fetal genotype should be carried out. # ====
### glycosuria correlation ====

# data
## extract rs13337037 dosage
module add apps/qctool-2.0
qctool -g ./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/dosage_bgen/data_chr16.bgen \
-s ./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/data.sample \
-incl-range 31478710-31478711  \
-assume-chromosome 16 \
-ofiletype dosage \
-og ./001_glycosuria/reviewer_analysis/revision_1/step1/all_rs13337037.dosage

# open R and load packages
R
library(dplyr)

# transpose dosage files and ensure each aln_qlet has correct dosage 
dosage <- read.table("./001_glycosuria/reviewer_analysis/revision_1/step1/all_rs13337037.dosage", header = F, fill = T, sep = "", quote = "")
head(dosage)
dim(dosage)
dosage <- t(dosage)
head(dosage, 10)
tail(dosage, 10)  
dim(dosage)
colnames(dosage) <- c("aln_qlet", "rs13337037_dosage")

write.table(dosage, "./001_glycosuria/reviewer_analysis/revision_1/step1/all_rs13337037_final.dosage", 
            row.names = FALSE, col.names = T, quote = FALSE, sep = "\t", )

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
write.table(genetic_data, "./001_glycosuria/reviewer_analysis/revision_1/step1/offspring.sample", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# merge offspring sample file with dosage file
head(genetic_data)
head(dosage)
colnames(dosage)[1] <- "join_col"
colnames(genetic_data)[5] <- "join_col"
dosage <- as.data.frame(dosage)
dosage$rs13337037_dosage <- as.character(dosage$rs13337037_dosage)
data1 <- left_join(genetic_data, dosage, by = "join_col")
offspring_dosage <- data1[,c(5,8)]
head(offspring_dosage,100)
dim(offspring_dosage)

## load mothers pheno_file
pheno_file <- read.table("./001_glycosuria/other_analysis/regressions/pheno_file.txt", header = T)
head(pheno_file)
dim(pheno_file)

## join mothers_phenofile with offspring_dosage
colnames(offspring_dosage)[1] <- "ID_1"
colnames(offspring_dosage)[2] <- "rs13337037_dosage_offspring"
pheno_file <- inner_join(pheno_file, offspring_dosage, by = "ID_1")
dim(pheno_file)
head(pheno_file)
str(pheno_file)
pheno_file$rs13337037_dosage_offspring <- as.integer(pheno_file$rs13337037_dosage_offspring)



# phenotype
phenotype <- "Self-reported glycosuria (yes/no)"
group <- "offspring genotype"
note <- "raw"

# n
table(pheno_file$e113)
n <- 3878
n_cases <- 922 

# regression
glm <- glm(pheno_file$e113 ~ pheno_file[,"rs13337037_dosage_offspring"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

b <- glm_coef[2,1]
se <- glm_coef[2,2]

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

# P
P <- glm_coef[2,4]

#OR
OR <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- data.frame("Trait" = phenotype, "Group" = group, "Note" = note, "n" = n, "n_cases" = n_cases, "b" = b, "se" = se, "OR" = OR, "CI lower" = ci_lower, "CI upper" = ci_upper, "P" = P)
dim(table)
head(table)
table$Trait <- as.character(table$Trait)
table$Note <- as.character(table$Note)
table$Group <- as.character(table$Group)

# phenotype
phenotype <- "Self-reported glycosuria (yes/no)"
group <- "offspring genotype"
note <- "doubled"

# regression
b <- (glm_coef[2,1] *2)
se <- (glm_coef[2,2] *2)

# ci
ci_lower <- (glm_coef[2,1]*2) - 1.96 * (glm_coef[2,2]*2) # 0.3228759
ci_upper <- (glm_coef[2,1]*2) + 1.96 * (glm_coef[2,2]*2) # 0.6751995

#OR
OR <- exp(glm_coef[2,1]*2)
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- (rbind(table,list(phenotype, group, note, n, n_cases, b, se, OR, ci_lower, ci_upper, P)))
dim(table)
head(table)

### glycosuria correlation conditional ====
# phenotype
phenotype <- "Self-reported glycosuria (yes/no)"
group <- "offspring genotype"
note <- "raw conditional on maternal genotype"

# n
table(pheno_file$e113, pheno_file$rs13337037_dosage_offspring)
n <- 3878
n_cases <- 922  

# regression
glm <- glm(data = pheno_file, formula = e113 ~ rs13337037_dosage_offspring + rs13337037_dosage, family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)
b <- glm_coef[2,1]
se <- glm_coef[2,2]

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

# P
P <- glm_coef[2,4]

#OR
OR <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- (rbind(table,list(phenotype, group, note, n, n_cases, b, se, OR, ci_lower, ci_upper, P)))
dim(table)
head(table)

# phenotype
phenotype <- "Self-reported glycosuria (yes/no)"
group <- "offspring genotype"
note <- "doubled conditional on maternal genotype"

# regression
b <- (glm_coef[2,1] *2)
se <- (glm_coef[2,2] *2)

# ci
ci_lower <- (glm_coef[2,1]*2) - 1.96 * (glm_coef[2,2]*2) # 0.3228759
ci_upper <- (glm_coef[2,1]*2) + 1.96 * (glm_coef[2,2]*2) # 0.6751995

#OR
OR <- exp(glm_coef[2,1]*2)
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- (rbind(table,list(phenotype, group, note, n, n_cases, b, se, OR, ci_lower, ci_upper, P)))
dim(table)
head(table)

### glycosuria correlation conditional ====
# phenotype
phenotype <- "Self-reported glycosuria (yes/no)"
group <- "maternal genotype"
note <- "conditional on offspring genotype"

# n
table(pheno_file$e113, pheno_file$rs13337037_dosage_offspring)
n <- 3878
n_cases <- 922  

# regression
glm <- glm(data = pheno_file, formula = e113 ~ rs13337037_dosage + rs13337037_dosage_offspring, family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)
b <- glm_coef[2,1]
se <- glm_coef[2,2]

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

# P
P <- glm_coef[2,4]

#OR
OR <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- (rbind(table,list(phenotype, group, note, n, n_cases, b, se, OR, ci_lower, ci_upper, P)))
dim(table)
head(table)



### glycosuria correlation mums ====
# phenotype
phenotype <- "Self-reported glycosuria (yes/no)"
group <- "maternal genotype"
note <- ""

# n
table(pheno_file$e113)
n <- 3878
n_cases <- 922  

# regression
glm <- glm(data = pheno_file, formula = e113 ~ rs13337037_dosage, family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)
b <- glm_coef[2,1]
se <- glm_coef[2,2]

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

# P
P <- glm_coef[2,4]

#OR
OR <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- (rbind(table,list(phenotype, group, note, n, n_cases, b, se, OR, ci_lower, ci_upper, P)))
dim(table)
head(table)


# phenotype
phenotype <- "Glycosuria determined by reagent strip (yes/no)"
group <- "offspring genotype"
note <- "raw"

# n
table(pheno_file$stix)
n <- 4911
n_cases <- 165 

# regression
glm <- glm(pheno_file$stix ~ pheno_file[,"rs13337037_dosage_offspring"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

b <- glm_coef[2,1]
se <- glm_coef[2,2]

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

# P
P <- glm_coef[2,4]

#OR
OR <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- (rbind(table,list(phenotype, group, note, n, n_cases, b, se, OR, ci_lower, ci_upper, P)))
dim(table)
head(table)

# phenotype
phenotype <- "Glycosuria determined by reagent strip (yes/no)"
group <- "offspring genotype"
note <- "doubled"

# regression
b <- (glm_coef[2,1] *2)
se <- (glm_coef[2,2] *2)

# ci
ci_lower <- (glm_coef[2,1]*2) - 1.96 * (glm_coef[2,2]*2) # 0.3228759
ci_upper <- (glm_coef[2,1]*2) + 1.96 * (glm_coef[2,2]*2) # 0.6751995

#OR
OR <- exp(glm_coef[2,1]*2)
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- (rbind(table,list(phenotype, group, note, n, n_cases, b, se, OR, ci_lower, ci_upper, P)))
dim(table)
head(table)


### glycosuria correlation conditional ====
# phenotype
phenotype <- "Glycosuria determined by reagent strip (yes/no)"
group <- "maternal genotype"
note <- "conditional on offspring genotype"

# n
table(pheno_file$stix, pheno_file$rs13337037_dosage_offspring)
n <- 4911
n_cases <- 165  

# regression
glm <- glm(data = pheno_file, formula = stix ~ rs13337037_dosage + rs13337037_dosage_offspring, family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)
b <- glm_coef[2,1]
se <- glm_coef[2,2]

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

# P
P <- glm_coef[2,4]

#OR
OR <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- (rbind(table,list(phenotype, group, note, n, n_cases, b, se, OR, ci_lower, ci_upper, P)))
dim(table)
head(table)



### glycosuria correlation mums ====
# phenotype
phenotype <- "Glycosuria determined by reagent strip (yes/no)"
group <- "maternal genotype"
note <- ""

# n
table(pheno_file$stix)
n <- 4911
n_cases <- 165  

# regression
glm <- glm(data = pheno_file, formula = stix ~ rs13337037_dosage, family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)
b <- glm_coef[2,1]
se <- glm_coef[2,2]

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

# P
P <- glm_coef[2,4]

#OR
OR <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- (rbind(table,list(phenotype, group, note, n, n_cases, b, se, OR, ci_lower, ci_upper, P)))
dim(table)
table

# write finished table
write.table(table, "./001_glycosuria/reviewer_analysis/offspring_dosage_glycosuria_regression.txt", 
            row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")



# Reviewer comment: Another type of analysis that could be useful and informative would be to take confirmed variants for fasting glucose, type 2 diabetes etc and draw scatter plots of their glucose or T2D effects vs. glycosuria effects. # ====
library(data.table)
library(dplyr)

## fasting glucose
magic_fg <- fread("./001_glycosuria/other_analysis/genetic_overlap/MAGIC/MAGIC_FastingGlucose.txt", header = T, sep = "\t")
head(magic_fg)
dim(magic_fg)
colnames(magic_fg)[1] <- "SNP"
colnames(magic_fg)[7] <- "P"
magic_fg <- subset(magic_fg, P <= 5e-8)
head(magic_fg)
dim(magic_fg)
magic_fg <- magic_fg[,c(1,7,5,6)]
magic_fg$trait_GWAS <- "magic_fg"
head(magic_fg)
colnames(magic_fg) <- c("SNP", "GWAS_P", "GWAS_BETA", "GWAS_SE", "trait_GWAS")
write.table(magic_fg, "./001_glycosuria/reviewer_analysis/magic_fg.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## HbA1c
magic_hba1c <- fread("./001_glycosuria/other_analysis/genetic_overlap/MAGIC/HbA1c_METAL_European.txt", header = T, sep = "\t")
head(magic_hba1c)
dim(magic_hba1c)
colnames(magic_hba1c)[1] <- "SNP"
colnames(magic_hba1c)[9] <- "P"
magic_hba1c <- subset(magic_hba1c, P <= 5e-8)
head(magic_hba1c)
dim(magic_hba1c)
magic_hba1c <- magic_hba1c[,c(1,9,7,8)]
magic_hba1c$trait_GWAS <- "magic_hba1c"
head(magic_hba1c)
colnames(magic_hba1c) <- c("SNP", "GWAS_P", "GWAS_BETA", "GWAS_SE", "trait_GWAS")
write.table(magic_hba1c, "./001_glycosuria/reviewer_analysis/magic_hba1c.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## type 2 diabetes
diagram_T2D <- fread("./001_glycosuria/other_analysis/genetic_overlap/DIAGRAM/METAANALYSIS_DIAGRAM_SE1.txt", header = T)
head(diagram_T2D)
dim(diagram_T2D)
colnames(diagram_T2D)[1] <- "SNPID"
colnames(diagram_T2D)[6] <- "P"
diagram_T2D <- subset(diagram_T2D, P <= 5e-8)
head(diagram_T2D)
dim(diagram_T2D)
diagram_T2D <- diagram_T2D[,c(1,6,4,5)]
diagram_T2D$trait_GWAS <- "diagram_T2D"
head(diagram_T2D)
colnames(diagram_T2D) <- c("SNPID", "GWAS_P", "GWAS_BETA", "GWAS_SE", "trait_GWAS")
write.table(diagram_T2D, "./001_glycosuria/reviewer_analysis/diagram_T2D.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## BMI
giant_BMI <- fread("./001_glycosuria/other_analysis/genetic_overlap/GIANT/SNP_gwas_mc_merge_nogc.tbl.uniq", header = T)
head(giant_BMI)
dim(giant_BMI)
colnames(giant_BMI)[1] <- "SNP"
colnames(giant_BMI)[7] <- "P"
giant_BMI <- subset(giant_BMI, P <= 5e-8)
head(giant_BMI)
dim(giant_BMI)
giant_BMI <- giant_BMI[,c(1,7,5,6)]
giant_BMI$trait_GWAS <- "giant_BMI"
head(giant_BMI)
colnames(giant_BMI) <- c("SNP", "GWAS_P", "GWAS_BETA", "GWAS_SE", "trait_GWAS")
write.table(giant_BMI, "./001_glycosuria/reviewer_analysis/giant_BMI.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## estimated glomerular filtration rate 
eGFR <- fread("./001_glycosuria/other_analysis/genetic_overlap/eGFR/CKDGen_1000Genomes_DiscoveryMeta_eGFRcrea_overall.csv")
head(eGFR)
dim(eGFR)
colnames(eGFR)[1] <- "SNPID"
colnames(eGFR)[7] <- "P"
eGFR <- subset(eGFR, P <= 5e-8)
head(eGFR)
dim(eGFR)
eGFR <- eGFR[,c(1,7,5,6)]
eGFR$trait_GWAS <- "eGFR"
head(eGFR)
colnames(eGFR) <- c("SNPID", "GWAS_P", "GWAS_BETA", "GWAS_SE", "trait_GWAS")
write.table(eGFR, "./001_glycosuria/reviewer_analysis/eGFR.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## rbind 
SNPID_traits <- rbind(eGFR, diagram_T2D)
head(SNPID_traits)
dim(SNPID_traits)

SNP <- rbind(magic_fg, magic_hba1c, giant_BMI)
head(SNP)
dim(SNP)

## glycosuria GWAS SNPs
data <- fread("./001_glycosuria/GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final.txt", header = T, sep = "\t")
head(data)
dim(data)
data$SNPID <- paste0(data$chr, ":", data$position)
head(data)
data <- data[,c(1,19,16,17,18)]
colnames(data) <- c("SNP", "SNPID", "glycosuria_P", "glycosuria_BETA", "glycosuria_SE")

## join
SNPID_traits <- left_join(SNPID_traits, data, by = "SNPID")
head(SNPID_traits)
dim(SNPID_traits)
SNPID_traits <- SNPID_traits[,c(6,1:5,7:9)]

SNP <- left_join(SNP, data, by = "SNP")
head(SNP)
dim(SNP)

final_data <- rbind(SNP, SNPID_traits)
head(final_data)
dim(final_data)
write.table(final_data, "./001_glycosuria/reviewer_analysis/associations_from_other_traits_and_glycosuria_estimates.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


final_data <- read.table("associations_from_other_traits_and_glycosuria_estimates.txt", header = T, sep = "\t")

final_data$direction <- 
  ifelse(final_data$GWAS_BETA < 0 & final_data$glycosuria_BETA < 0, "Same", 
         ifelse(final_data$GWAS_BETA > 0 & final_data$glycosuria_BETA > 0, "Same", 
                ifelse(final_data$GWAS_BETA > 0 & final_data$glycosuria_BETA < 0, "Opposite", 
                       ifelse(final_data$GWAS_BETA < 0 & final_data$glycosuria_BETA > 0, "Opposite", "Same"))))
table(final_data$direction)
table(final_data$trait_GWAS)
str(final_data)
library(plyr)
final_data$trait_GWAS <- revalue(final_data$trait_GWAS, c("diagram_T2D"="Type 2 diabetes", "eGFR"="eGFR", "giant_BMI"="BMI","magic_fg"="Fasting glucose", "magic_hba1c"="HbA1c"))


source("ggplot_my_theme.R")

head(final_data)
library(ggplot2)
library(ggpubr)
final_data <- na.omit(final_data)
png("manuscript/HMG_submission/revision/revision_figures/scatter_plot.png",
    height = 300,
    width = 300,
    units = "mm",
    res = 1000)
p <- ggplot(data = final_data, 
       aes(x = glycosuria_BETA, y = GWAS_BETA)) +
  geom_point(aes(x = glycosuria_BETA, y = GWAS_BETA)) +
  geom_smooth(method = "lm",
              alpha = 1,
              show.legend = FALSE,
              se = F) +  
  facet_grid(trait_GWAS ~ .) +
  my_theme() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept=0) +
  labs(colour = "Direction of effect estimates") +
  xlab("Glycosuria effect estimate") +
  ylab("GWAS effect estimate")
print(p)

p + stat_cor(method = "pearson", label.x = -0.2, label.y = 0.4, label.sep = "")
dev.off()

# Pearson correlation coefficient for each trait GWAS beta and glycosuria beta
table(final_data$trait_GWAS)
t2d <- subset(final_data, trait_GWAS == "Type 2 diabetes")
egfr <- subset(final_data, trait_GWAS == "eGFR")
bmi <- subset(final_data, trait_GWAS == "BMI")
fg <- subset(final_data, trait_GWAS == "Fasting glucose")
hba1c <- subset(final_data, trait_GWAS == "HbA1c")

cor.test(t2d$GWAS_BETA, t2d$glycosuria_BETA)
cor.test(egfr$GWAS_BETA, egfr$glycosuria_BETA)
cor.test(bmi$GWAS_BETA, bmi$glycosuria_BETA)
cor.test(fg$GWAS_BETA, fg$glycosuria_BETA)
cor.test(hba1c$GWAS_BETA, hba1c$glycosuria_BETA)






# Reviewer comment: 9.1 Could MR analyses test whether fasting glucose or eGFR causally influences glycosuria in pregnancy? I think these are currently missed opportunities. ====
library(MRInstruments) 
library(TwoSampleMR)
library(data.table)

## source environment/data etc. ====
### MR Base data
ao <- available_outcomes(access_token=NULL)

### methods
methods <- mr_method_list()
methods_heterogeneity <- subset(methods, heterogeneity_test == TRUE)
methods_heterogeneity <- subset(methods_heterogeneity, obj != "mr_ivw_radial")
methods <- subset(methods, obj != "mr_ivw_radial")
methods <- subset(methods, obj != "mr_raps")
methods <- methods$obj
methods_heterogeneity <- methods_heterogeneity$obj

#install_github("WSpiller/RadialMR")
library(RadialMR)
#install_github("qingyuanzhao/mr.raps")
library(mr.raps)

### colours
#install.packages("wesanderson")
library(wesanderson)
d1 <- wes_palette("Royal1", type = "discrete")
d2 <- wes_palette("GrandBudapest2", type = "discrete")
d3 <- wes_palette("Cavalcanti1", type = "discrete")
d4 <- wes_palette("Rushmore1", type = "discrete")
discrete_wes_pal <- c(d1, d2, d3, d4)
rm(d1,d2,d3,d4)

### source other scripts
source("my_mr_scatter_plot.R")

## extract exposure instruments ====
exposure_data <- read_exposure_data("GWAS_output_final_MR_rs13337037.txt",
                                    clump = F,
                                    sep = "\t",
                                    snp_col = "SNP",
                                    beta_col = "BETA",
                                    se_col = "SE",
                                    eaf_col = "EAF",
                                    effect_allele_col = "EA",
                                    other_allele_col = "NEA",
                                    pval_col = "P",
                                    samplesize_col = "all_total",
                                    min_pval = 5e-8)

exposure_data$exposure <- "Glycosuria"
exposure_data$id.exposure <- "Glycosuria"

dim(exposure_data)
head(exposure_data)


## extract outcome data ====
outcome_data <- read_outcome_data(
  snps = exposure_data$SNP,
  filename = "Maternal_Effect_European_meta_NG2019.txt",
  sep = " ",
  snp_col = "RSID",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "ea",
  other_allele_col = "nea",
  eaf_col = "eaf",
  pval_col = "p"
)

write.table(outcome_data, "offspring_bw_egg.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


outcome_data_egg_fetal <- read.table("offspring_bw_egg_fetalstats.txt", header = T, sep = "\t")
outcome_data_egg_maternal <- read.table("offspring_bw_egg_maternalstats.txt", header = T, sep = "\t")
outcome_data_egg <- read.table("egg_bw.txt", header = T, sep= "\t")

outcome_data_t2d <- extract_outcome_data(snps = exposure_data$SNP, outcomes = c(23,24,25,26,976,1090), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token = NULL)

## harmonize data ====
harmonise_data <- harmonise_data(exposure_data, outcome_data, action = 2)
harmonise_data_egg <- harmonise_data(exposure_data, outcome_data_egg, action = 2)
harmonise_data_t2d <- harmonise_data(exposure_data, outcome_data_t2d, action = 2)

## MR ====
mr_results <- mr(harmonise_data)
mr_results_egg <- mr(harmonise_data_egg)
mr_results_t2d <- mr(harmonise_data_t2d)



