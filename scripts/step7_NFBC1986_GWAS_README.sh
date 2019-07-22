# NFBC GWAS conversion to mother genotype effect
# to convert to mothers genotype effect Peter Joshi multiplies the estimate and SE by 2

cd ./GWAS/NFBC

R

library(data.table)

data <- fread("./GWAS/NFBC/GLYCOSURIA_NFBC1986_MODEL1.txt", header = T)
dim(data)
names(data)
head(data)
summary(data$P)
summary(data$BETA)
summary(data$INFO) #18260796
data2 <- subset(data, INFO >= 0.3)
summary(data2$INFO)
dim(data2) #17449460

18260796 - 17449460

summary(data2$EAF)
data2 <- subset(data2, EAF >= 0.01)
summary(data2$EAF)
dim(data2) #9323831

18260796 - 9873828

### final check of data - ensure values make sense
summary(data2$Sample_Size)
summary(data$BETA)

# calculate OR 
data2$OR <- exp(data2$BETA)
data2$ci_lower <- data2$OR - 1.96 * data2$SE
data2$ci_upper <- data2$OR + 1.96 * data2$SE
head(data2)

write.table(data2, "NFBC_GWAS_EAF_INFO.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# calculate parental beta and se
data2$parental_BETA <- data2$BETA * 2
data2$parental_SE <- data2$SE * 2

# calculate parental OR
data2$parental_OR <- exp(data2$parental_BETA)
data2$parental_ci_lower <- data2$parental_OR - 1.96 * data2$parental_SE
data2$parental_ci_upper <- data2$parental_OR + 1.96 * data2$parental_SE

# save final NFBC GWAS file
write.table(data2, "NFBC_GWAS_final.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

###### SNPs reaching GWAS threshol 5e-8
data_sig <- subset(data2, P <= 5e-8)
dim(data_sig) #55
head(data_sig)
data_sig <- data_sig[order(P),] 

write.table(data_sig, "NFBC_GWAS_final_sig.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data <- fread("./GWAS/NFBC/NFBC_GWAS_final.txt", header = T)
dim(data)
head(data)
names(data)

data_parental <- data[,c(1,2,3,4,5,6,7,10,11,12,16,17,18,19,20)] 
names(data_parental)
colnames(data_parental) <- c("SNP", "choromosome", "Position", "EA", "NEA", "Sample_Size", "EAF", "P", "INFO", "HWE_P", "BETA", "SE", "OR", "ci_lower", "ci_upper") 
head(data_parental)

write.table(data_parental, "NFBC_GWAS_final_parental.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")








