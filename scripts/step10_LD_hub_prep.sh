#LD Hub data prep

cd ./GWAS/ALSPAC/step10_LD_hub

#make file LD Hub ready - 
R

library(data.table)
data <- fread("./GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final.txt", header = T)
head(data)

#snpid	A1	A2	Zscore	N	P-value
data_new <- data[,c(1,4,5,7,16,17)]
head(data_new)
data_new <- data_new[,c(1,3,2,6,4,5)]
head(data_new)
colnames(data_new) <- c("snpid", "A1", "A2", "beta", "N", "P-value")
str(data_new)
nrow(data_new)
summary(data_new$beta)

############ calculate Z scores
sd <- sd(data_new$beta, na.rm = T)
sd

mean <- mean(data_new$beta)
mean

data_new$Zscore <- data_new$beta - mean /sd
head(data_new)

data_new <- data_new[,c("snpid", "A1", "A2", "Zscore", "N", "P-value")]
head(data_new)

data_new <- data_new[,c(1,3,2,4,5,6)]
colnames(data_new) <- c("snpid", "A1", "A2", "Zscore", "N", "P-value")
head(data_new)

rs13337037 <- subset(data_new, snpid == "rs13337037")
head(rs13337037)


rs10991823 <- subset(data_new, snpid == "rs10991823")
head(rs10991823)


write.table(data_new, "./GWAS/ALSPAC/step10_LD_hub/LD_hub.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(rs13337037, "./GWAS/ALSPAC/step10_LD_hub/rs13337037_ld_hub.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(rs10991823, "./GWAS/ALSPAC/step10_LD_hub/rs10991823_ld_hub.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# zip files
zip rs13337037_ld_hub.zip rs13337037_ld_hub.txt
zip rs10991823_ld_hub.zip rs10991823_ld_hub.txt
zip ld_hub.zip LD_hub.txt

# download ld_hub.zip and upload it to LD Hub











