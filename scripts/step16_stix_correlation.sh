# association between self-report glycosuria and reagent strip glycosuria
# association between lead SNPs and reagent-strip glycosuria

cd ./other_analysis/regressions

R

data <- read.table("./other_analysis/regressions/pheno_file.txt", header = T)

head(data)
dim(data)

# association between self-report glycosuria and reagent strip glycosuria
chisq <- chisq.test(data$e113, data$stix)
chisq
chisq[["observed"]]
sqrt(chisq.test(data$e113, data$stix)$statistic/length(data$e113))

# association between lead SNPs and reagent-strip glycosuria
chisq <- chisq.test(data$rs13337037_dosage, data$stix)
chisq
sqrt(chisq.test(data$rs13337037_dosage, data$stix)$statistic/length(data$e113))

chisq <- chisq.test(data$rs10991823_dosage, data$stix)
chisq
sqrt(chisq.test(data$rs10991823_dosage, data$stix)$statistic/length(data$e113))

