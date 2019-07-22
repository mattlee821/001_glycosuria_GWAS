#################################################################################
######### correlations for SNP rs10991823

R
library("dplyr")

data <- read.table("./other_analysis/regressions/pheno_file.txt", header = T)

head(data)
dim(data)

### stix correlation
# phenotype
phenotype <- "Stix test (yes/no)"

# n
table(data$stix)
n <- 6639
n_cases <- 227 

# regression
glm <- glm(data$stix ~ data[,"rs10991823_dosage"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

# P
P <- glm_coef[2,4]

#OR
b <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

#add to table
table <- data.frame("Trait" = phenotype, "n" = n, "n_cases" = n_cases, "OR" = b, "CI lower" = ci_lower, "CI upper" = ci_upper, "P" = P)
dim(table)
head(table)
table$Trait <- as.character(table$Trait)

###############################################
### bmi correlation
head(data)
# phenotype
phenotype <- "Pre-pregnancy BMI (kg/m2)"

# n
summary(data$dw042)
n <- 6866 - 416
n_cases <- NA

# regression
glm <- glm(data$dw042 ~ data[,"rs10991823_dosage"]) #, family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)
str(table)

###############################################
### water_drinking_freq_j - sum of j217 and j218 
head(data)
# phenotype
phenotype <- "Frequency of water drunk in one week"

# n
summary(data$water_drinking_freq_j)
n <- 6866 - 1330
n_cases <- NA

# regression
glm <- glm(data$water_drinking_freq_j ~ data[,"rs10991823_dosage"]) #, family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)

###############################################
# n1010	thrush ever - A1k: Mother has ever had thrush	
head(data)
# phenotype
phenotype <- "Mother has ever had thrush"

# n
table(data$n1010)
n <- 1211
n_cases <- 3327

# regression
glm <- glm(data$n1010 ~ data[,"rs10991823_dosage"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

#OR
b <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)

###############################################
# e108	thursh - Thrush in PREG from 7MTHS on
head(data)
# phenotype
phenotype <- "Mother had thrush in pregnancy from 7 months onwards"

# n
table(data$e108)
n <- 5615
n_cases <- 744

# regression
glm <- glm(data$e108 ~ data[,"rs10991823_dosage"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

#OR
b <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)

###############################################
# n1009	urinary infection ever - A1j: Mother has ever had a urinary infection, cystitis or pyelitis	
head(data)
# phenotype
phenotype <- "Mother has ever had a urinary infection, cystitis or pyelitis"

# n
table(data$n1009)
n <- 1735
n_cases <- 2751

# regression
glm <- glm(data$n1009 ~ data[,"rs10991823_dosage"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

#OR
b <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)

###############################################
# e105	urinary infection in preg from 7 months on - Urinary infection in PREG from 7MTHS on	
head(data)
# phenotype
phenotype <- "Mother had urinary infection in pregnancy from 7 months onwards"

# n
table(data$e105)
n <- 6048
n_cases <- 341

# regression
glm <- glm(data$e105 ~ data[,"rs10991823_dosage"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

#OR
b <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)

###############################################
# r3001	urinate at night - C1b: Frequency respondent urinates throughout the night	
head(data)
# phenotype
phenotype <- "Frequency respondent urinates throught the night"

# n
summary(data$r3001)
n <- 6866 - 2337
n_cases <- NA

# regression
glm <- glm(data$r3001 ~ data[,"rs10991823_dosage"]) #, family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)

###############################################
# d046 	hypertension ever - History of hypertension	
head(data)
# phenotype
phenotype <- "Mother has history of hypertension"

# n
table(data$d046)
n <- 5766
n_cases <- 966

# regression
glm <- glm(data$d046 ~ data[,"rs10991823_dosage"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

#OR
b <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)

###############################################
## d047	hypertension in preg - Hypertension in PREG only 
head(data)
# phenotype
phenotype <- "Mother had hypertension during pregnancy only"

# n
table(data$d047)
n <- 272
n_cases <- 687

# regression
glm <- glm(data$d047 ~ data[,"rs10991823_dosage"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

#OR
b <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)

###############################################
## gesthyp - gestational hypertension
head(data)
# phenotype
phenotype <- "Mother had gestational hypertension during pregnancy"

# n
table(data$gesthyp)
n <- 5762
n_cases <- 966

# regression
glm <- glm(data$gesthyp ~ data[,"rs10991823_dosage"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

#OR
b <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)


###############################################
## preeclampsia - 
head(data)
# phenotype
phenotype <- "Mother had preeclampsia during pregnancy"

# n
table(data$preeclampsia)
n <- 6728
n_cases <- 138

# regression
glm <- glm(data$preeclampsia ~ data[,"rs10991823_dosage"], family=binomial(link="logit"))
glm_coef <- summary(glm)$coefficients
glm_coef
summary(glm_coef)

# ci
ci_lower <- glm_coef[2,1] - 1.96 * glm_coef[2,2] # 0.3228759
ci_upper <- glm_coef[2,1] + 1.96 * glm_coef[2,2] # 0.6751995

#b/log_odds/OR
b <- glm_coef[2,1]

#OR
b <- exp(glm_coef[2,1])
ci_lower <- exp(ci_lower)
ci_upper <- exp(ci_upper)

# P
P <- glm_coef[2,4]

#add to table
table <- (rbind(table,list(phenotype, n, n_cases, b, ci_lower, ci_upper, P)))
head(table)

#############
# check table 
dim(table)
# save table
write.table(table, "./other_analysis/regressions/regression_rs10991823_to_trait.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

