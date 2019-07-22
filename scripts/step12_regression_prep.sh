## regressions

## do 0,1,2 for the lead SNP on glycosuria and you should get the same as the beta for GWAS
	# use the dosage files 

## phenotypes
# ID_1
# ID_2
# missing
# father
# mother
# sex
# plink_pheno
# genotype_rs13337037 - dosage
# e113 glycosuria from questionnaire - Sugar in urine in PREG from 7MTHS on 
# stix 	glycosuria from stix test - Glycosuria (2+ or more) on two occasions
# rs13337037_dosage - genotype_rs13337037 dosage for lead SNP
# dw042 BMI
# underweight
# normalweight
# overweight
# obese1
# obese2
# j217	bottled water - FREQ MUM drinks Bottled Drinks	
# j218	tap water - FREQ MUM drinks Tap Water	
# n7074	bottled water - H4e: Frequency mother drinks bottled water	
# n7075	tap water - H4f: Frequency mother drinks tap water	
# water_drinking_freq_j - sum of j217 and j218
# water_drinking_freq_n - sum of n7074 and n7075
# n1010	thrush ever - A1k: Mother has ever had thrush	
# e108	thursh - Thrush in PREG from 7MTHS on
# r3001	urinate at night - C1b: Frequency respondent urinates throughout the night	
# n1009	urinary infection ever - A1j: Mother has ever had a urinary infection, cystitis or pyelitis	
# e105	urinary infection in preg from 7 months on - Urinary infection in PREG from 7MTHS on	
# d046 	hypertension ever - History of hypertension	
# d047	hypertension in preg - Hypertension in PREG only	


## create file with genotype_rs13337037 data
module add apps/qctool-2.0
qctool -g ./GWAS/ALSPAC/step4_filter_SNPs/output/step2_SNP2KEEP_filter/output/data_chr16_sample_filter_subsetted.bgen -s ./GWAS/ALSPAC/step2_subsample_genetic_data/output/data.sample -incl-range :31478711-31478711  -assume-chromosome 16 -ofiletype dosage -og ./other_analysis/regressions/rs13337037.dosage
qctool -g ./GWAS/ALSPAC/step4_filter_SNPs/output/step2_SNP2KEEP_filter/output/data_chr09_sample_filter_subsetted.bgen -s ./GWAS/ALSPAC/step2_subsample_genetic_data/output/data.sample -incl-range :93957251-93957251  -assume-chromosome 9 -ofiletype dosage -og ./other_analysis/regressions/rs10991823.dosage

cd ./other_analysis/regressions

# transpose dosage files and ensure each aln_qlet has correct dosage 
R
data <- read.table("./other_analysis/regressions/rs13337037.dosage", header = F, fill = T, sep = "", quote = "")
head(data)
dim(data)
head(data)[,1:6]
tail(data)
data2 <- t(data)
head(data2, 10)
tail(data2, 10)
colnames(data2) <- c("aln_qlet", "rs13337037_dosage")

write.table(data2, "./other_analysis/regressions/rs13337037_final.dosage", 
            row.names = FALSE, col.names = T, quote = FALSE, sep = "\t", )
  ## because i dont know how to shit values in a column i open the dosage file in excel and do it manually (i shift the dosage column down by one)

data <- read.table("./other_analysis/regressions/rs10991823.dosage", header = F, fill = T, sep = "", quote = "")
head(data)
dim(data)
head(data)[,1:6]
tail(data)
data2 <- t(data)
head(data2, 10)
tail(data2, 10)
colnames(data2) <- c("ID_1", "rs10991823_dosage")

write.table(data2, "./other_analysis/regressions/rs10991823_final.dosage", 
            row.names = FALSE, col.names = T, quote = FALSE, sep = "\t", )
  ## because i dont know how to shit values in a column i open the dosage file in excel and do it manually (i shift the dosage column down by one)


#####################################################################################################################################
cd ./other_analysis/regressions

R

library(foreign)
library(readstata13)
library(dplyr)

# create new phenofile
data <- read.table("./GWAS/ALSPAC/step1_phenofile/pheno_file_include-list-only_regressions.txt", header = T)
head(data)
head(data)
dim(data)

# remove unnecessary columns
names(data)
data <- data[,c(1,12,13,14)]
head(data)

## add genotype_rs13337037 data for each individual
genotype_rs13337037 <- read.table("./other_analysis/regressions/rs13337037_final.dosage", header = T)
head(genotype_rs13337037)
nrow(genotype_rs13337037)
## remove dupliates
inds = unique(genotype_rs13337037[,1])
w = sapply(inds, function(x){   
  a = which(genotype_rs13337037[,1] %in% x)[1]
  return(a)
})
nrow(genotype_rs13337037)

### make a new dataframe
genotype_rs13337037 = genotype_rs13337037[w,]
rownames(genotype_rs13337037) = genotype_rs13337037[,1]
w = which(table(genotype_rs13337037[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(genotype_rs13337037[,1])[w])
w = which(genotype_rs13337037[,1] %in% n)
genotype_rs13337037[w,]

#summary(stix$stix_glycosuria)
summary(genotype_rs13337037$rs13337037_dosage)

#join
colnames(genotype_rs13337037)[1] <- "ID_1"
rownames(genotype_rs13337037) = genotype_rs13337037[,1]
data <- left_join(data, genotype_rs13337037, by = "ID_1", all = T)
head(data)
nrow(data)
str(data)
summary(data$rs13337037_dosage)


## add genotype_rs10991823 data for each individual
genotype_rs10991823 <- read.table("./other_analysis/regressions/rs10991823_final.dosage", header = T)
head(genotype_rs10991823)
nrow(genotype_rs10991823)
## remove dupliates
inds = unique(genotype_rs10991823[,1])
w = sapply(inds, function(x){   
  a = which(genotype_rs10991823[,1] %in% x)[1]
  return(a)
})
nrow(genotype_rs10991823)

### make a new dataframe
genotype_rs10991823 = genotype_rs10991823[w,]
rownames(genotype_rs10991823) = genotype_rs10991823[,1]
w = which(table(genotype_rs10991823[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(genotype_rs10991823[,1])[w])
w = which(genotype_rs10991823[,1] %in% n)
genotype_rs10991823[w,]

#summary(stix$stix_glycosuria)
summary(genotype_rs10991823$rs10991823_dosage)

#join
colnames(genotype_rs10991823)[1] <- "ID_1"
rownames(genotype_rs10991823) = genotype_rs10991823[,1]
data <- left_join(data, genotype_rs10991823, by = "ID_1", all = T)
head(data)
nrow(data)
str(data)
summary(data$rs10991823_dosage)

#####################################

# check data
names(data)
colnames(data)[3:4] <- c("e113", "stix")
head(data, 20)
dim(data) #6866
str(data)

data$pregnancy_diabetes <- as.character(data$pregnancy_diabetes)
data$pregnancy_diabetes <- as.numeric(data$pregnancy_diabetes)
summary(data$pregnancy_diabetes)

data$e113 <- as.character(data$e113)
data$e113 <- as.numeric(data$e113)
summary(data$e113)

data$stix <- as.character(data$stix)
data$stix <- as.numeric(data$stix)
summary(data$stix)
str(data)
######################################################################################################################################

## add phenotype data on BMI (actual BMI)
#open data
bmi <- read.dta("./data/d_3c.dta")
bmi <- bmi[,c("aln", "dw042")]
nrow(bmi)
head(bmi)

#organise data
bmi$qlet <- "M"
bmi$aln_qlet <- paste0(bmi$aln, "M")
bmi <- bmi[c("aln_qlet", "dw042")]
summary(bmi$dw042)

## remove dupliates
inds = unique(bmi[,1])
w = sapply(inds, function(x){   
  a = which(bmi[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
bmi = bmi[w,]
rownames(bmi) = bmi[,1]
w = which(table(bmi[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(bmi[,1])[w])
w = which(bmi[,1] %in% n)
bmi[w,]

#summary(stix$stix_glycosuria)
summary(bmi$dw042)

#join
colnames(bmi)[1] <- "ID_1"
rownames(data) = data[,1]
data2 <- left_join(data, bmi, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$dw042)

data2$dw042[data2$dw042 <= 0] <- NA
summary(data2$dw042)
data2$bmi_yes_no[data2$dw042 >= 1] <- 1
data2$bmi_yes_no[data2$dw042 < 1] <- 0
data2$bmi_yes_no[data2$dw042 == NA] <- 0
summary(data2$bmi_yes_no)
table(data2$bmi_yes_no)

######################################################################################################################################


## add phenotype data on amount of water drunk during the week

# bottled water j217
#open data
bottled_water <- read.dta13("./data/j_5a.dta")
bottled_water <- bottled_water[,c("aln", "j217")]
nrow(bottled_water)
head(bottled_water)
summary(bottled_water$j217)
table(bottled_water$j217)

#organise data
bottled_water$qlet <- "M"
bottled_water$aln_qlet <- paste0(bottled_water$aln, "M")
bottled_water <- bottled_water[c("aln_qlet", "j217")]
summary(bottled_water$j217)

## remove dupliates
inds = unique(bottled_water[,1])
w = sapply(inds, function(x){   
  a = which(bottled_water[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
bottled_water = bottled_water[w,]
rownames(bottled_water) = bottled_water[,1]
w = which(table(bottled_water[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(bottled_water[,1])[w])
w = which(bottled_water[,1] %in% n)
bottled_water[w,]

#summary
table(bottled_water$j217)

#join
colnames(bottled_water)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, bottled_water, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$j217)
table(data2$j217)

#recode values
data2$j217[data2$j217 == -1] <- NA #missing data
table(data2$j217)
#####################################

# bottled water n7074
#open data
bottled_water <- read.dta13("./data/n_2a.dta")
bottled_water <- bottled_water[,c("aln", "n7074")]
nrow(bottled_water)
head(bottled_water)
summary(bottled_water$n7074)
table(bottled_water$n7074)

#organise data
bottled_water$qlet <- "M"
bottled_water$aln_qlet <- paste0(bottled_water$aln, "M")
bottled_water <- bottled_water[c("aln_qlet", "n7074")]
summary(bottled_water$n7074)

## remove dupliates
inds = unique(bottled_water[,1])
w = sapply(inds, function(x){   
  a = which(bottled_water[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
bottled_water = bottled_water[w,]
rownames(bottled_water) = bottled_water[,1]
w = which(table(bottled_water[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(bottled_water[,1])[w])
w = which(bottled_water[,1] %in% n)
bottled_water[w,]

#summary
table(bottled_water$n7074)

#join
colnames(bottled_water)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, bottled_water, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$n7074)
table(data2$n7074)

#recode values
data2$n7074[data2$n7074 == -11] <- NA #missing data
data2$n7074[data2$n7074 == -10] <- NA #missing data
data2$n7074[data2$n7074 == -1] <- NA #missing data
table(data2$n7074)
#####################################

# tap water
#open data
tap_water <- read.dta13("./data/j_5a.dta")
tap_water <- tap_water[,c("aln", "j218")]
nrow(tap_water)
head(tap_water)
summary(tap_water$j218)
table(tap_water$j218)

#organise data
tap_water$qlet <- "M"
tap_water$aln_qlet <- paste0(tap_water$aln, "M")
tap_water <- tap_water[c("aln_qlet", "j218")]
summary(tap_water$j218)

## remove dupliates
inds = unique(tap_water[,1])
w = sapply(inds, function(x){   
  a = which(tap_water[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
tap_water = tap_water[w,]
rownames(tap_water) = tap_water[,1]
w = which(table(tap_water[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(tap_water[,1])[w])
w = which(tap_water[,1] %in% n)
tap_water[w,]

#summary
table(tap_water$j218)

#join
colnames(tap_water)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, tap_water, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$j218)
table(data2$j218)

#recode values
data2$j218[data2$j218 == -1] <- NA #missing data
table(data2$j218)
#####################################

# tap water
#open data
tap_water <- read.dta13("./data/n_2a.dta")
tap_water <- tap_water[,c("aln", "n7075")]
nrow(tap_water)
head(tap_water)
summary(tap_water$n7075)
table(tap_water$n7075)

#organise data
tap_water$qlet <- "M"
tap_water$aln_qlet <- paste0(tap_water$aln, "M")
tap_water <- tap_water[c("aln_qlet", "n7075")]
summary(tap_water$n7075)

## remove dupliates
inds = unique(tap_water[,1])
w = sapply(inds, function(x){   
  a = which(tap_water[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
tap_water = tap_water[w,]
rownames(tap_water) = tap_water[,1]
w = which(table(tap_water[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(tap_water[,1])[w])
w = which(tap_water[,1] %in% n)
tap_water[w,]

#summary
table(tap_water$n7075)

#join
colnames(tap_water)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, tap_water, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$n7075)


table(data2$n7075)

#recode values
data2$n7075[data2$n7075 == -11] <- NA #missing data
data2$n7075[data2$n7075 == -10] <- NA #missing data
data2$n7075[data2$n7075 == -1] <- NA #missing data
table(data2$n7075)
#####################################


#####################################
# sum the two water phenotypes
# water_drinking_freq_j - sum of j217 and j218
data2$water_drinking_freq_j <- data2$j217 + data2$j218
table(data2$water_drinking_freq_j)
# water-drinking_freq_n - sum of n7074 and n7075
data2$water_drinking_freq_n <- data2$n7074 + data2$n7075
table(data2$water_drinking_freq_n)

head(data2)
######################################################################################################################################

## add in phenotype data on thrush

# n1010 thrush - A1k: Mother has ever had thrush	
thrush_ever <- read.dta13("./data/n_2a.dta")
thrush_ever <- thrush_ever[,c("aln", "n1010")]
nrow(thrush_ever)
head(thrush_ever)
summary(thrush_ever$n1010)
table(thrush_ever$n1010)

#organise data
thrush_ever$qlet <- "M"
thrush_ever$aln_qlet <- paste0(thrush_ever$aln, "M")
thrush_ever <- thrush_ever[c("aln_qlet", "n1010")]
summary(thrush_ever$n1010)

## remove dupliates
inds = unique(thrush_ever[,1])
w = sapply(inds, function(x){   
  a = which(thrush_ever[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
thrush_ever = thrush_ever[w,]
rownames(thrush_ever) = thrush_ever[,1]
w = which(table(thrush_ever[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(thrush_ever[,1])[w])
w = which(thrush_ever[,1] %in% n)
thrush_ever[w,]

#summary
table(thrush_ever$n1010)

#join
colnames(thrush_ever)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, thrush_ever, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$n1010)
table(data2$n1010)

#recode values
data2$n1010[data2$n1010 == -11] <- NA #missing data
data2$n1010[data2$n1010 == -10] <- NA #missing data
data2$n1010[data2$n1010 == -1] <- NA #missing data
data2$n1010[data2$n1010 == 3] <- NA #missing data
data2$n1010[data2$n1010 == 2] <- 0 
table(data2$n1010)
#####################################

# e108  thursh - Thrush in PREG from 7MTHS on
thrush_preg <- read.dta("./data/e_4d.dta")
thrush_preg <- thrush_preg[,c("aln", "e108")]
nrow(thrush_preg)
head(thrush_preg)
summary(thrush_preg$e108)
table(thrush_preg$e108)

#organise data
thrush_preg$qlet <- "M"
thrush_preg$aln_qlet <- paste0(thrush_preg$aln, "M")
thrush_preg <- thrush_preg[c("aln_qlet", "e108")]
summary(thrush_preg$e108)

## remove dupliates
inds = unique(thrush_preg[,1])
w = sapply(inds, function(x){   
  a = which(thrush_preg[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
thrush_preg = thrush_preg[w,]
rownames(thrush_preg) = thrush_preg[,1]
w = which(table(thrush_preg[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(thrush_preg[,1])[w])
w = which(thrush_preg[,1] %in% n)
tap_water[w,]

#summary
table(thrush_preg$e108)
str(thrush_preg)
thrush_preg$e108 = as.character(thrush_preg$e108)

#recode values
thrush_preg$e108[thrush_preg$e108 == "Missing"] <- NA
#recode 'No' to 0
thrush_preg$e108[thrush_preg$e108 == "N"] <- 0
#recode 'Yes' to 1
thrush_preg$e108[thrush_preg$e108 == "Y"] <- 1

thrush_preg$e108 = as.numeric(thrush_preg$e108)

#join
colnames(thrush_preg)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, thrush_preg, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$e108)
table(data2$e108)
######################################################################################################################################

## add in phenotype data on urinating more often than normal, particulalry at night
# r3001 urinate at night - C1b: Frequency respondent urinates throughout the night	
r3001 <- read.dta("./data/r_r1b.dta")
r3001 <- r3001[,c("aln", "r3001")]
nrow(r3001)
head(r3001)
summary(r3001$r3001)
table(r3001$r3001)

#organise data
r3001$qlet <- "M"
r3001$aln_qlet <- paste0(r3001$aln, "M")
r3001 <- r3001[c("aln_qlet", "r3001")]
summary(r3001$r3001)

## remove dupliates
inds = unique(r3001[,1])
w = sapply(inds, function(x){   
  a = which(r3001[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
r3001 = r3001[w,]
rownames(r3001) = r3001[,1]
w = which(table(r3001[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(r3001[,1])[w])
w = which(r3001[,1] %in% n)
tap_water[w,]

#summary
table(r3001$r3001)
str(r3001)

#join
colnames(r3001)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, r3001, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$r3001)
table(data2$r3001)

#recode values
data2$r3001 = as.character(data2$r3001)
data2$r3001[data2$r3001 == "Once"] <- 1
data2$r3001[data2$r3001 == "Twice"] <- 2
data2$r3001[data2$r3001 == "Three times"] <- 3
data2$r3001[data2$r3001 == "Four times or more"] <- 4
data2$r3001[data2$r3001 == "None"] <- 0
data2$r3001[data2$r3001 == "No response "] <- "NA"
data2$r3001[data2$r3001 == "Not completed"] <- "NA"
data2$r3001[data2$r3001 == "Triplet / quadruplet"] <- "NA"
summary(data2$r3001)
table(data2$r3001)
data2$r3001 = as.numeric(data2$r3001)
summary(data2$r3001)
table(data2$r3001)
######################################################################################################################################

## add in phenotype data on urinary infection 
# n1009 urinary infection ever - A1j: Mother has ever had a urinary infection, cystitis or pyelitis	
n1009 <- read.dta13("./data/n_2a.dta")
n1009 <- n1009[,c("aln", "n1009")]
nrow(n1009)
head(n1009)
summary(n1009$n1009)
table(n1009$n1009)

#organise data
n1009$qlet <- "M"
n1009$aln_qlet <- paste0(n1009$aln, "M")
n1009 <- n1009[c("aln_qlet", "n1009")]
summary(n1009$n1009)

## remove dupliates
inds = unique(n1009[,1])
w = sapply(inds, function(x){   
  a = which(n1009[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
n1009 = n1009[w,]
rownames(n1009) = n1009[,1]
w = which(table(n1009[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(n1009[,1])[w])
w = which(n1009[,1] %in% n)
tap_water[w,]

#summary
table(n1009$n1009)
str(n1009)

#join
colnames(n1009)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, n1009, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$n1009)
table(data2$n1009)

#recode values
data2$n1009[data2$n1009 == -11] <- NA #missing data
data2$n1009[data2$n1009 == -10] <- NA #missing data
data2$n1009[data2$n1009 == -1] <- NA #missing data
data2$n1009[data2$n1009 == 3] <- NA #missing data
data2$n1009[data2$n1009 == 2] <- 0 #missing data
table(data2$n1010)
######################################################################################################################################

# e105 urinary infection in preg from 7 months on - Urinary infection in PREG from 7MTHS on	
e105 <- read.dta("./data/e_4d.dta")
e105 <- e105[,c("aln", "e105")]
nrow(e105)
head(e105)
summary(e105$e105)
table(e105$e105)

#organise data
e105$qlet <- "M"
e105$aln_qlet <- paste0(e105$aln, "M")
e105 <- e105[c("aln_qlet", "e105")]
summary(e105$e105)

## remove dupliates
inds = unique(e105[,1])
w = sapply(inds, function(x){   
  a = which(e105[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
e105 = e105[w,]
rownames(e105) = e105[,1]
w = which(table(e105[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(e105[,1])[w])
w = which(e105[,1] %in% n)
tap_water[w,]

#summary
table(e105$e105)
str(e105)

#join
colnames(e105)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, e105, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$e105)
table(data2$e105)

#recode values
data2$e105 = as.character(data2$e105)
data2$e105[data2$e105 == "Y"] <- 1
data2$e105[data2$e105 == "N"] <- 0
data2$e105[data2$e105 == "Missing"] <- "NA"
data2$e105[data2$e105 == "N"] <- 0
data2$e105 = as.numeric(data2$e105)
summary(data2$e105)
table(data2$e105)
######################################################################################################################################

## add phenotype data on hypertension 

# d046 hypertension ever - History of hypertension	
hyper_preg <- read.dta("./data/d_3c.dta")
hyper_preg <- hyper_preg[,c("aln", "d046")]
nrow(hyper_preg)
head(hyper_preg)

#organise data
hyper_preg$qlet <- "M"
hyper_preg$aln_qlet <- paste0(hyper_preg$aln, "M")
hyper_preg <- hyper_preg[c("aln_qlet", "d046")]
summary(hyper_preg$d046)

## remove dupliates
inds = unique(hyper_preg[,1])
w = sapply(inds, function(x){   
  a = which(hyper_preg[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
hyper_preg = hyper_preg[w,]
rownames(hyper_preg) = hyper_preg[,1]
w = which(table(hyper_preg[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(hyper_preg[,1])[w])
w = which(hyper_preg[,1] %in% n)
hyper_preg[w,]

#summary
table(hyper_preg$d046)

#join
colnames(hyper_preg)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, hyper_preg, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$d046)
table(data2$d046)

#recode values
data2$d046[data2$d046 == -1] <- NA #missing data
data2$d046[data2$d046 == 2] <- 0 #missing data
table(data2$d046)
#####################################

# d047 hypertension in preg only - Hypertension in PREG only	
#open data
hyper_preg <- read.dta("./data/d_3c.dta")
hyper_preg <- hyper_preg[,c("aln", "d047")]
nrow(hyper_preg)
head(hyper_preg)

#organise data
hyper_preg$qlet <- "M"
hyper_preg$aln_qlet <- paste0(hyper_preg$aln, "M")
hyper_preg <- hyper_preg[c("aln_qlet", "d047")]
summary(hyper_preg$d047)

## remove dupliates
inds = unique(hyper_preg[,1])
w = sapply(inds, function(x){   
  a = which(hyper_preg[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
hyper_preg = hyper_preg[w,]
rownames(hyper_preg) = hyper_preg[,1]
w = which(table(hyper_preg[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(hyper_preg[,1])[w])
w = which(hyper_preg[,1] %in% n)
hyper_preg[w,]

#summary
table(hyper_preg$d047)

#join
colnames(hyper_preg)[1] <- "ID_1"
rownames(data2) = data2[,1]
data2 <- left_join(data2, hyper_preg, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$d047)
table(data2$d047)

#recode values
data2$d047[data2$d047 == -1] <- NA #missing data
data2$d047[data2$d047 == -2] <- NA #missing data
data2$d047[data2$d047 == 2] <- 0 #missing data
table(data2$d047)
#####################################

# gesthyp - gestational hypertenstion	
#open data
gesthyp <- read.dta("./data/OA_r1b.dta")
gesthyp <- gesthyp[,c("aln", "gesthyp")]
nrow(gesthyp)
head(gesthyp)

#organise data
gesthyp$qlet <- "M"
gesthyp$aln_qlet <- paste0(gesthyp$aln, "M")
gesthyp <- gesthyp[c("aln_qlet", "gesthyp")]
summary(gesthyp$gesthyp)

## remove dupliates
inds = unique(gesthyp[,1])
w = sapply(inds, function(x){   
  a = which(gesthyp[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
gesthyp = gesthyp[w,]
rownames(gesthyp) = gesthyp[,1]
w = which(table(gesthyp[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(gesthyp[,1])[w])
w = which(gesthyp[,1] %in% n)
gesthyp[w,]

#summary
table(gesthyp$gesthyp)

#join
colnames(gesthyp)[1] <- "ID_1"
rownames(gesthyp) = gesthyp[,1]
data2 <- left_join(data2, gesthyp, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$gesthyp)
table(data2$gesthyp)

#recode values
data2$gesthyp <- as.character(data2$gesthyp)
data2$gesthyp[data2$gesthyp == "No"] <- 0 
data2$gesthyp[data2$gesthyp == "Yes"] <- 1
table(data2$gesthyp)
#####################################

preeclampsia <- read.dta("./data/OA_r1b.dta")
preeclampsia <- preeclampsia[,c("aln", "preeclampsia")]
nrow(preeclampsia)
head(preeclampsia)
table(preeclampsia$preeclampsia)

#organise data
preeclampsia$qlet <- "M"
preeclampsia$aln_qlet <- paste0(preeclampsia$aln, "M")
preeclampsia <- preeclampsia[c("aln_qlet", "preeclampsia")]
summary(preeclampsia$preeclampsia)

## remove dupliates
inds = unique(preeclampsia[,1])
w = sapply(inds, function(x){   
  a = which(preeclampsia[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
preeclampsia = preeclampsia[w,]
rownames(preeclampsia) = preeclampsia[,1]
w = which(table(preeclampsia[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(preeclampsia[,1])[w])
w = which(preeclampsia[,1] %in% n)
preeclampsia[w,]

#summary
table(preeclampsia$preeclampsia)

#join
colnames(preeclampsia)[1] <- "ID_1"
rownames(preeclampsia) = preeclampsia[,1]
data2 <- left_join(data2, preeclampsia, by = "ID_1", all = T)
head(data2)
nrow(data2)
str(data2)
summary(data2$preeclampsia)
table(data2$preeclampsia)

#recode values
data2$preeclampsia <- as.character(data2$preeclampsia)
data2$preeclampsia[data2$preeclampsia == "No"] <- 0 
data2$preeclampsia[data2$preeclampsia == "Yes"] <- 1
table(data2$preeclampsia)
#################################################################################

head(data2)

write.table(data2, "./other_analysis/regressions/pheno_file.txt", 
	row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
#################################################################################

#################################################################################
cd ./other_analysis/regressions

R

data <- read.table("./other_analysis/regressions/pheno_file.txt", header = T)
head(data)
dim(data)

## check all phenotypes
str(data)

# e113 glycosuria from questionnaire - Sugar in urine in PREG from 7MTHS on 
summary(data$e113)
table(data$e113)

# stix 	glycosuria from stix test - Glycosuria (2+ or more) on two occasions
summary(data$stix)
table(data$stix)

# dw042 BMI
summary(data$dw042)
table(data$dw042)
summary(data$bmi_yes_no)
table(data$bmi_yes_no)

# j217	bottled water - FREQ MUM drinks Bottled Drinks	
summary(data$j217)
table(data$j217)

# n7074	bottled water - H4e: Frequency mother drinks bottled water	
summary(data$n7074)
table(data$n7074)

# j218	tap water - FREQ MUM drinks Tap Water	
summary(data$j218)
table(data$j218)

# n7075	tap water - H4f: Frequency mother drinks tap water	
summary(data$n7075)
table(data$n7075)

# water_drinking_freq_j - sum of j217 and j218
summary(data$water_drinking_freq_j)
table(data$water_drinking_freq_j)

# water_drinking_freq_n - sum of n7074 and n7075
summary(data$water_drinking_freq_n)
table(data$water_drinking_freq_n)

# n1010	thrush ever - A1k: Mother has ever had thrush	
summary(data$n1010)
table(data$n1010)

# e108	thursh - Thrush in PREG from 7MTHS on
summary(data$e108)
table(data$e108)

# r3001	urinate at night - C1b: Frequency respondent urinates throughout the night	
summary(data$r3001)
table(data$r3001)

# n1009	urinary infection ever - A1j: Mother has ever had a urinary infection, cystitis or pyelitis	
summary(data$n1009)
table(data$n1009)

# e105	urinary infection in preg from 7 months on - Urinary infection in PREG from 7MTHS on	
summary(data$e105)
table(data$e105)

# d046 	hypertension ever - History of hypertension	
summary(data$d046)
table(data$d046)

# d047	hypertension in preg - Hypertension in PREG only
summary(data$d047)
table(data$d047)

# gestational hypertension
summary(data$gesthyp)
table(data$gesthyp)

# preeclampsia
summary(data$preeclampsia)
table(data$preeclampsia)
#################################################################################







