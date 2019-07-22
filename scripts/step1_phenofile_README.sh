##########################
## step1
## make phenofile
##########################

# some of this file can only be created once you have made principal components and identified related individuals - if using premade PCs and relatedness then it should work straight away

#this script will prepare a dataset of all ALSPAC individuals who have genetic data as well as whether they have measurements for the phenotypes of interest
#and if they have withdrawn consent or need excluding based on cryptic relatedness. 

#this will output two .txt files, one with no exclusions, and one with exclusions
#it will also output two further .txt files, one with all individuals, and one with individuals to be excluded

#1 is always YES (yes this person has the phenotype, yes this person has withdrawn consent, yes this person needs excluding on relatedness)
#0 is always NO (NO this person does not have the phenotype, NO this person has not withdrawn consent, NO this person does not need excluding on relatedness)
#NA is present where a measurement for the phenotype is not available

cd ./GWAS/ALSPAC/step1_phenofile

R

library(foreign)
library(dplyr)

######################################
## genetic data
######################################

#aln file for all individuals
genetic_data <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/data.sample", skip=2)
colnames(genetic_data) <- c("ID_1","ID_2","missing","father","mother","sex","plink_pheno")
colnames(genetic_data)[1] <- "aln_qlet"

#split genetic_data aln and qlet column into two
genetic_data$qlet <- gsub("[0-9]", "", genetic_data$aln_qlet)
genetic_data$aln <- gsub("[a-zA-Z]", "", genetic_data$aln_qlet)
genetic_data$qlet <- as.factor(genetic_data$qlet)
table(genetic_data$qlet)
head(genetic_data)

######################################
## withdrawl of consent
######################################

#withdrawl of consent file
WOC <- read.csv("./GWAS/ALSPAC/data/ALSPAC_CONSENT_WITHDRAWN_LIST_PLUS_TripQuads.csv", header = T)
WOC <- WOC[4]
WOC$WOC <- 1

#create new dataframe with aln and WOC. where WOC individuals have 1 and everyone else has 0
dataframe <- left_join(genetic_data, WOC, by="aln_qlet")
dataframe$WOC <- as.factor(dataframe$WOC)
summary(dataframe$WOC)
dataframe$WOC <- as.character(as.numeric(dataframe$WOC))
dataframe$WOC[is.na(dataframe$WOC)] <- 0
summary(dataframe$WOC)
table(dataframe$qlet, dataframe$WOC)
head(dataframe)

######################################
## related individuals
######################################
#open and organise data
related_exclusion1 <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/unrelated_ids/mothers/exclusion_list.txt", header = F) #mothers
related_exclusion2 <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/unrelated_ids/children/exclusion_list.txt", header = F) #children
related_exclusion <- rbind(related_exclusion1, related_exclusion2)
rm(related_exclusion1, related_exclusion2)
related_exclusion <- related_exclusion[1]
colnames(related_exclusion) <- "aln_qlet"
related_exclusion$related_exclusion <- 1

#join
dataframe1 <- left_join(dataframe, related_exclusion, by = "aln_qlet", all = T)
dataframe1$related_exclusion[is.na(dataframe1$related_exclusion)] <- 0

#check join
table(dataframe1$related_exclusion)
head(dataframe1)
rm(related_exclusion)

######################################
## gestational diabetes - diabetes in pregnancy - pregnancy_diabetes
######################################
#open data
oa_r1b <- read.dta("./GWAS/ALSPAC/data/OA_r1b.dta")
pregnancy_diabetes <- oa_r1b[c("aln", "pregnancy_diabetes")]
pregnancy_diabetes$pregnancy_diabetes = as.character(pregnancy_diabetes$pregnancy_diabetes)
rm(oa_r1b)

#organise data
pregnancy_diabetes$qlet <- "M"
pregnancy_diabetes$aln_qlet <- paste0(pregnancy_diabetes$aln, "M")
pregnancy_diabetes <- pregnancy_diabetes[c("aln_qlet", "pregnancy_diabetes")]

## remove dupliates
inds = unique(pregnancy_diabetes[,1])
w = sapply(inds, function(x){   
  a = which(pregnancy_diabetes[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
pregnancy_diabetes = pregnancy_diabetes[w,]
rownames(pregnancy_diabetes) = pregnancy_diabetes[,1]

w = which(table(pregnancy_diabetes[,1]) > 1 )
## all the aln ids that are present more than once
n = names(table(pregnancy_diabetes[,1])[w])
##
w = which(pregnancy_diabetes[,1] %in% n)
pregnancy_diabetes[w,]
dim(pregnancy_diabetes)

#join
dataframe2 <- left_join(dataframe1, pregnancy_diabetes, by = "aln_qlet", all = T)
dim(dataframe2)

#check join
table(pregnancy_diabetes$pregnancy_diabetes)
table(dataframe2$pregnancy_diabetes)
head(dataframe2)

######################################
## glycosuria questionnaire
######################################

#open data
e_4d <- read.dta("./GWAS/ALSPAC/data/e_4d.dta")
glycosuria <- e_4d[c("aln", "e113")]
glycosuria$e113 = as.character(glycosuria$e113)
rm(e_4d)

#organise data
glycosuria$qlet <- "M"
glycosuria$aln_qlet <- paste0(glycosuria$aln, "M")
glycosuria <- glycosuria[c("aln_qlet", "e113")]
summary(glycosuria$e113)

## remove dupliates
inds = unique(glycosuria[,1])
w = sapply(inds, function(x){   
  a = which(glycosuria[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
glycosuria = glycosuria[w,]

rownames(glycosuria) = glycosuria[,1]

w = which(table(glycosuria[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(glycosuria[,1])[w])
##
w = which(glycosuria[,1] %in% n)
glycosuria[w,]

table(glycosuria$e113)

# check coding of pregnancy_diabetes and glycosuria variables in the different data varaibles from ALSPAC
glycosuria_diab <- left_join(glycosuria, pregnancy_diabetes, by = "aln_qlet", all = T)
dim(glycosuria_diab)
head(glycosuria_diab)
table(glycosuria_diab$e113, glycosuria_diab$pregnancy_diabetes)

#join
rownames(dataframe2) = dataframe2[,1]
dataframe3 <- left_join(dataframe2, glycosuria, by = "aln_qlet", all = T)
dim(dataframe3)
head(dataframe3)

#check join
table(dataframe3$e113)

#recode 'Missing' to 'NA'
dataframe3$e113[dataframe3$e113 == "Missing"] <- NA
#recode 'No' to 0
dataframe3$e113[dataframe3$e113 == "N"] <- 0
#recode 'Yes' to 1
dataframe3$e113[dataframe3$e113 == "Y"] <- 1

#check recode
table(dataframe3$e113)
table(dataframe3$e113, dataframe3$pregnancy_diabetes)
summary(dataframe3$e113)
head(dataframe3)

######################################
## stix
######################################
#open data
oa_r1b <- read.dta("./GWAS/ALSPAC/data/OA_r1b.dta")
stix <- oa_r1b[c("aln", "glycosuria")]
colnames(stix)[2] <- "stix_glycosuria"
stix$stix_glycosuria = as.character(stix$stix_glycosuria)
rm(oa_r1b)

#organise data
stix$qlet <- "M"
stix$aln_qlet <- paste0(stix$aln, "M")
stix <- stix[c("aln_qlet", "stix_glycosuria")]

## remove dupliates
inds = unique(stix[,1])
w = sapply(inds, function(x){   
  a = which(stix[,1] %in% x)[1]
  return(a)
})

### make a new dataframe
stix = stix[w,]
rownames(stix) = stix[,1]


w = which(table(stix[,1]) > 1 )
## all the aln ids that are present more than once
n = names(table(stix[,1])[w])
##
w = which(stix[,1] %in% n)
stix[w,]

table(stix$stix_glycosuria)

# check coding of pregnancy_diabetes and stix variables in the different data varaibles from ALSPAC
stix_diab <- left_join(stix, pregnancy_diabetes, by = "aln_qlet", all = T)
dim(stix_diab)
head(stix_diab)
table(stix_diab$stix_glycosuria, stix_diab$pregnancy_diabetes)

#join
rownames(dataframe3) = dataframe3[,1]
dataframe4 <- left_join(dataframe3, stix, by = "aln_qlet", all = T)

#check join
table(dataframe4$stix_glycosuria)
table(dataframe4$stix_glycosuria, dataframe4$pregnancy_diabetes)

#recode 'No' to 0
dataframe4$stix_glycosuria[dataframe4$stix_glycosuria == "No"] <- 0
#recode 'Yes' to 1
dataframe4$stix_glycosuria[dataframe4$stix_glycosuria == "Yes"] <- 1

#check recode
table(dataframe4$stix_glycosuria)
head(dataframe4)

######################################
## check and make final .txt file with/without exclusions
######################################

#check genetic and phenotype data
table(dataframe4$e113) #5910 1519 = 7429
1519/7429*100
table(dataframe4$pregnancy_diabetes) #7391 261 = 7652
261/7652*100
table(dataframe4$stix) #8022 323 = 8345
323/8345*100
dim(dataframe4)

write.table(dataframe4, "no_exclusions.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# remove individuals with existing diabetes and gestational diabetes
table(dataframe4$pregnancy_diabetes)
dataframe5 <- subset(dataframe4, pregnancy_diabetes != "Existing diabetes")
table(dataframe5$pregnancy_diabetes)
dataframe5 <- subset(dataframe5, pregnancy_diabetes != "Gestational diabetes")
table(dataframe5$pregnancy_diabetes)
table(dataframe5$e113) #5690 1399 = 7089
1399/7089*100
table(dataframe5$pregnancy_diabetes) #7391 261 = 7652
261/7652*100
table(dataframe5$stix) #7391 261 = 7652
261/7652*100
dim(dataframe5)

table(dataframe5$pregnancy_diabetes)
#recode 'yes' to 1
dataframe5$pregnancy_diabetes[dataframe5$pregnancy_diabetes == "Glycosuria"] <- 1
#recode 'no' to 0
dataframe5$pregnancy_diabetes[dataframe5$pregnancy_diabetes == "No glycosuria or diabetes"] <- 0
head(dataframe5)

# make exclusion dataframe
dataframe6 <- subset(dataframe5, WOC ==0)
dataframe6 <- subset(dataframe6, related_exclusion ==0)
table(dataframe6$e113) #5140 1249 = 6389
1249/6389*100
table(dataframe6$pregnancy_diabetes) #6639 227 = 6866
227/6866*100
table(dataframe6$stix) #6639 227 = 6866
227/6866*100

write.table(dataframe6, "yes_exclusions.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

######################################
## make  .txt file for all and for exclusions
######################################

all <- dataframe4[,1]
write.table(all, "all.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

exclude <- subset(dataframe4, WOC ==1)
exclude2 <- subset(dataframe4, related_exclusion ==1)
exclude3 <- subset(dataframe4, pregnancy_diabetes == "Existing diabetes")
exclude4 <- subset(dataframe4, pregnancy_diabetes == "Gestational diabetes")

exclude <- rbind(exclude, exclude2, exclude3, exclude4)
rm(exclude2, exclude3, exclude4)
exclude <- exclude[,1]

write.table(exclude, "exclude.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

include <- dataframe6[,1]
write.table(include, "include.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



######################################
## principle components for all
######################################

mother_PC <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/principal_components/mothers/data.eigenvec", header = F)
children_PC <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/principal_components/children/data.eigenvec", header = F)

PC <- rbind(mother_PC, children_PC)
colnames(PC)[1] <- "ID_1"
colnames(dataframe4)[1] <- "ID_1"

#join
rownames(PC) = PC[,1]
dataframe7 <- left_join(dataframe4, PC, by = "ID_1", all = T)

######################################
## finish phenofile for all
######################################
names(dataframe7)
pheno_file <- dataframe7
colnames(pheno_file) <- c("ID_1","ID_2","missing","father","mother","sex","plink_pheno","qlet","aln","WOC","related_exclusion","pregnancy_diabetes","e113","stix_glycosuria", "aln_qlet", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
head(pheno_file)
write.table(pheno_file, "all_pheno_file.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#pheno_file <- read.table("/Users/ml16847/OneDrive - University of Bristol/001_projects/001_Glycosuria/GWAS/ALSPAC_GWAS/pheno_file.txt", header = T)


######################################
## principle components for include list only
######################################

mother_PC <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/principal_components/mothers/data.eigenvec", header = F)
children_PC <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/principal_components/children/data.eigenvec", header = F)

PC <- rbind(mother_PC, children_PC)
colnames(PC)[1] <- "ID_1"

#join
rownames(PC) = PC[,1]
colnames(dataframe6)[1] <- "ID_1"
dataframe8 <- left_join(dataframe6, PC, by = "ID_1", all = T)

######################################
## finish phenofile for inlcude list
######################################

pheno_file <- dataframe8
colnames(pheno_file) <- c("ID_1","ID_2","missing","father","mother","sex","plink_pheno","qlet","aln","WOC","related_exclusion","pregnancy_diabetes","e113","stix_glycosuria", "aln_qlet", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
head(pheno_file)
names(pheno_file)

write.table(pheno_file, "pheno_file_include-list-only_regressions.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


pheno_file_new <- pheno_file[,c("ID_1","ID_2","missing","father","mother","sex","plink_pheno", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "e113")]
head(pheno_file_new)

write.table(pheno_file_new, "pheno_file_include-list-only.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# add in row at top of pheno_file_include-list-only.txt for PLINK



