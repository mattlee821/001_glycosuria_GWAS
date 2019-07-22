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

module add languages/R-3.5-ATLAS-gcc-7.1.0

R

library(foreign)
library(dplyr)

######################################
## genetic data
######################################

#aln file for all individuals
dataframe <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/data.sample", skip=2)
colnames(dataframe) <- c("ID_1","ID_2","missing","father","mother","sex","plink_pheno")
#dataframe <- dataframe[1]
colnames(dataframe)[1] <- "aln_qlet"

#split dataframe aln and qlet column into two
dataframe$qlet <- gsub("[0-9]", "", dataframe$aln_qlet)
dataframe$aln <- gsub("[a-zA-Z]", "", dataframe$aln_qlet)
dataframe$qlet <- as.factor(dataframe$qlet)
table(dataframe$qlet)


######################################
## withdrawl of consent
######################################

#withdrawl of consent file
WOC <- read.csv("./GWAS/ALSPAC/data/ALSPAC_CONSENT_WITHDRAWN_LIST_PLUS_TripQuads.csv", header = T)
WOC <- WOC[4]
WOC$WOC <- 1

#create new dataframe with aln and WOC. where WOC individuals have 1 and everyone else has 0
rownames(dataframe) = dataframe[,1]
dataframe2 <- left_join(dataframe, WOC, by="aln_qlet")
dataframe2$WOC <- as.factor(dataframe2$WOC)
summary(dataframe2$WOC)
dataframe2$WOC <- as.character(as.numeric(dataframe2$WOC))
dataframe2$WOC[is.na(dataframe2$WOC)] <- 0
summary(dataframe2$WOC)
rm(WOC)

table(dataframe2$qlet, dataframe2$WOC)

######################################
## glycosuria questionnaire
######################################

#open data
e_4d <- read.dta("./GWAS/ALSPAC/data/e_4d.dta")
head(e_4d)
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

#summary(stix$stix_glycosuria)
table(glycosuria$e113)

#join
rownames(dataframe2) = dataframe2[,1]
dataframe3 <- left_join(dataframe2, glycosuria, by = "aln_qlet", all = T)

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
summary(dataframe3$e113)


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

#summary(stix$stix_glycosuria)
table(stix$stix_glycosuria)

#join
rownames(dataframe3) = dataframe3[,1]
dataframe4 <- left_join(dataframe3, stix, by = "aln_qlet", all = T)

#check join
table(dataframe4$stix)

#recode 'No' to 0
dataframe4$stix[dataframe4$stix == "No"] <- 0
#recode 'Yes' to 1
dataframe4$stix[dataframe4$stix == "Yes"] <- 1

#check recode
table(dataframe4$stix)


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
dataframe5 <- left_join(dataframe4, related_exclusion, by = "aln_qlet", all = T)
dataframe5$related_exclusion[is.na(dataframe5$related_exclusion)] <- 0

#check join
table(dataframe5$related_exclusion)
rm(related_exclusion)


######################################
## check and make final .txt file with/without exclusions
######################################

#check genetic and phenotype data
table(dataframe5$e113) #5910 1519 = 7429
table(dataframe5$stix) #8022 323 = 8345

write.table(dataframe5, "no_exclusions.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#make exclusion dataframe
dataframe6 <- subset(dataframe5, WOC ==0)
dataframe6 <- subset(dataframe6, related_exclusion ==0)

#check genetic and phenotype data
table(dataframe6$e113) #5336 1354 = 6690
table(dataframe6$stix) #7180 283 = 7463

write.table(dataframe6, "yes_exclusions.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

######################################
## make  .txt file for all and for exclusions
######################################

all <- dataframe5[,1]
exclude <- subset(dataframe5, WOC ==1)
exclude2 <- subset(dataframe5, related_exclusion ==1)
exclude <- rbind(exclude, exclude2)
rm(exclude2)
exclude <- exclude[,1]

write.table(all, "all.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(exclude, "exclude.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

######################################
## principle components for all
######################################

mother_PC <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/principal_components/mothers/data.eigenvec", header = F)
children_PC <- read.table("./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/derived/principal_components/children/data.eigenvec", header = F)

PC <- rbind(mother_PC, children_PC)
colnames(PC)[1] <- "ID_1"
colnames(dataframe5)[1] <- "ID_1"

#join
rownames(PC) = PC[,1]
dataframe7 <- left_join(dataframe5, PC, by = "ID_1", all = T)

######################################
## finish phenofile for all
######################################

pheno_file <- dataframe7[,c(1,2,3,4,5,6,7,11,13,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
colnames(pheno_file) <- c("ID_1","ID_2","missing","father","mother","sex","plink_pheno", "glycosuria_questionnaire_e113", "glycosuria_stix_clinic", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

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

pheno_file <- dataframe8[,c(1,2,3,4,5,6,7,11,13,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
colnames(pheno_file) <- c("ID_1","ID_2","missing","father","mother","sex","plink_pheno", "glycosuria_questionnaire_e113", "glycosuria_stix_clinic", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
head(pheno_file)

pheno_file_new <- pheno_file[,c("ID_1","ID_2","missing","father","mother","sex","plink_pheno", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "glycosuria_questionnaire_e113", "glycosuria_stix_clinic")]
head(pheno_file_new)

write.table(pheno_file_new, "pheno_file_include-list-only.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")




