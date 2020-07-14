# step1 make phenofile

module add languages/R-3.5.1-ATLAS-gcc-6.1

R

library(foreign)
library(dplyr)
library(readstata13)

######################################
## genetic data
######################################

#aln file for all individuals
genetic_data <- read.table("data/data.sample", skip=2)
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
WOC <- read.csv("./001_glycosuria/GWAS/ALSPAC/data/ALSPAC_CONSENT_WITHDRAWN_LIST_PLUS_TripQuads.csv", header = T)
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
related_exclusion1 <- read.table("data/derived/unrelated_ids/mothers/exclusion_list.txt", header = F) #mothers
related_exclusion2 <- read.table("data/derived/unrelated_ids/children/exclusion_list.txt", header = F) #children
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
## hair colour
######################################

#open data
t_2a <- read.dta13("./001_glycosuria/reviewer_analysis/GWAS_revision2/data/t_2a.dta")
head(t_2a)
hair_colour <- t_2a[c("aln", "t6001")]
rm(t_2a)

#organise data
hair_colour$qlet <- "M"
hair_colour$aln_qlet <- paste0(hair_colour$aln, "M")
hair_colour <- hair_colour[c("aln_qlet", "t6001")]
summary(hair_colour$t6001)

## remove dupliates
inds = unique(hair_colour[,1])
w = sapply(inds, function(x){   
  a = which(hair_colour[,1] %in% x)[1]
  return(a)
})


### make a new dataframe
hair_colour = hair_colour[w,]

rownames(hair_colour) = hair_colour[,1]

w = which(table(hair_colour[,1]) > 1 )

## all the aln ids that are present more than once
n = names(table(hair_colour[,1])[w])
##
w = which(hair_colour[,1] %in% n)
hair_colour[w,]

# hair colour available on:
table(hair_colour$t6001)
hair_colour$t6001[hair_colour$t6001 == "-10"] <- NA
hair_colour$t6001[hair_colour$t6001 == "-1"] <- NA
#recode everything but blond hair to 0
hair_colour$t6001[hair_colour$t6001 == "2"] <- 0
hair_colour$t6001[hair_colour$t6001 == "3"] <- 0
hair_colour$t6001[hair_colour$t6001 == "4"] <- 0
hair_colour$t6001[hair_colour$t6001 == "5"] <- 0
table(hair_colour$t6001)
table(hair_colour$t6001)[2] / (table(hair_colour$t6001)[1] + table(hair_colour$t6001)[2]) * 100

#join
rownames(dataframe1) = dataframe1[,1]
dataframe2 <- left_join(dataframe1, hair_colour, by = "aln_qlet", all = T)

######################################
## check and make final .txt file with/without exclusions
######################################

#check genetic and phenotype data
table(dataframe2$t6001) 
table(dataframe2$t6001)[2] / (table(dataframe2$t6001)[1] + table(dataframe2$t6001)[2]) * 100

# save no exclusion data frame
write.table(dataframe2, "pheno_file/no_exclusions.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# make exclusion dataframe
dataframe3 <- subset(dataframe2, WOC ==0)
dataframe3 <- subset(dataframe3, related_exclusion ==0)
table(dataframe3$t6001) 
table(dataframe3$t6001)[2] / (table(dataframe3$t6001)[1] + table(dataframe3$t6001)[2]) * 100
write.table(dataframe3, "pheno_file/yes_exclusions.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

######################################
## make  aln_qlet .txt file for all and for exclusions
######################################
all <- dataframe2[,1]
write.table(all, "pheno_file/all.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# exclusions
exclude <- subset(dataframe2, WOC ==1)
exclude2 <- subset(dataframe2, related_exclusion ==1)
exclude <- rbind(exclude, exclude2)
rm(exclude2)
exclude <- exclude[,1]
# exclude thes epeople from GWAS
write.table(exclude, "pheno_file/exclude.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# include the people in the gWAS
include <- dataframe3[,1]
write.table(include, "pheno_file/include.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

######################################
## principle components for all
######################################

mother_PC <- read.table("data/derived/principal_components/mothers/data.eigenvec", header = F)
children_PC <- read.table("data/derived/principal_components/children/data.eigenvec", header = F)

PC <- rbind(mother_PC, children_PC)
colnames(PC)[1] <- "ID_1"
colnames(dataframe2)[1] <- "ID_1"

#join
rownames(PC) = PC[,1]
dataframe4 <- left_join(dataframe2, PC, by = "ID_1", all = T)

######################################
## finish phenofile for all
######################################
names(dataframe4)
pheno_file <- dataframe4
colnames(pheno_file) <- c("ID_1","ID_2","missing","father","mother","sex","plink_pheno","qlet","aln","WOC","related_exclusion","t6001", "aln_qlet", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
head(pheno_file)
write.table(pheno_file, "pheno_file/all_pheno_file.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#pheno_file <- read.table("/Users/ml16847/OneDrive - University of Bristol/001_projects/001_Glycosuria/GWAS/ALSPAC_GWAS/pheno_file.txt", header = T)


######################################
## principle components for include list only
######################################

mother_PC <- read.table("data/derived/principal_components/mothers/data.eigenvec", header = F)
children_PC <- read.table("data/derived/principal_components/children/data.eigenvec", header = F)

PC <- rbind(mother_PC, children_PC)
colnames(PC)[1] <- "ID_1"

#join
rownames(PC) = PC[,1]
colnames(dataframe3)[1] <- "ID_1"
dataframe5 <- left_join(dataframe3, PC, by = "ID_1", all = T)

######################################
## finish phenofile for inlcude list
######################################

pheno_file <- dataframe5
colnames(pheno_file) <- c("ID_1","ID_2","missing","father","mother","sex","plink_pheno","qlet","aln","WOC","related_exclusion","t6001", "aln_qlet", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
head(pheno_file)
names(pheno_file)

write.table(pheno_file, "pheno_file/pheno_file_include-list-only_regressions.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

pheno_file_new <- pheno_file[,c("ID_1","ID_2","missing","father","mother","sex","plink_pheno", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "t6001")]
head(pheno_file_new)

write.table(pheno_file_new, "pheno_file/pheno_file_include-list-only.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# add in row at top of pheno_file_include-list-only.txt for PLINK
0       0       0       D       D       D       D       C       C       C       C       C       C       C       C       C       C       C       C       C       C       C       C       C       C       C       C       B


