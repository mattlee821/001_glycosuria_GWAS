module add languages/R-3.5-ATLAS-gcc-7.1.0
R

data <- read.table("./other_analysis/regressions/pheno_file.txt", header = T)

head(data)

###############################################
## startify by gestational hypertension
summary(data$gesthyp)
table(data$gesthyp)
gesthyp_yes <- subset(data, gesthyp == 1)
gesthyp_no <- subset(data, gesthyp == 0)
table(gesthyp_yes$gesthyp)
table(gesthyp_no$gesthyp)

table(gesthyp_yes$e113)
table(gesthyp_no$e113)

## startify by hypertension pregnancy 
summary(data$preeclampsia)
table(data$preeclampsia)
preeclampsia_yes <- subset(data, preeclampsia == 1)
preeclampsia_no <- subset(data, preeclampsia == 0)
table(preeclampsia_yes$preeclampsia)
table(preeclampsia_no$preeclampsia)

table(preeclampsia_yes$e113)
table(preeclampsia_no$e113)

###

#gesthyp YES
table(gesthyp_yes$e113)
nrow(gesthyp_yes)
e113_gesthyp <- glm(gesthyp_yes$e113 ~ gesthyp_yes[,"rs13337037_dosage"], family=binomial(link="logit"))
e113_gesthyp <- summary(e113_gesthyp)$coefficients
e113_gesthyp
summary(e113_gesthyp)
#ci
0.178257 - 1.96 * 0.1131731 # -0.04356228
0.178257 + 1.96 * 0.1131731 # 0.4000763
exp(0.178257)
exp(-0.04356228)
exp(0.4000763)
### 

#gesthyp NO
table(gesthyp_no$e113)
nrow(gesthyp_no)
e113_gesthyp <- glm(gesthyp_no$e113 ~ gesthyp_no[,"rs13337037_dosage"], family=binomial(link="logit"))
e113_gesthyp <- summary(e113_gesthyp)$coefficients
e113_gesthyp
summary(e113_gesthyp)
#ci
0.3885649 - 1.96 * 0.05307484 # 0.2845382
0.3885649 + 1.96 * 0.05307484 # 0.4925916
exp(0.3885649)
exp(0.2845382)
exp(0.4925916)

######### gesthyp_interaction
gesthyp_interaction <- glm(data$e113 ~ data$rs13337037_dosage * data$gesthyp, family=binomial(link="logit"))
gesthyp_interaction <- summary(gesthyp_interaction)$coefficients
gesthyp_interaction

#SNP
0.3885649 - 1.96 * 0.05307484 # 
0.3885649 + 1.96 * 0.05307484 # 
exp(0.3885649)
exp(0.2845382)
exp(0.4925916)

#gesthyp
0.3885649 - 1.96 * 0.05307484 # 
0.3885649 + 1.96 * 0.05307484 # 
exp(0.3885649)
exp(0.2845382)
exp(0.4925916)

#interaction
0.3885649 - 1.96 * 0.05307484 # 
0.3885649 + 1.96 * 0.05307484 # 
exp(0.3885649)
exp(0.2845382)
exp(0.4925916)

#preeclampsia YES
table(preeclampsia_yes$e113)
nrow(preeclampsia_yes)
e113_preeclampsia <- glm(preeclampsia_yes$e113 ~ preeclampsia_yes[,"rs13337037_dosage"], family=binomial(link="logit"))
e113_preeclampsia <- summary(e113_preeclampsia)$coefficients
e113_preeclampsia
summary(e113_preeclampsia)
#ci
0.7339117 - 1.96 * 0.3000565 # 0.145801
0.7339117 + 1.96 * 0.3000565 # 1.322022
exp(0.7339117)
exp(0.145801)
exp(1.322022)
### 

#preeclampsia NO
table(preeclampsia_no$e113)
nrow(preeclampsia_no)
e113_preeclampsia <- glm(preeclampsia_no$e113 ~ preeclampsia_no[,"rs13337037_dosage"], family=binomial(link="logit"))
e113_preeclampsia <- summary(e113_preeclampsia)$coefficients
e113_preeclampsia
summary(e113_preeclampsia)
#ci
0.3457821 - 1.96 * 0.04789057 # 0.2519166
0.3457821 + 1.96 * 0.04789057 # 0.4396476
exp(0.3457821)
exp(0.2519166)
exp(0.4396476)


######### preeclampsia_interaction
preeclampsia_interaction <- glm(data = data, e113 ~ rs13337037_dosage * preeclampsia, family=binomial(link="logit"))
preeclampsia_interaction <- summary(preeclampsia_interaction)$coefficients
preeclampsia_interaction

#SNP
0.3885649 - 1.96 * 0.05307484 # 
0.3885649 + 1.96 * 0.05307484 # 
exp(0.3885649)
exp(0.2845382)
exp(0.4925916)

#gesthyp
0.3885649 - 1.96 * 0.05307484 # 
0.3885649 + 1.96 * 0.05307484 # 
exp(0.3885649)
exp(0.2845382)
exp(0.4925916)

#interaction
0.3885649 - 1.96 * 0.05307484 # 
0.3885649 + 1.96 * 0.05307484 # 
exp(0.3885649)
exp(0.2845382)