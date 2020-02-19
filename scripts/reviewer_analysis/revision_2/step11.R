data <- read.table("reviewer_analysis/GWAS_revision2/mother_GWAS/sig_snps_mum-offspring.txt", header = T, sep = "\t")
data <- data[,c("SNP","chr","EA","NEA","P","BETA","SE","group")]

# make data frame including doubled child estimates
mum <- subset(data, group == "mum")
child <- subset(data, group == "offspring")
child_doubled <- child
child_doubled$BETA <- child_doubled$BETA * 2
child_doubled$SE <- child_doubled$SE * 2
child_doubled$group <- "offspring_doubled" 
data <- rbind(data, child_doubled)

# calculate ratio of raw child beta to raw mum beta
colnames(mum) <- paste(colnames(mum), sep = "", "_mum")
colnames(child) <- paste(colnames(child), sep = "", "_child")
colnames(child_doubled) <- paste(colnames(child_doubled), sep = "", "_child2")

colnames(mum)[1] <- "SNP"
colnames(child)[1] <- "SNP"
colnames(child_doubled)[1] <- "SNP"

library(dplyr)
data1 <- left_join(mum, child, by = "SNP")
data1 <- left_join(data1, child_doubled)

# overall
mean(data1$BETA_child / data1$BETA_mum)
mean(data1$BETA_child2 / data1$BETA_mum)
png("reviewer_analysis/GWAS_revision2/histogram_all_ratio.png")
hist(data1$BETA_child / data1$BETA_mum, main = "All chromosomes, 124 SNPs", breaks = nrow(data1))
dev.off()

# chromosome 6
chr6 <- subset(data1, chr_mum == 6)
mean(chr6$BETA_child / chr6$BETA_mum)
mean(chr6$BETA_child2 / chr6$BETA_mum)
png("reviewer_analysis/GWAS_revision2/histogram_chr6_ratio.png")
hist(chr6$BETA_child / chr6$BETA_mum, main = "chromsome 6, 3 SNPs", breaks = nrow(chr6))
dev.off()

# chromosome 14
chr14 <- subset(data1, chr_mum == 14)
mean(chr14$BETA_child / chr14$BETA_mum)
mean(chr14$BETA_child2 / chr14$BETA_mum)
png("reviewer_analysis/GWAS_revision2/histogram_chr14_ratio.png")
hist(chr14$BETA_child / chr14$BETA_mum, main = "chromsome 14, 15 SNPs", breaks = nrow(chr14))
dev.off()

# chromosome 15
chr15 <- subset(data1, chr_mum == 15)
mean(chr15$BETA_child / chr15$BETA_mum)
mean(chr15$BETA_child2 / chr15$BETA_mum)
png("reviewer_analysis/GWAS_revision2/histogram_chr15_ratio.png")
hist(chr15$BETA_child / chr15$BETA_mum, main = "chromsome 15, 106 SNPs", breaks = nrow(chr15))
dev.off()

# glycosuria histogram
glycosuria_child <- read.table("reviewer_analysis/GWAS_revision2/mother_GWAS/glm_logit_glycosuria_child-dosage.txt", header = T, sep = "\t")
colnames(glycosuria_child)[1] <- "SNP"
glycosuria_mums <- read.table("reviewer_analysis/GWAS_revision2/mother_GWAS/GWAS_output_final_sig.txt", header = T, sep = "\t")
glycosuria_mums <- subset(glycosuria_mums, chr == 16)
glycosuria <- left_join(glycosuria_mums, glycosuria_child, by = "SNP")

mean(glycosuria$b / glycosuria$BETA)
mean(glycosuria$b*2 / glycosuria$BETA)
png("reviewer_analysis/GWAS_revision2/histogram_glycosuria_ratio.png")
hist(glycosuria$b / glycosuria$BETA, main = "Glycosuria: chromsome 16, 54 SNPs", breaks = nrow(glycosuria))
dev.off()


# correlation between raw child beta to raw mum beta
cor(data1$BETA_child, data1$BETA_mum)
cor(data1$BETA_child2, data1$BETA_mum)

# forest plot of all SNP estimates ====
##### Packages
library(ggplot2)
library(ggforestplot)

##### Plotting variables
ci <- 0.95
psignif <- 5e-8

data <- data[order(data$chr, data$P),]
data$chr <- as.factor(data$chr)
data <- subset(data, SNP != "rs2346098")

##### Plot
pdf("reviewer_analysis/GWAS_revision2/mother_GWAS/SNP_forestplot.pdf", height = 20)
forestplot(df = data, 
           name = SNP, 
           estimate = BETA, 
           se = SE, 
           pvalue = P, 
           colour = group, 
           shape = chr, 
           logodds = FALSE, 
           psignif = psignif, 
           ci = ci, 
           title = "Forestplot",
           xlab = "Forestplot with 95% CI")
dev.off()

## Plotting variables
ci <- 0.95
psignif <- 5e-8

## plot
xmin <- min(data$BETA * 1.96 - data$SE)
xmax <- max(data$BETA * 1.96 + data$SE)

library(tidyr)
library(purrr)
data_sig_nested <- 
  data %>%
  group_by(chr) %>%
  nest()

svg("reviewer_analysis/GWAS_revision2/mother_GWAS/SNP_forestplot_groups.pdf", height = 30)
forestplot_multi <-
  data_sig_nested %>%
  mutate(gg_groups = map2(data, chr, ~ forestplot(df = .x, 
                                                       name = SNP, 
                                                       estimate = BETA, 
                                                       se = SE, 
                                                       pvalue = P, 
                                                       colour = group, 
                                                       shape = NULL, 
                                                       logodds = FALSE, 
                                                       psignif = psignif,
                                                       ci = ci,
                                                       title = .y,
                                                       xlab = "Effect estimate with 95% CI") +
                            theme(
                              legend.position = "bottom",
                              legend.title.align = 0,
                              legend.text.align = 0,
                              legend.title = element_blank(),
                              axis.title.x = element_text(hjust = 0.5)) +
                            ggplot2::coord_cartesian(xlim = c(-3.5, 3))
  ),
  # Optional: remove x-axis and legend for all plots except the bottom one
  gg_groups = ifelse(
    test = row_number() != n(),
    yes =
      purrr::map(gg_groups, ~ . +
                   theme(
                     axis.text.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     plot.margin = unit(c(1, 2, 1, 2), "mm")
                   ) +
                   ggplot2::theme(legend.position = "none")),
    no = gg_groups
  ),
  rel_heights = purrr::map(
    data,  ~ nrow(.) 
  ) %>% unlist()
  )

patchwork::wrap_plots(
  forestplot_multi$gg_groups,
  ncol = 1, 
  heights = forestplot_multi$rel_heights
)
dev.off()




