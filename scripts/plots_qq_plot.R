library(ggplot2)
library(data.table)
data <- fread("./GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final.txt", header = T)

gg.qqplot <- function(ps, title, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), alpha = 1, size = 0.5, colour = "skyblue") +
    geom_abline(intercept = 0, slope = 1, alpha = 1, colour = "grey", size = 0.5) +
    geom_line(aes(expected, cupper), linetype = "solid", size = 1, colour = "black") +
    geom_line(aes(expected, clower), linetype = "solid", size = 1, colour = "black") +
    xlab(log10Pe) +
    ylab(log10Po) +
    
    # add plot title
    ggtitle(paste0(title)) +
    
    
    # Customise the theme:
    theme_bw(base_size = 11) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank() 
    )
}


png("./GWAS/ALSPAC/plots/Figure 2. QQ plot of GWAS of reported glycosuria in the third trimester in ALSPAC.png",
    height = 100,
    width = 100,
    units = "mm",
    res = 1000)
gg.qqplot(data$P,
          title = "λ = 1.007438",
          ci=NA) #remove confidence intervals with 'NA'
dev.off()
