# Manhattan plot combinado
Rscript -e '
  # Load qqman package
  library(qqman)

  # Load GWAS data
  gwas_alzheimer <- read.table("gwas_ad_results.assoc", header=TRUE)
  gwas_diabetes_type2 <- read.table("gwas_dia_results.assoc", header=TRUE)

  # Load pleiotropic SNPs
  snps_pleiotropicos <- read.table("snps_pleiotropicos.txt", header=FALSE)$V1

  # Create PDF with two Manhattan plots
  pdf("manhattan_plots_combined.pdf", width=12, height=10)
  par(mfrow=c(2,1))

  # Manhattan plot for Alzheimer
  manhattan(gwas_alzheimer, highlight=snps_pleiotropicos, col=c("blue4", "skyblue"),
           main="Manhattan Plot - Alzheimer")

  # Manhattan plot for Type 2 diabetes
  manhattan(gwas_diabetes_type2, highlight=snps_pleiotropicos, col=c("red4", "indianred1"),
           main="Manhattan Plot - Type 2 diabetes")

  dev.off()

  # Notify user
  cat("Plots saved to manhattan_plots_combined.pdf\n")
'
