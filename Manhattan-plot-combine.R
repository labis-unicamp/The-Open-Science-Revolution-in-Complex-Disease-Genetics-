# Manhattan plot combinado
Rscript -e '
  # Carregar pacote qqman
  library(qqman)

  # Carregar dados GWAS
  gwas_alzheimer <- read.table("gwas_ad_results.assoc", header=TRUE)
  gwas_diabetes_type2 <- read.table("gwas_dia_results.assoc", header=TRUE)

  # Carregar SNPs pleiotrópicos
  snps_pleiotropicos <- read.table("snps_pleiotropicos.txt", header=FALSE)$V1

  # Criar PDF com dois Manhattan plots
  pdf("manhattan_plots_combined.pdf", width=12, height=10)
  par(mfrow=c(2,1))

  # Manhattan plot para Alzheimer
  manhattan(gwas_alzheimer, highlight=snps_pleiotropicos, col=c("blue4", "skyblue"),
           main="Manhattan Plot - Alzheimer")

  # Manhattan plot para Diabetes tipo 2
  manhattan(gwas_diabetes_type2, highlight=snps_pleiotropicos, col=c("red4", "indianred1"),
           main="Manhattan Plot - Type 2 diabetes")

  dev.off()

  # Informar ao usuário
  cat("Plots gerados no arquivo manhattan_plots_combined.pdf\n")
'
