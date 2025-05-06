# Scatter plot of p-values:
Rscript -e '
  data <- read.table("pleiotropy_results.txt", header=TRUE)
  png("pvalue_scatter.png", width=800, height=800, res=100)
  plot(-log10(data$P1), -log10(data$P2),
       xlab="-log10(P-valor Alzheimer)", ylab="-log10(P-valor Type2 Diabetes)",
       main="Significance Comparison between Phenotypes",
       pch=19, col=ifelse(data$P_pleio < 0.05, "red", "black"))
  abline(h=-log10(0.05), v=-log10(0.05), lty=2, col="blue")
  legend("topright", legend=c("Pleiotropy (p<0.05)", "Not significant"),
         col=c("red", "black"), pch=19)
  dev.off()
'
