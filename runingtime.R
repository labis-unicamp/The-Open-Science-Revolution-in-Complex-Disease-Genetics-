# Carregar pacotes necessários
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Criar dataframe com os dados
data <- data.frame(
  Sample = c("SRR1279767", "SRR1279773", "SRR1279777", "average"),
  GATK = c(52.875, 39.545, 55.516, 49.3119),
  freebayes = c(25.398, 20.456, 25.936, 25.5965),
  DeepVariant = c(6.821, 4.995, 5.555, 5.7903),
  bcftools = c(20.160, 14.997, 17.472, 17.5431)
)

# Converter para formato longo (tidy)
data_long <- data %>%
  pivot_longer(cols = -Sample, names_to = "Software", values_to = "Minutes")

# Ordenar os softwares e amostras
data_long$Software <- factor(data_long$Software, 
                            levels = c("GATK", "freebayes", "bcftools", "DeepVariant"))
data_long$Sample <- factor(data_long$Sample, 
                          levels = c("SRR1279767", "SRR1279773", "SRR1279777", "average"))

# Criar o gráfico
ggplot(data_long, aes(x = Software, y = Minutes, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", Minutes)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")) +
  labs(title = "Runtime comparison between software",
       x = "Software",
       y = "Execution time (minutes)",
       fill = "Sample") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Salvar o gráfico
ggsave("runtime_comparison.png", width = 8, height = 6, dpi = 300)
