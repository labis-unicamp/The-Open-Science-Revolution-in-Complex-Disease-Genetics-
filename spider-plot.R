#!/usr/bin/env Rscript

# Verificar e instalar pacotes se necessário
required_packages <- c("fmsb", "dplyr", "tidyr", "tibble")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

library(fmsb)
library(dplyr)
library(tidyr)
library(tibble)

# Dados fornecidos (convertendo "XmYs" para minutos decimais)
dados <- list(
  "bwa_mem" = data.frame(
    row.names = c("real", "user", "sys"),
    SRR1279767 = c(139 + 35.495/60, 129 + 28.981/60, 19 + 18.789/60),
    SRR1279773 = c(96 + 5.992/60, 89 + 12.153/60, 13 + 11.416/60),
    SRR1279777 = c(127 + 47.054/60, 117 + 37.493/60, 19 + 25.048/60)
  ),
  "bwa_mem2" = data.frame(
    row.names = c("real", "user", "sys"),
    SRR1279767 = c(87 + 45.322/60, 83 + 59.037/60, 12 + 9.096/60),
    SRR1279773 = c(62 + 8.247/60, 60 + 27.872/60, 8 + 30.130/60),
    SRR1279777 = c(75 + 38.656/60, 75 + 55.222/60, 9 + 31.863/60)
  ),
  "bowtie" = data.frame(
    row.names = c("real", "user", "sys"),
    SRR1279767 = c(194 + 58.700/60, 392 + 1.445/60, 12 + 19.679/60),
    SRR1279773 = c(129 + 31.212/60, 260 + 59.470/60, 8 + 4.367/60),
    SRR1279777 = c(183 + 30.620/60, 370 + 54.058/60, 10 + 11.078/60)
  )
)

# Converter para horas (opcional)
dados <- lapply(dados, function(x) x / 60)

# --- Gráfico 1: Médias por alinhador (user, real, sys) ---
# Calcular médias e converter para data.frame sem nomes de linha
medias <- lapply(dados, function(df) rowMeans(df)) %>% 
  bind_rows(.id = "Alinhador") %>% 
  as.data.frame()

# Definir nomes de linha a partir da coluna "Alinhador"
rownames(medias) <- medias$Alinhador
medias <- medias[, -1, drop = FALSE]

# Preparar dados para o radar
radar_data <- rbind(
  rep(ceiling(max(medias) * 1.1), ncol(medias)),  # Máximo
  rep(0, ncol(medias)),                           # Mínimo
  medias
)

# Plotar
png("grafico1_medias.png", width = 800, height = 600)
radarchart(
  radar_data,
  title = "Average Time per Aligner (hours)",
  pcol = c("#1B9E77", "#D95F02", "#7570B3"),
  plwd = 2,
  pfcol = scales::alpha(c("#1B9E77", "#D95F02", "#7570B3"), 0.3),
  cglcol = "grey",
  vlcex = 0.9
)
legend("topright", legend = rownames(medias), fill = c("#1B9E77", "#D95F02", "#7570B3"))
dev.off()

# --- Gráfico 2: Tempo real por amostra ---
# Extrair dados de tempo real e converter para data.frame sem nomes de linha
real_data <- lapply(dados, function(df) df["real", ]) %>% 
  bind_rows(.id = "Alinhador") %>% 
  as.data.frame()

# Definir nomes de linha a partir da coluna "Alinhador"
rownames(real_data) <- real_data$Alinhador
real_data <- real_data[, -1, drop = FALSE]

# Preparar dados para o radar
radar_data_real <- rbind(
  rep(ceiling(max(real_data) * 1.1), ncol(real_data)),
  rep(0, ncol(real_data)),
  real_data
)

# Plotar
png("grafico2_tempo_real.png", width = 800, height = 600)
radarchart(
  radar_data_real,
  title = "Real Time Per Sample (hours)",
  pcol = c("#1B9E77", "#D95F02", "#7570B3"),
  plwd = 2,
  pfcol = scales::alpha(c("#1B9E77", "#D95F02", "#7570B3"), 0.3),
  cglcol = "grey",
  vlcex = 0.9
)
legend("topright", legend = rownames(real_data), fill = c("#1B9E77", "#D95F02", "#7570B3"))
dev.off()
