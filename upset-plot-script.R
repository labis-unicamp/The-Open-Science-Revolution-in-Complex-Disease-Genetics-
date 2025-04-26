# Script para criar um UpSet plot com dados de variantes

# Primeiro, verificar e instalar os pacotes necessários se não estiverem presentes
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

# Carregando os pacotes necessários
library(UpSetR)
library(dplyr)

# Função para ler as variantes de um arquivo
read_variants <- function(file_path) {
  variants <- readLines(file_path)
  return(variants)
}

# Lendo os arquivos de variantes
bcftools_variants <- read_variants("var_bcftools_all_variants.txt")
deepvariants_variants <- read_variants("var_deepvariants_all_variants.txt")
freebayes_variants <- read_variants("var_freebayes_all_variants.txt")
gatk_variants <- read_variants("var_gatk_all_variants.txt")

# Criando uma lista com os conjuntos de variantes
variant_sets <- list(
  BCFTOOLS = bcftools_variants,
  DeepVariant = deepvariants_variants,
  FreeBayes = freebayes_variants,
  GATK = gatk_variants
)

# Calculando o número total de variantes em cada conjunto
variant_counts <- sapply(variant_sets, length)
print("Número de variantes por ferramenta:")
print(variant_counts)

# Calculando o número total de variantes únicas em todos os conjuntos
total_unique_variants <- length(unique(unlist(variant_sets)))
print(paste("Total de variantes únicas:", total_unique_variants))

# Criando uma tabela de presença/ausência para cada variante em cada conjunto
all_variants <- unique(unlist(variant_sets))
presence_matrix <- matrix(0, nrow = length(all_variants), ncol = length(variant_sets))
colnames(presence_matrix) <- names(variant_sets)
rownames(presence_matrix) <- all_variants

for (i in 1:length(variant_sets)) {
  set_name <- names(variant_sets)[i]
  presence_matrix[variant_sets[[i]], set_name] <- 1
}

# Convertendo a matriz para um dataframe
variant_df <- as.data.frame(presence_matrix)

# Calculando estatísticas
# Interseções entre todas as ferramentas
all_intersect <- Reduce(intersect, variant_sets)
print(paste("Variantes detectadas por todas as ferramentas:", length(all_intersect)))

# Variantes exclusivas de cada ferramenta
exclusive_variants <- list()
for (name in names(variant_sets)) {
  other_sets <- variant_sets[names(variant_sets) != name]
  other_variants <- unique(unlist(other_sets))
  exclusive_variants[[name]] <- setdiff(variant_sets[[name]], other_variants)
}

# Imprimindo o número de variantes exclusivas por ferramenta
exclusive_counts <- sapply(exclusive_variants, length)
print("Número de variantes exclusivas por ferramenta:")
print(exclusive_counts)

# Calculando as porcentagens relativas ao total
variant_percentages <- round((variant_counts / total_unique_variants) * 100, 1)
print("Porcentagem de variantes por ferramenta (em relação ao total único):")
print(variant_percentages)

exclusive_percentages <- round((exclusive_counts / total_unique_variants) * 100, 1)
print("Porcentagem de variantes exclusivas por ferramenta:")
print(exclusive_percentages)

# Salvando as estatísticas em um arquivo de texto
sink("estatisticas_variantes_upset.txt")
cat("Estatísticas de Variantes por Ferramenta\n")
cat("========================================\n\n")
cat("Número total de variantes por ferramenta:\n")
print(variant_counts)
cat("\nPorcentagem de variantes por ferramenta:\n")
print(variant_percentages)
cat("\nNúmero de variantes exclusivas por ferramenta:\n")
print(exclusive_counts)
cat("\nPorcentagem de variantes exclusivas por ferramenta:\n")
print(exclusive_percentages)
cat("\nVariantes detectadas por todas as ferramentas:", length(all_intersect), "\n")
cat("\nTotal de variantes únicas:", total_unique_variants, "\n")
sink()

# Criando o UpSet plot básico
png("variantes_upset_plot.png", width = 1000, height = 800, res = 120)
upset_plot <- upset(
  variant_df,
  sets = colnames(presence_matrix),
  order.by = "freq",
  mainbar.y.label = "Tamanho da Interseção",
  sets.x.label = "Tamanho do Conjunto",
  text.scale = c(1.3, 1.3, 1, 1, 1.3, 1),
  point.size = 3,
  line.size = 1
)
dev.off()

# Criando o UpSet plot com mais detalhes
png("variantes_upset_plot_detalhado.png", width = 1200, height = 800, res = 120)
upset_plot_detailed <- upset(
  variant_df,
  sets = colnames(presence_matrix),
  order.by = "freq",
  mb.ratio = c(0.55, 0.45),
  number.angles = 0,
  text.scale = c(1.3, 1.3, 1, 1, 1.3, 1),
  point.size = 3.5,
  line.size = 1,
  mainbar.y.label = "Número de Variantes Compartilhadas",
  sets.x.label = "Variantes por Ferramenta",
  set_size.show = TRUE,
  set_size.numbers_size = 8,
  set_size.scale_max = max(variant_counts) * 1.1,
  sets.bar.color = c("red", "blue", "green", "purple")
)
dev.off()

# Criando um gráfico adicional com anotações
png("variantes_upset_plot_anotado.png", width = 1200, height = 800, res = 120)

# Definindo as cores para as ferramentas
set_colors <- c("red", "blue", "green", "purple")
names(set_colors) <- names(variant_sets)

# Criando o terceiro plot com anotações mais simples
upset_plot_anotado <- upset(
  variant_df,
  sets = colnames(presence_matrix),
  order.by = "freq",
  mb.ratio = c(0.55, 0.45),
  number.angles = 0,
  text.scale = c(1.3, 1.3, 1, 1, 1.3, 1),
  point.size = 3.5,
  line.size = 1,
  mainbar.y.label = "Número de Variantes Compartilhadas",
  sets.x.label = "Variantes por Ferramenta",
  set_size.show = TRUE,
  set_size.scale_max = max(variant_counts) * 1.1,
  sets.bar.color = set_colors,
  keep.order = TRUE
)
dev.off()

print("Análise completa! Os UpSet plots foram salvos como:")
print("1. variantes_upset_plot.png (básico)")
print("2. variantes_upset_plot_detalhado.png (detalhado)")
print("3. variantes_upset_plot_anotado.png (com anotações)")
print("As estatísticas detalhadas foram salvas em 'estatisticas_variantes_upset.txt'")
