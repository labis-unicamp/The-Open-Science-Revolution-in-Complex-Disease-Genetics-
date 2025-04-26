# Script para criar diagrama de Venn com variantes
# Carregando os pacotes necessários
library(VennDiagram)
library(gridExtra)
library(grid)

# Definindo o diretório de saída para os arquivos temporários do VennDiagram
dir.create("venn_temp", showWarnings = FALSE)
venn.dir <- "./venn_temp"

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

# Função para calcular as porcentagens
calculate_percentages <- function(counts, total) {
  percentages <- round((counts / total) * 100, 1)
  return(percentages)
}

# Calculando as porcentagens
variant_percentages <- calculate_percentages(variant_counts, total_unique_variants)
print("Porcentagem de variantes por ferramenta (em relação ao total único):")
print(variant_percentages)

# Configurando as cores para o diagrama de Venn
venn_colors <- c("red", "blue", "green", "purple")

# Calculando todas as interseções necessárias para o diagrama de Venn
A <- variant_sets[[1]]  # BCFTOOLS
B <- variant_sets[[2]]  # DeepVariant
C <- variant_sets[[3]]  # FreeBayes
D <- variant_sets[[4]]  # GATK

# Interseções de dois conjuntos
n12 <- length(intersect(A, B))
n13 <- length(intersect(A, C))
n14 <- length(intersect(A, D))
n23 <- length(intersect(B, C))
n24 <- length(intersect(B, D))
n34 <- length(intersect(C, D))

# Interseções de três conjuntos
n123 <- length(Reduce(intersect, list(A, B, C)))
n124 <- length(Reduce(intersect, list(A, B, D)))
n134 <- length(Reduce(intersect, list(A, C, D)))
n234 <- length(Reduce(intersect, list(B, C, D)))

# Interseção de todos os quatro conjuntos
n1234 <- length(Reduce(intersect, list(A, B, C, D)))

# Criando o diagrama de Venn manualmente com as interseções calculadas
futile.logger::flog.threshold(futile.logger::ERROR)

venn_plot <- draw.quad.venn(
  area1 = length(A),
  area2 = length(B),
  area3 = length(C),
  area4 = length(D),
  n12 = n12,
  n13 = n13,
  n14 = n14,
  n23 = n23,
  n24 = n24,
  n34 = n34,
  n123 = n123,
  n124 = n124,
  n134 = n134,
  n234 = n234,
  n1234 = n1234,
  category = paste0(
    names(variant_sets), 
    "\n(", variant_counts, ", ", variant_percentages, "%)"
  ),
  fill = venn_colors,
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 0.8,
  cat.fontface = "bold",
  cat.pos = c(0, 0, 180, 180),
  cat.dist = c(0.055, 0.055, 0.055, 0.055),
  cat.default.pos = "outer",
  margin = 0.1
)

# Salvando o diagrama de Venn em um arquivo PNG
png("variantes_venn_diagram.png", width = 800, height = 800, res = 100)
grid.draw(venn_plot)
dev.off()

# Exibindo o diagrama de Venn na tela
grid.newpage()
grid.draw(venn_plot)

# Calculando estatísticas adicionais para análise (opcional)
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

# Calculando as porcentagens de variantes exclusivas
exclusive_percentages <- calculate_percentages(exclusive_counts, total_unique_variants)
print("Porcentagem de variantes exclusivas por ferramenta:")
print(exclusive_percentages)

# Salvando as estatísticas em um arquivo de texto
sink("estatisticas_variantes.txt")
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

print("Análise completa! O diagrama de Venn foi salvo como 'variantes_venn_diagram.png'")
print("As estatísticas detalhadas foram salvas em 'estatisticas_variantes.txt'")
