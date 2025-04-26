# Script R para gerar o gráfico UpSet

# Instalar e carregar pacotes necessários
if (!require("UpSetR")) install.packages("UpSetR")
if (!require("tidyverse")) install.packages("tidyverse")
library(UpSetR)
library(tidyverse)

# Função para ler e processar cada arquivo
read_variants <- function(file_path, tool_name) {
  read_tsv(file_path, col_names = "variant", show_col_types = FALSE) %>% 
    distinct(variant) %>% 
    mutate(tool = tool_name)
}

# Carregar e combinar todos os dados
all_data <- bind_rows(
  read_variants("var_deepvariants_all_variants.txt", "DeepVariant"),
  read_variants("var_gatk_all_variants.txt", "GATK"),
  read_variants("var_bcftools_all_variants.txt", "bcftools"),
  read_variants("var_freebayes_all_variants.txt", "freebayes")
) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = tool, values_from = present, values_fill = 0) %>%
  as.data.frame()

# Verificar a estrutura dos dados
str(all_data)

# Criar o gráfico UpSet com ajustes visuais
upset_plot <- upset(all_data, 
                    nsets = 4,
                    sets = c("DeepVariant", "GATK", "bcftools", "freebayes"),
                    order.by = "freq",
                    empty.intersections = "on",
                    mainbar.y.label = "Intersection Size",
                    sets.x.label = "Set Size",
                    text.scale = c(1.3, 1.3, 1, 1, 1.5, 1),
                    mb.ratio = c(0.7, 0.3),
                    point.size = 3.5,
                    line.size = 1.5,
                    keep.order = TRUE)

# Salvar o gráfico em alta resolução
png("upset_plot.png", width = 1200, height = 900, res = 120)
upset_plot
dev.off()

# Mostrar o gráfico
upset_plot

