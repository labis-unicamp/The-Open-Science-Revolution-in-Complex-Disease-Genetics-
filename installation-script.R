# Script para instalar os pacotes necessários
# Execute este script apenas uma vez para instalar os pacotes

# Verificando e instalando o pacote VennDiagram
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}

# Verificando e instalando o pacote gridExtra (para manipulação de gráficos)
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}

# Verificando e instalando o pacote ggplot2 (opcional, para personalização adicional)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Verificando e instalando o pacote dplyr (para manipulação de dados)
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

print("Todos os pacotes necessários foram instalados com sucesso!")
