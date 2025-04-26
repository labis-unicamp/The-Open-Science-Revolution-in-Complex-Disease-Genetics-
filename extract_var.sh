#!/bin/bash

# Lista de arquivos VCF de entrada
vcf_files=(
  "var_bcftools.vcf.gz"
  "var_deepvariants.vcf.gz"
  "var_freebayes.vcf.gz"
  "var_gatk.vcf.gz"
)

# Verificar se o bcftools está instalado
if ! command -v bcftools &> /dev/null; then
  echo "Erro: bcftools não está instalado ou não está no PATH"
  echo "Instale com: sudo apt install bcftools (Ubuntu) ou conda install -c bioconda bcftools"
  exit 1
fi

# Processar cada arquivo VCF
for vcf in "${vcf_files[@]}"; do
  # Verificar se o arquivo existe e não está vazio
  if [[ ! -s "$vcf" ]]; then
    echo -e "\n\033[31mERRO: Arquivo $vcf não existe ou está vazio\033[0m"
    continue
  fi
  
  # Verificar se o arquivo é um VCF válido
  if ! bcftools view "$vcf" &> /dev/null; then
    echo -e "\n\033[31mERRO: $vcf não é um arquivo VCF válido ou está corrompido\033[0m"
    continue
  fi

  # Nome do arquivo de saída
  output_txt="${vcf%.vcf.gz}_all_variants.txt"
  
  echo -e "\n\033[32mProcessando $vcf -> $output_txt\033[0m"
  
  # Extrair todas as variantes (incluindo múltiplos ALTs se existirem)
  bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "$vcf" > "$output_txt"
  
  # Verificar se a extração foi bem-sucedida
  if [[ $? -ne 0 ]]; then
    echo -e "\033[31mFalha na extração de $vcf\033[0m"
  else
    # Contar variantes extraídas
    count=$(wc -l < "$output_txt")
    if [[ $count -eq 0 ]]; then
      echo -e "\033[33mAVISO: Nenhuma variante encontrada em $vcf\033[0m"
      echo "Possíveis causas:"
      echo "1. Arquivo VCF realmente sem variantes"
      echo "2. Filtros aplicados inadvertidamente"
      echo "3. Formato incorreto do VCF"
      
      # Verificar se há variantes no arquivo original
      echo -e "\nVerificando conteúdo original:"
      bcftools view "$vcf" | grep -v "^#" | head -n 5
    else
      echo -e "\033[32mSucesso: $count variantes extraídas\033[0m"
      echo -e "Exemplo de variantes:"
      head -n 3 "$output_txt" | sed 's/^/  /'
    fi
  fi
done

echo -e "\n\033[1mProcessamento concluído. Verifique os arquivos de saída.\033[0m"
