#!/bin/bash

echo "Iniciando processamento de arquivos BAM..."

for bam_file in *.bam; do
    # Pula arquivos que já são _sorted
    if [[ "$bam_file" == *_sorted.bam ]]; then
        echo " Pulando arquivo já ordenado: $bam_file"
        continue
    fi

    base_name="${bam_file%.bam}"
    sorted_bam="${base_name}_sorted.bam"
    
    echo "--------------------------------------------------"
    echo "Processando: $bam_file"
    
    # 1. Ordenar o BAM
    echo "Ordenando arquivo..."
    if ! samtools sort -@ $(nproc) "$bam_file" -o "$sorted_bam"; then
        echo "Falha ao ordenar $bam_file"
        continue  # Pula para o próximo arquivo em caso de erro
    fi
    
    # 2. Indexar o BAM ordenado
    echo "Indexando arquivo ordenado..."
    if ! samtools index "$sorted_bam"; then
        echo "Falha ao indexar $sorted_bam"
        continue
    fi
    
    # 3. Verificar integridade
    echo "Ordenação e indexação concluídas para:"
    echo "   - Arquivo ordenado: $sorted_bam"
    echo "   - Índice gerado:   ${sorted_bam}.bai"
    
    # 4. Opcional: remover o original não ordenado
    # read -p "Remover arquivo original não ordenado? (s/n) " -n 1 -r
    # echo
    # if [[ $REPLY =~ ^[Ss]$ ]]; then
    #     rm "$bam_file"
    #     echo "Arquivo original removido"
    # fi
done

echo "--------------------------------------------------"
echo "Processamento concluído! Arquivos gerados:"
ls -lh *_sorted.bam *_sorted.bam.bai
