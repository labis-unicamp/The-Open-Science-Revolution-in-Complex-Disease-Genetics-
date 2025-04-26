#!/bin/bash

# Configurações
OUTDIR="bam_benchmark_results/combined"
FLAGSTAT_DIR="bam_benchmark_results/flagstat"
STATS_DIR="bam_benchmark_results/stats"
OUTPUT_CSV="${OUTDIR}/comprehensive_metrics.csv"

# Cabeçalho do CSV
echo "Sample,Aligner,Total_Reads,Mapped_Reads,Mapped_Pct,Properly_Paired_Pct,Singletons_Pct,Avg_Insert_Size,Avg_Coverage,Mismatch_Rate,Avg_Quality,Error_Rate" > "$OUTPUT_CSV"

# Processar cada arquivo flagstat
for flagstat_file in "$FLAGSTAT_DIR"/*_flagstat.txt; do
    # Extrair nome da amostra e alinhador
    filename=$(basename "$flagstat_file")
    sample=$(echo "$filename" | cut -d'_' -f1)
    aligner=$(echo "$filename" | cut -d'_' -f2)
    
    # Extrair métricas do flagstat
    total_reads=$(grep "in total" "$flagstat_file" | awk '{print $1}')
    mapped_reads=$(grep "mapped (" "$flagstat_file" | awk '{print $1}')
    mapped_pct=$(grep "mapped (" "$flagstat_file" | awk '{print $5}' | tr -d '()%')
    properly_paired_pct=$(grep "properly paired" "$flagstat_file" | awk '{print $6}' | tr -d '()%')
    singletons_pct=$(grep "singletons" "$flagstat_file" | awk '{print $5}' | tr -d '()%')
    
    # Encontrar arquivo stats correspondente
    stats_file="$STATS_DIR/${sample}_${aligner}_sorted_stats.txt"
    
    # Extrair métricas do stats (se existir)
    if [ -f "$stats_file" ]; then
        avg_insert_size=$(grep "insert size average" "$stats_file" | awk '{print $5}')
        avg_coverage=$(grep "average length" "$stats_file" | awk '{print $5}')
        mismatch_rate=$(grep "mismatch rate" "$stats_file" | awk '{print $5}')
        avg_quality=$(grep "average quality" "$stats_file" | awk '{print $5}')
        error_rate=$(grep "error rate" "$stats_file" | awk '{print $5}')
    else
        avg_insert_size="NA"
        avg_coverage="NA"
        mismatch_rate="NA"
        avg_quality="NA"
        error_rate="NA"
    fi
    
    # Adicionar linha ao CSV
    echo "$sample,$aligner,$total_reads,$mapped_reads,$mapped_pct,$properly_paired_pct,$singletons_pct,$avg_insert_size,$avg_coverage,$mismatch_rate,$avg_quality,$error_rate" >> "$OUTPUT_CSV"
done

echo "Tabela de métricas gerada em: $OUTPUT_CSV"
