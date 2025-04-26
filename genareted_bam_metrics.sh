#!/bin/bash

OUTDIR="bam_benchmark_results"
mkdir -p "$OUTDIR/flagstat" "$OUTDIR/stats" "$OUTDIR/qualimap_joint" "$OUTDIR/combined"

# 1. Gerar lista de BAMs (verificando se existem)
find . -maxdepth 1 -name "*.bam" > "$OUTDIR/bam_list.txt"
if [ ! -s "$OUTDIR/bam_list.txt" ]; then
    echo "ERRO: Nenhum arquivo .bam encontrado no diretório atual!"
    exit 1
fi

# 2. Rodar Qualimap multi-bamqc (com verificação de erro)
if command -v qualimap &>/dev/null; then
    echo "Executando Qualimap multi-bamqc..."
    qualimap multi-bamqc \
        -d "$OUTDIR/bam_list.txt" \
        -outdir "$OUTDIR/qualimap_joint" \
        -nt $(nproc) \
        --java-mem-size=8G 2>&1 | tee "$OUTDIR/qualimap_joint/qualimap_log.txt"
    
    if [ $? -ne 0 ]; then
        echo "AVISO: Qualimap falhou. Verifique o log em $OUTDIR/qualimap_joint/qualimap_log.txt"
    fi
else
    echo "AVISO: Qualimap não encontrado. Pulando análise Qualimap."
fi

# 3. Processar cada BAM individualmente
echo "filename,total_reads,mapped_reads,mapped_pct,properly_paired_pct,singletons_pct,avg_insert_size,avg_coverage,mismatch_rate,avg_quality" > "$OUTDIR/combined/summary.csv"

for bam_file in *.bam; do
    if [[ -f "$bam_file" ]]; then
        base_name=$(basename "$bam_file" .bam)
        echo "Processando $bam_file..."
        
        # Samtools flagstat
        echo "Rodando flagstat..."
        samtools flagstat "$bam_file" > "$OUTDIR/flagstat/${base_name}_flagstat.txt"
        
        # Samtools stats
        echo "Rodando samtools stats..."
        samtools stats "$bam_file" > "$OUTDIR/stats/${base_name}_stats.txt"
        
        # Extrair métricas
        echo "Extraindo métricas..."
        total_reads=$(grep "in total" "$OUTDIR/flagstat/${base_name}_flagstat.txt" | awk '{print $1}')
        mapped_reads=$(grep "mapped (" "$OUTDIR/flagstat/${base_name}_flagstat.txt" | awk '{print $1}')
        mapped_pct=$(grep "mapped (" "$OUTDIR/flagstat/${base_name}_flagstat.txt" | awk '{print $5}' | tr -d '()%')
        properly_paired_pct=$(grep "properly paired" "$OUTDIR/flagstat/${base_name}_flagstat.txt" | awk '{print $6}' | tr -d '()%')
        singletons_pct=$(grep "singletons" "$OUTDIR/flagstat/${base_name}_flagstat.txt" | awk '{print $5}' | tr -d '()%')
        
        avg_insert_size=$(grep "^SN" "$OUTDIR/stats/${base_name}_stats.txt" | grep "insert size average" | awk '{print $5}')
        avg_coverage=$(grep "^SN" "$OUTDIR/stats/${base_name}_stats.txt" | grep "average length" | awk '{print $5}')
        mismatch_rate=$(grep "^SN" "$OUTDIR/stats/${base_name}_stats.txt" | grep "mismatch rate" | awk '{print $5}')
        avg_quality=$(grep "^SN" "$OUTDIR/stats/${base_name}_stats.txt" | grep "average quality" | awk '{print $5}')
        
        # Adicionar ao CSV
        echo "$bam_file,$total_reads,$mapped_reads,$mapped_pct,$properly_paired_pct,$singletons_pct,$avg_insert_size,$avg_coverage,$mismatch_rate,$avg_quality" >> "$OUTDIR/combined/summary.csv"
    fi
done

echo "✅ Processamento concluído!"
echo "📋 Resultados em: $OUTDIR/combined/summary.csv"
if [ -d "$OUTDIR/qualimap_joint" ]; then
    echo "📊 Relatório Qualimap em: $OUTDIR/qualimap_joint/"
fi
