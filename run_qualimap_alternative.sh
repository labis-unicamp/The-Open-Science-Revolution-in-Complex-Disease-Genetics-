#!/bin/bash
# run_qualimap_alternative.sh

OUTDIR="bam_benchmark_results/qualimap_individual"
mkdir -p "$OUTDIR"

# Processar cada BAM individualmente
for bam in *_sorted.bam; do
    sample=$(basename "$bam" _sorted.bam)
    echo "Processando $sample..."
    
    qualimap bamqc \
        -bam "$bam" \
        -outdir "$OUTDIR/$sample" \
        --java-mem-size=2000M 2>&1 | tee "$OUTDIR/${sample}_log.txt"
done

# Consolidar resultados manualmente
echo "Sample,Mapped_Reads,Mapped_Pct,Avg_Coverage,Avg_Quality" > "$OUTDIR/summary.csv"
for report in "$OUTDIR"/*/genome_results.txt; do
    sample=$(echo "$report" | cut -d'/' -f3)
    mapped=$(grep "number of mapped reads" "$report" | cut -d'=' -f2 | xargs)
    mapped_pct=$(grep "percentage of mapped reads" "$report" | cut -d'=' -f2 | xargs)
    coverage=$(grep "mean coverageData" "$report" | cut -d'=' -f2 | xargs)
    quality=$(grep "mean quality" "$report" | cut -d'=' -f2 | xargs)
    
    echo "$sample,$mapped,$mapped_pct,$coverage,$quality" >> "$OUTDIR/summary.csv"
done

echo "Processamento concluído!"
echo "Relatórios individuais em: $OUTDIR/"
echo "Sumário consolidado em: $OUTDIR/summary.csv"
