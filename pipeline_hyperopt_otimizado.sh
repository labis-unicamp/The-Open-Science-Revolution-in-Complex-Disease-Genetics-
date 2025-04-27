#!/bin/bash -l
#$ -N ngs_hyperopt
#$ -o logs/${JOB_NAME}_${TASK_ID}.out
#$ -e logs/${JOB_NAME}_${TASK_ID}.err
#$ -t 1-1
#$ -tc 2
#$ -l h_vmem=24G
#$ -pe smp 12

# ======== ATIVAÇÃO DO CONDA ========
source /home/kaira/miniconda/bin/activate
conda activate /home/kaira/miniconda/envs/bioinfo2

echo "================================================"
echo "Pipeline HyperOpt - Análise NGS"
echo "Início: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Host: $(hostname)"
echo "Job ID: $JOB_ID"
echo "Task ID: $SGE_TASK_ID"
echo "================================================"

# ======== CONFIGURAÇÕES OTIMIZADAS ========
SRA_LIST="/home/kaira/powerjob/grupo_unico.txt"
REF_GENOME="/home/kaira/powerjob/genome.fa"
NGC_FILE="/home/kaira/powerjob/prj_32024.ngc"
WORKDIR="/home/kaira/pipeline_results"
THREADS=$NSLOTS
THREADS_SORT=$((NSLOTS/2))
THREADS_FREEBAYES=$((NSLOTS/4))

# Parâmetros de qualidade otimizados
MIN_ALT_COUNT=1
MIN_ALT_FRACTION=0.01
BWA_MEM_PARAMS="-K 100000000 -v 3 -t $THREADS"
SAMTOOLS_PARAMS="--output-fmt BAM -@ $THREADS_SORT"
FREEBAYES_PARAMS="--standard-filters --min-alternate-count $MIN_ALT_COUNT --min-alternate-fraction $MIN_ALT_FRACTION -j $THREADS_FREEBAYES"

# Controle de execução
CLEANUP=${CLEANUP:-1}  # 1=limpeza ativada, 0=modo debug

# ======== FUNÇÕES AUXILIARES ========
validate_dependencies() {
    echo "=== Validando dependências ==="
    local deps=("prefetch" "fastq-dump" "bwa-mem2" "samtools" "freebayes" "bgzip" "tabix")
    for cmd in "${deps[@]}"; do
        if ! command -v "$cmd" &>/dev/null; then
            echo "[ERRO] '$cmd' não encontrado no PATH" >&2
            exit 1
        fi
        echo "[OK] $cmd: $(which $cmd)"
    done
}

validate_files() {
    echo "=== Validando arquivos de entrada ==="
    [ -f "$SRA_LIST" ] || { echo "[ERRO] Arquivo SRA_LIST não encontrado" >&2; exit 1; }
    [ -f "$REF_GENOME" ] || { echo "[ERRO] Genoma de referência não encontrado" >&2; exit 1; }
    [ -f "$REF_GENOME.fai" ] || { echo "[ERRO] Índice do genoma faltando. Execute: samtools faidx $REF_GENOME" >&2; exit 1; }
    echo "[OK] Todos os arquivos de entrada validados"
}

# ======== FUNÇÕES DO PIPELINE OTIMIZADAS ========
download_sra() {
    local sra=$1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ETAPA 1/4 - DOWNLOAD: $sra"
    
    if [ -f "${sra}.sra" ]; then
        echo "[INFO] Arquivo ${sra}.sra já existe. Pulando download."
        return 0
    fi

    if ! prefetch --ngc "$NGC_FILE" \
              --max-size 100G \
              --progress \
              --verify yes \
              --type sra \
              --output-file "${WORKDIR}/${sra}.sra" \
              "$sra"; then
        echo "[ERRO CRÍTICO] Prefetch falhou para $sra" >&2
        exit 1
    fi
    echo "[SUCESSO] Download concluído: $(du -h ${sra}.sra)"
}

convert_to_fastq() {
    local sra=$1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ETAPA 2/4 - CONVERSÃO: $sra"
    
    if [ -f "${sra}_1.fastq.gz" ] && [ -f "${sra}_2.fastq.gz" ]; then
        echo "[INFO] Arquivos FASTQ já existem. Pulando conversão."
        return 0
    fi
    
    if ! fastq-dump --gzip --split-files --skip-technical --clip "${sra}.sra"; then
        echo "[ERRO CRÍTICO] Conversão FASTQ falhou" >&2
        exit 1
    fi
    echo "[SUCESSO] FASTQs gerados: $(ls -lh ${sra}_?.fastq.gz)"
}

align_reads() {
    local sra=$1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ETAPA 3/4 - ALINHAMENTO: $sra"
    
    if [ -f "${sra}_sorted.bam" ]; then
        echo "[INFO] BAM alinhado já existe. Pulando alinhamento."
        return 0
    fi
    
    # Pipeline otimizado: alinhamento e conversão para BAM em um passo
    bwa-mem2 mem $BWA_MEM_PARAMS "$REF_GENOME" \
        "${sra}_1.fastq.gz" "${sra}_2.fastq.gz" | \
    samtools sort $SAMTOOLS_PARAMS -o "${sra}_sorted.bam" && \
    samtools index "${sra}_sorted.bam" || {
        echo "[ERRO] Alinhamento falhou para $sra" >&2
        exit 1
    }
    
    echo "[SUCESSO] BAM gerado: $(du -h ${sra}_sorted.bam)"
}

call_variants() {
    local sra=$1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ETAPA 4/4 - CHAMADA DE VARIANTES: $sra"
    
    mkdir -p vcf
    if [ -f "vcf/${sra}.vcf.gz" ]; then
        echo "[INFO] VCF já existe. Pulando chamada de variantes."
        return 0
    fi
    
    # Comando corrigido (sem -j e com redirecionamento de erros):
    freebayes --standard-filters \
              --min-alternate-count "$MIN_ALT_COUNT" \
              --min-alternate-fraction "$MIN_ALT_FRACTION" \
              -f "$REF_GENOME" "${sra}_sorted.bam" 2> "vcf/${sra}.freebayes.err" | \
    bgzip -@ $THREADS > "vcf/${sra}.vcf.gz" && \
    tabix -p vcf "vcf/${sra}.vcf.gz" || {
        echo "[ERRO CRÍTICO] FreeBayes falhou" >&2
        exit 1
    }
    echo "[SUCESSO] VCF gerado: $(du -h vcf/${sra}.vcf.gz)"
}

cleanup_files() {
    local sra=$1
    if [ "$CLEANUP" -eq 1 ]; then
        echo "[INFO] Limpando arquivos intermediários..."
        rm -f \
            "${sra}.sra" \
            "${sra}_1.fastq.gz" "${sra}_2.fastq.gz" \
            "${sra}_bwa2.bam" \
            "${sra}".*.err \
            "${sra}_regions.txt"
    else
        echo "[INFO] Modo debug: Arquivos preservados"
        echo "Arquivos restantes: $(ls -lh ${sra}*)"
    fi
}

# ======== EXECUÇÃO PRINCIPAL OTIMIZADA ========
main() {
    # Configuração inicial
    mkdir -p "${WORKDIR}/vcf" "${WORKDIR}/logs"
    cd "$WORKDIR"
    SRA_ID=$(sed -n "${SGE_TASK_ID}p" "$SRA_LIST")
    
    echo "Processando SRA ID: $SRA_ID"
    echo "Diretório de trabalho: $(pwd)"
    echo "Threads disponíveis: $THREADS"
    echo "Threads para sort: $THREADS_SORT"
    echo "Threads para freebayes: $THREADS_FREEBAYES"
    
    validate_dependencies
    validate_files
    
    # Pipeline principal com verificações de arquivos existentes
    download_sra "$SRA_ID"
    convert_to_fastq "$SRA_ID"
    align_reads "$SRA_ID"
    call_variants "$SRA_ID"
    
    # Limpeza condicional
    cleanup_files "$SRA_ID"
    
    echo "================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] CONCLUSÃO: Pipeline finalizado para $SRA_ID"
    echo "================================================"
}

main "$@"
