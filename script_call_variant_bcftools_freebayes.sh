#!/bin/bash

# Parte 1: Processamento com bcftools
#echo "Iniciando processamento com bcftools..."
#time (bcftools mpileup -f genome.fa SRR1279767_bwa2_sorted.bam | bcftools call -mv | bgzip > SRR1279767_bcftools.vcf.gz)
#time (bcftools mpileup -f genome.fa SRR1279773_bwa2_sorted.bam | bcftools call -mv | bgzip > SRR1279773_bcftools.vcf.gz)
#time (bcftools mpileup -f genome.fa SRR1279777_bwa2_sorted.bam | bcftools call -mv | bgzip > SRR1279777_bcftools.vcf.gz)

# Parte 2: Processamento com freebayes
echo "Iniciando processamento com freebayes..."
time (freebayes -f genome.fa SRR1279767_bwa2_sorted.bam | bgzip > SRR1279767_bwa2_sorted_freebayes.vcf.gz)
time (freebayes -f genome.fa SRR1279773_bwa2_sorted.bam | bgzip > SRR1279773_bwa2_sorted_freebayes.vcf.gz)
time (freebayes -f genome.fa SRR1279777_bwa2_sorted.bam | bgzip > SRR1279777_bwa2_sorted_freebayes.vcf.gz)

echo "Processamento concluído!"
