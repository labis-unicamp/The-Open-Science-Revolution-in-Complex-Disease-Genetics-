#!/bin/bash
set -euo pipefail

# 1. First verify all required files exist
echo "=== Verifying input files ==="
[ -f "genome.fa" ] || { echo "Missing genome.fa"; exit 1; }
[ -f "genome.fa.fai" ] || { echo "Indexing genome.fa"; samtools faidx genome.fa; }
[ -f "genome.dict" ] || { echo "Creating genome.dict"; gatk CreateSequenceDictionary -R genome.fa -O genome.dict; }

[ -f "SRR1279767_bwa2_sorted.bam" ] || { echo "Missing BAM file"; exit 1; }
[ -f "SRR1279767_bwa2_sorted.bam.bai" ] || { echo "Indexing BAM file"; samtools index SRR1279767_bwa2_sorted.bam; }

# 2. Verify BAM read groups
echo -e "\n=== Checking BAM read groups ==="
if ! samtools view -H SRR1279767_bwa2_sorted.bam | grep -q '^@RG.*SM:'; then
    echo "Fixing read groups..."
    gatk AddOrReplaceReadGroups \
        I=SRR1279767_bwa2_sorted.bam \
        O=SRR1279767_fixed.bam \
        RGID=SRR1279767 \
        RGSM=SRR1279767 \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=unit1 \
        CREATE_INDEX=true
    INPUT_BAM="SRR1279767_fixed.bam"
else
    INPUT_BAM="SRR1279767_bwa2_sorted.bam"
fi

# 3. Run HaplotypeCaller with full debugging
echo -e "\n=== Running GATK HaplotypeCaller ==="
time gatk --java-options "-Xmx4G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    HaplotypeCaller \
    -R genome.fa \
    -I "$INPUT_BAM" \
    -O SRR1279767_gatk.vcf.gz

# 4. Index if successful
if [ $? -eq 0 ]; then
    echo -e "\n=== Indexing VCF ==="
    time tabix -p vcf SRR1279767_gatk.vcf.gz
else
    echo "GATK failed, check above errors"
    exit 1
fi

echo -e "\n=== Processing complete ==="
ls -lh SRR1279767_gatk.vcf.gz*
