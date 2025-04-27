# O link para download do hg38 está disponível em https://support.illumina.com/sequencing/sequencing_software/igenome.html
%%time
!wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
!tar -xf Homo_sapiens_UCSC_hg38.tar.gz

# codigo index
bowtie2-build referencia.fasta nome_do_indice_bowtie2
bwa index referencia.fasta
bwa-mem2 index referencia.fasta


# download from dbgap
while read sra_id; do prefetch --ngc "prj_32024.ngc" "$sra_id"; done < "arquivo_unico.txt"

# fastq dump
fastq-dump --gzip --split-files SRR1279777.sra

# alinhamento com o bwa-mem 2 e .sam to .bam
bwa-mem2 mem /home/kaira/powerjob/genome.fa SRR1279767_1.fastq.gz SRR1279767_2.fastq.gz | samtools view -b -o SRR1279767_bwa2.bam

# alinhamento com o bwa-mem e .sam to .bam
bwa mem /home/fciamponi/Kaira/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa SRR32660680_1.fastq.gz SRR32660680_2.fastq.gz | samtools view -bS -1 - > SRR32660680_mapped.bam

# alinhamento com o Bowtie2 e .sam to .bam
bowtie2 -x /home/fciamponi/Kaira/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.fa -1 SRR32660680_1.fastq.gz -2 SRR32660680_2.fastq.gz | samtools view -b -o SRR32660680_bowtie.bam

# sort e index com o samtools
samtools sort SRR1279777_bwa2.bam -o SRR1279777_bwa2_sorted.bam

samtools index SRR1279767_bwa2_sorted.bam
