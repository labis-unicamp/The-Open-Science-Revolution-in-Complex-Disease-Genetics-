# Unveiling Polygenic Pleiotropy with an Open-Source Genomics Pipeline: Decoding Shared Metabolic and Neurodegenerative Pathways in Alzheimer's Disease and Type 2 Diabetes through Open Science.
Open-source code for a bioinformatics pipeline analyzing genomic interactions between two polygenic diseases.

This repository contains all the code required to process raw data (such as FASTQ files), align it to the human reference genome using three different aligners, and generate variant files using four distinct variant callers. The accompanying manuscript for this repository will be published in the Journal of the X-Meeting 2025 Congress, critically analyzing each software’s performance and highlighting differences in their results. It also provides a step-by-step guide for conducting a GWAS (Genome-Wide Association Study) and a pleiotropy analysis, focusing on two diseases: Type 2 diabetes and Alzheimer’s disease.

## Downloading SRA files using prefetch and converting them to FASTQ format:
**Note:** For complete instructions on downloading dbGaP files, please refer to: [https://www.ncbi.nlm.nih.gov/sra/docs/sra-dbgap-download/]. Depending on the file type (whether they require an access key or not), they should be downloaded using one of the following methods:
``` 
#!/bin/bash

# Input file containing the list of SRA IDs
input_file="file_containing_the_sra_file_ids.txt"

# Directory where the SRA files will be downloaded
sra_dir="/path_to_the_destination_folder_of_the_generated_files"

# Directory where the FASTQ files will be generated
fastq_dir="/path_to_the_destination_folder_of_the_generated_files"

# NGC file for the prefetch command
ngc_file="/the_directory_path_where_your_access_key_is_stored_and_the_filename_of_your_access_key/prj_xxxxx.ngc"

# Loop through the input file
while read sra_id; do
    echo "Processing SRA ID: $sra_id"

    # Download the SRA file (aumentando o limite de tamanho para 100GB)
    prefetch --ngc "$ngc_file" --max-size 100G "$sra_id"

    # Convert the SRA file to FASTQ
    fastq-dump --gzip --split-files -O "$fastq_dir" "$sra_dir/$sra_id/$sra_id.sra"
    
    # Delete the SRA file
    if [ -f "$sra_dir/$sra_id/$sra_id.sra" ]; then
        rm -r "$sra_dir/$sra_id"
    fi

    echo "Processing complete for SRA ID: $sra_id"
done < "$input_file"
``` 
 
## Quality check:
This code performs quality control analysis on each FASTQ file in the directory, generating individual quality reports for every sample. It then consolidates all individual reports into a single comprehensive quality report in the output directory.

```
fastqc *.fastq *.fastq.gz --outdir=./fastqc_results
# inside the fastqc_results directory, the following commands were executed:
cd ./fastqc_results
multiqc fastqc_results/ -o multiqc_report --filename multiqc_report.html
```

## If the read quality of the samples falls below the required threshold, quality trimming and adapter removal can be performed using Trimmomatic:
```
trimmomatic PE -threads 4 input_R1.fastq.gz input_R2.fastq.gz output_R1_paired.fq.gz output_R1_unpaired.fq.gz output_R2_paired.fq.gz output_R2_unpaired.fq.gz ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
``` 

## Downloading the human reference genome (GRCh38/hg38) from official sources (https://support.illumina.com/sequencing/sequencing_software/igenome.html)
``` 
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
tar -xf Homo_sapiens_UCSC_hg38.tar.gz
``` 
## Generating index files for each alignment tool using the GRCh38/hg38 reference genome:
**Note:** Select an alignment tool and run the corresponding command for your chosen aligner.

### Bowtie index
```
bowtie2-build /Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
```
### BWA-MEM index
```
bwa index /Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
```
### BWA_MEM2 index
```
bwa-mem2 /Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
```

## Performing sequence alignment against the reference genome:
**Note:** Theses codes automatically converts SAM files to compressed BAM format during alignment.

### Bowtie alignment
```
bowtie2 -x /Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.fa -1 your_file_read_1.fastq.gz -2 your_file_read__2.fastq.gz | samtools view -b -o output_bowtie.bam
```
### BWA-MEM index
```
bwa mem /Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa your_file_read_1.fastq.gz your_file_read_2.fastq.gz | samtools view -bS -1 - > output_bwamem.bam
```
### BWA_MEM2 index
```
bwa-mem2 mem /Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa your_file_read_1.fastq.gz your_file_read_2.fastq.gz  | samtools view -b -o output_bwamem2.bam
```
## Visualizing the computational time required by each alignment tool:
**Note:** This code requires you to extract the execution time generated in the command line output and insert it into the script.
Refer to the available code at: /github_pipeline1/pipeline1/spider-plot.R
Two similar graphs to these ones will be generated. !github_pipeline1/pipeline1/grafico1_medias.png and !github_pipeline1/pipeline1/grafico2_tempo_real.png

## Post-alignment processing: Sorting and indexing BAM files using samtools:
```
samtools sort your_file.bam -o output_sorted.bam

samtools index output_sorted.bam
```
**Note:** For complete implementation details, refer to the alignment pipeline script: /github_pipeline1/pipeline1/from_downloading_the_reference_genome_to_bam_sorted.sh

## Variant Calling Workflow
**Here**, you can choose which variant calling tools you prefer. Tip: Read the article that generated this pipeline to make a better choice.
**Pay attention:** Remember to change the name and the way of your file.

### For variant call processing with **GATK** :
Refer to the available code at: /github_pipeline1/pipeline1/call_variant_gatk.sh
### For variant call  processing with **bcftools**:
Refer to the available code at: /github_pipeline1/pipeline1/script_call_variant_bcftools_freebayes.sh
### For variant call processing with **Freebayes**:
Refer to the available code at: /github_pipeline1/pipeline1/script_call_variant_bcftools_freebayes.sh
### For variant call processing with **DeepVariant**:
Refer to the available code at: /github_pipeline1/pipeline1/call_variant_deepvariant.sh

## Visualizing the computational time required by each variant calling tool:
**Note:** This code requires you to extract the execution time generated in the command line output and insert it into the script.
Refer to the available code at: /github_pipeline1/pipeline1/runingtime.R
A similar graph to this one will be generated. !github_pipeline1/pipeline1/runtime_comparison.png

