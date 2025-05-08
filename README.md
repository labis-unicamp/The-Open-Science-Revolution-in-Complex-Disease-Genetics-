# Unveiling Polygenic Pleiotropy with an Open-Source Genomics Pipeline: Decoding Shared Metabolic and Neurodegenerative Pathways in Alzheimer's Disease and Type 2 Diabetes through Open Science.
Open-source code for a bioinformatics pipeline analyzing genomic interactions between two polygenic diseases.

This repository contains all the code required to process raw data (such as FASTQ files), align it to the human reference genome using three different aligners, and generate variant files using four distinct variant callers. The accompanying manuscript for this repository will be published in the Journal of the X-Meeting 2025 Congress, [https://labis.cbmeg.unicamp.br/labis/publicacoes/131-the-open-science-revolution-in-complex-disease-genetics-an-integrated-pipeline-from-fastq-to-gwas-and-functional-pleiotropy], critically analyzing each software’s performance and highlighting differences in their results. It also provides a step-by-step guide for conducting a GWAS (Genome-Wide Association Study) and a pleiotropy analysis, focusing on two diseases: Type 2 diabetes and Alzheimer’s disease.

The code snippets explicitly shown in this file are written in Bash. The remaining scripts are implemented in:
        
* R (files with .R extension)

* Python (.py or .ipynb files)

* Bash (.sh files)

## Download the bioinfo2_env.yml file to your local machine and create and activate the Conda environment using the following commands:
```
conda env create -f bioinfo2_env.yml
```

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

## Generating the BAI index file:
Refer to the available code at: /github_pipeline1/pipeline1/sorted_bai_bam_files.sh


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
## Generating quality metrics for BAM files:
Refer to the available code at: /github_pipeline1/pipeline1/genareted_bam_metrics.sh

## Building a table of key metrics:
Refer to the available code at: /github_pipeline1/pipeline1/table_bam_metrics.sh

## Visualizing the computational time required by each alignment tool:
**Note:** This code requires you to extract the execution time generated in the command line output and insert it into the script.
Refer to the available code at: /github_pipeline1/pipeline1/spider-plot.R
Two similar graphs to these ones will be generated. 

<div align="center">
<img src="https://github.com/user-attachments/assets/4e6c3a8c-5f62-4827-94ed-80d6e36deffc" width="300px" />
</div> 

<div align="center">
<img src="https://github.com/user-attachments/assets/9732a444-c773-4c83-8234-524f3d657310" width="300px" />
</div> 

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
A similar graph to this one will be generated. 
<div align="center">
<img src="https://github.com/user-attachments/assets/cf57fd15-7525-47cc-a8fe-177408352f76" width="300px" />
</div> 

## To merge multiple sample-specific VCF files, use bcftools with the following command:

Before merging, it's necessary to perform the indexing:
```
for file in *.vcf.gz; do
    bcftools index "$file"
done
```

For VCF files generated by FreeBayes or DeepVariant:
```
bcftools merge --force-samples -o output_name.vcf.gz -O z file1.vcf.gz file2.vcf.gz file3.vcf.gz
```
For VCF files generated by GATK or bcftools:
```
bcftools merge -o output_name.vcf.gz -O z file1.vcf.gz file2.vcf.gz file3.vcf.gz
```

## Extracting variants from VCF files:
Refer to the available code at: /github_pipeline1/pipeline1/extract_var.sh

## Upset plot to compare variant counts across different callers:
Refer to the available code at: /github_pipeline1/pipeline1/upset-plot-script.R
Remember to replace filenames and paths in the command. This will generate a plot similar to:
<div align="center">
<img src="https://github.com/user-attachments/assets/ea691bd5-51c9-4d25-9cc2-cada2e9b0fe3" width="300px" />
</div> 

## Inter-caller concordance was additionally assessed using **Venn diagrams** 
See implementation at: /github_pipeline1/pipeline1/venn-diagram-script.R

<div align="center">
<img src="https://github.com/user-attachments/assets/a5411cd4-0652-4cd2-8e07-66cca43a32b1" width="300px" />
</div> 

## **Additional option:** Optimized pipeline for cluster environments that processes downloaded SRA files into VCFs using BWA-MEM2 and FreeBayes as core tools:
Refer to the available code at: /github_pipeline1/pipeline1/pipeline_hyperopt_otimizado.sh

## GWAS
First of all, in order to process GWAS, you must have a vcf file containing the phenotype samples of interest and the control samples and .txt file with the names of the samples positioned in two equal columns, with the number '2' in front of the phenotype sample and the number '1' in front of the control sample. See the example .txt files:

* /github_pipeline1/pipeline1/dia_pheno_example.txt
* /github_pipeline1/pipeline1/ad_pheno_example.txt

### After that, you should filter out the variants with statistical power equal to or greater than the 5% allelic frequency:
```
bcftools view -S samples.txt -i 'MAF[0] > 0.0005' -Oz -o output_file.vcf.gz your_input_file.vcf.gz
```
### Convert to PLINK format and additional filters:
```
plink --vcf your_input_file.vcf.gz --make-bed --maf 0.05 --geno 0.05 --hwe 1e-6 --allow-extra-chr --out your_file_clean
```
### Execute the GWAS analysis:
```
plink --bfile your_file_clean --pheno file_phenotype.txt --assoc --allow-no-sex --out gwas_results
```

### Checking the results:
```
# Checking results:
#1. Check first lines of the file:
head -n 5 gwas_results.assoc
#2. Verify column names:
awk 'NR==1 {print; exit}' gwas_results.assoc
#3. Count total SNPs:
wc -l gwas_results.assoc  # Total lines
awk 'NR>1' gwas_results.assoc | wc -l  # SNPs only (excluding header)
#4. Check p-value distribution:
# For non-missing p-values:
awk '$9 != "NA" && NR>1 {print $9}' gwas_results.assoc | sort -g | uniq -c | head
# Check problematic values
# Non-numeric values in P column (usually column 9):
awk 'NR>1 && $9 !~ /^[0-9.eE-]+$/ {print $9}' gwas_results.assoc | sort | uniq -c
```

### Generating individual Manhattan and QQ plots using R:
```
#!/usr/bin/env Rscript

# Install package if needed
if (!require("qqman")) {
    install.packages("qqman", repos="https://cloud.r-project.org")
    library(qqman)
}

# Read data
gwas_data <- read.table("gwas_results.assoc", header=TRUE)

# Create unique SNP IDs (since SNP column contains ".")
gwas_data$SNP_ID <- paste0("chr", gwas_data$CHR, ":", gwas_data$BP)

# Check and clean p-values
gwas_data$P <- as.numeric(gwas_data$P)
valid_data <- gwas_data[!is.na(gwas_data$P) & is.finite(gwas_data$P), ]

# 1. Manhattan Plot
png("manhattan_final.png", width=1200, height=600)
manhattan(valid_data,
          chr = "CHR",
          bp = "BP",
          p = "P",
          snp = "SNP_ID",  # Use created IDs
          main = "Manhattan Plot - GWAS Results",
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8))
dev.off()

# 2. QQ Plot
png("qq_final.png", width=600, height=600)
qq(valid_data$P, main = "QQ Plot - GWAS Results")
dev.off()

print("Plots generated successfully! Check: manhattan_final.png and qq_final.png")
```
## Pleiotropy analysis:
### Processe o GWAS para duas doenças em questão e então com os resultados obtidos realizar as análises:
```
# Step 1: Identify significant SNPs in both diseases
# SNPs with p < 0.05 in both (adjust threshold as needed)
awk 'NR==FNR && $9 < 0.05 {print $2; next}
     NR!=FNR && $9 < 0.05 && ($2 in snps) {print $2}' \
     gwas_disease1_results.assoc gwas_disease2_results.assoc > snps_pleiotropicos.txt
```

```
# Step 2: Statistical analysis of pleiotropy
# Option 1: Heterogeneity test (Cochran's Q)
# Combine results (requires columns: SNP, P_DISEASE1, P_DISEASE2)
paste <(awk 'NR>1 {print $2, $9}' gwas_disease1_results.assoc) \
      <(awk 'NR>1 {print $9}' gwas_disease2_results.assoc) > combined_pvalues.txt
# Remove rows that don't have 3 columns
# Create filtered version (keeps only rows with 3 columns)
awk 'NF==3' combined_pvalues.txt > combined_pvalues_filtered.txt

# Run directly in terminal
Rscript -e '
  data <- read.table("combined_pvalues_filtered.txt", col.names=c("SNP", "P1", "P2"))
  data$Q <- -2*(log(data$P1) + log(data$P2))  # Q statistic
  data$P_pleio <- pchisq(data$Q, df=2, lower.tail=FALSE)
  write.table(data, "pleiotropy_results.txt", quote=F, row.names=F)'

```
```
# Option 2: Effect correlation (β)
# Extract effect sizes (OR or β)
paste <(awk 'NR>1 {print $2, $10}' gwas_disease1_results.assoc) \
      <(awk 'NR>1 {print $10}' gwas_disease2_results.assoc) > combined_effects.txt

# Run directly in terminal

# Filter rows with exactly 3 columns (SNP + both effects)
awk 'NF==3' combined_effects.txt > combined_effects_clean.txt

# Calculate Spearman correlation in R
Rscript -e '
  data <- read.table("combined_effects_clean.txt", col.names=c("SNP", "OR1", "OR2"))
  cor_test <- cor.test(data$OR1, data$OR2, method="spearman")
  print(paste("Effect correlation (rho):", cor_test$estimate))'
```

## Manhattan plot combine:
Refer to the available code at: /github_pipeline1/pipeline1/Manhattan-plot-combine.R

<div align="center">
<img src="https://github.com/user-attachments/assets/3ada9047-aa21-4018-b47a-213419be71e0" width="300px" />
</div> 

## Scatter plot of p-values:
Refer to the available code at: /github_pipeline1/pipeline1/Scatter-plot.R

<div align="center">
<img src="https://github.com/user-attachments/assets/3cf3a956-cd8e-41c8-8c61-eddf4eee5f29" width="300px" />
</div> 

## The variants were annotated in the VEP via a graphical interface on the website: 
https://www.ensembl.org/Homo_sapiens/Tools/VEP?tl=k2dBMwVpFcRHNgnD-10806356

### the variants were filtered using the python script available at: 
github_pipeline1/pipeline1/find_biomarkers.ipynb

## REVIGO analysis:
The onthology gene codes of each variant were collected in the previous script and the analysis in the code was processed via the graphical interface. To generate the main image contained in the article, the R script was extracted from the journal and processed.
Refer to the available code at: /github_pipeline1/pipeline1/scriptREVIGO_treemapBP.R

<div align="center">
<img src="https://github.com/user-attachments/assets/9496d4a6-2684-4e67-b3b5-9baf05a0f5f6" width="300px" />
</div> 


# References:
* deepseek : https://www.deepseek.com/
* chatGPT: https://openai.com/
* claude.ai : https://claude.ai/new
* https://www.youtube.com/watch?v=oMFiGEZ6UlQ
* VEP : https://www.ensembl.org/Homo_sapiens/Tools/VEP?tl=k2dBMwVpFcRHNgnD-10806356
* REVIGO : http://revigo.irb.hr/

