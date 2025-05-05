# Unveiling Polygenic Pleiotropy with an Open-Source Genomics Pipeline: Decoding Shared Metabolic and Neurodegenerative Pathways in Alzheimer's Disease and Type 2 Diabetes through Open Science.
Open-source code for a bioinformatics pipeline analyzing genomic interactions between two polygenic diseases.

This repository contains all the code required to process raw data (such as FASTQ files), align it to the human reference genome using three different aligners, and generate variant files using four distinct variant callers. The accompanying manuscript for this repository will be published in the Journal of the X-Meeting 2025 Congress, critically analyzing each software’s performance and highlighting differences in their results. It also provides a step-by-step guide for conducting a GWAS (Genome-Wide Association Study) and a pleiotropy analysis, focusing on two diseases: Type 2 diabetes and Alzheimer’s disease.

## Downloading SRA files using prefetch and converting them to FASTQ format:
### **Note:** For complete instructions on downloading dbGaP files, please refer to: [https://www.ncbi.nlm.nih.gov/sra/docs/sra-dbgap-download/]. Depending on the file type (whether they require an access key or not), they should be downloaded using one of the following methods:
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
