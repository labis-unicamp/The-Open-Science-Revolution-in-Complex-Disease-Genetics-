#!/bin/bash
BIN_VERSION="1.6.1"
HOST_DIR="/str1/home/mariana.cavalheiro/Kaira"  # Confirmed path

docker run \
  -v "${HOST_DIR}:/data" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="/data/genome.fa" \
  --reads="/data/SRR1279767_bwa2_sorted.bam" \
  --output_vcf="/data/SRR1279767_deepvariant.vcf.gz" \
  --num_shards=$(nproc)
