bcftools merge --force-samples -o var_deepvariants.vcf.gz -O z SRR1279767_deepvariant.vcf.gz SRR1279773_deepvariant.vcf.gz SRR1279777_deepvariant.vcf.gz
bcftools merge --force-samples -o var_freebayes.vcf.gz -O z SRR1279767_bwa2_sorted_freebayes.vcf.gz  SRR1279773_bwa2_sorted_freebayes.vcf.gz  SRR1279777_bwa2_sorted_freebayes.vcf.gz
bcftools merge -o var_bcftools.vcf.gz -O z SRR1279767_bcftools.vcf.gz SRR1279773_bcftools.vcf.gz SRR1279777_bcftools.vcf.gz
bcftools merge -o var_gatk.vcf.gz -O z SRR1279767_gatk.vcf.gz SRR1279773_gatk.vcf.gz SRR1279777_gatk.vcf.gz
