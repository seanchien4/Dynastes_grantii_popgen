#!/bin/bash
# Sean Chien
# GEA 

module purge
module load GCC/12.2.0 GCC/12.3.0 VCFtools/0.1.16

# vcf filteting

# 1 . remove individual and non-biallelic SNPs
vcftools --gzvcf $unfilter_vcf \
    --keep keep.list \
    --min-alleles 2 --max-alleles 2 \
    --remove-indels \
    --recode --stdout | bgzip -c > step1.g.vcf.gz


# 2. filter read depth and missing data 

vcftools --gzvcf ../Dxy/step1.g.vcf.gz \
--maf 0.01 --max-missing 0.90 \
--minQ 20 \
--min-meanDP 10 --max-meanDP 30 \
--minDP 10 --maxDP 30 \
--recode --stdout | bgzip -c > final_geno_GEA.g.vcf.gz


ml purge
ml GCC/12.3.0 PLINK/2.00a3.7
# --recod for .map
# --recod A for .raw

plink --vcf final_geno_GEA.g.vcf.gz \
--maf 0.01 --geno 0.1 \
--recode A \
--output-missing-genotype 9 \
--out genotype_data \
--allow-extra-chr \
--double-id 



# extract SNPs that is specific for a variable 


# extract SNPs gene IDs from braker.gff3
awk '$3=="gene"' ../BRAKER/braker.gff3 | \
awk 'NR==FNR {
    snp_chr[$1 ":" $2]
    next
}
{
    split($1, chr_parts, "_")
    gff_chr = chr_parts[1]

    for (s in snp_chr) {
        split(s, snp_parts, ":")
        if (snp_parts[1] == gff_chr &&
            snp_parts[2] >= $4 &&
            snp_parts[2] <= $5) {

            match($9, /ID=([^;]+)/, id)
            print snp_parts[1] "\t" snp_parts[2] "\t" id[1]
        }
    }
}' PC_SNP.txt - | awk '{print $3}' | sort | uniq > PC_gene.list

