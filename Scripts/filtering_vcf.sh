#!/bin/bash
# Sean Chien
# checking ID & variants 
module load GCC/13.2.0 BCFtools/1.19
# Show all individuals IDs
unfilter_vcf=38cohort_geno.g.vcf.gz
OUT='38cohort'
bcftools query -l $unfilter_vcf

# creating a stats file to record
echo $unfilter_vcf > vcf_stats.txt
echo unfilter_variant >> vcf_stats.txt
bcftools view -H $unfilter_vcf | wc -l >> vcf_stats.txt & 

module purge
module load GCC/12.2.0 GCC/12.3.0 VCFtools/0.1.16

# checking read depth for each individual
vcftools --gzvcf $unfilter_vcf --depth --out $OUT & 
awk '{sum += $3; count++} END {if (count >0) print sum / count}' $OUT.idepth
# missing data
vcftools --gzvcf $unfilter_vcf --missing-indv --out $OUT &
# heterozygosity
vcftools --gzvcf $unfilter_vcf --het --out $OUT &


# 1. Calculate mean depth per individual
vcftools --gzvcf $unfilter_vcf --depth -c > depth_per_individual.tsv &

# 2. Calculate mean depth per site
vcftools --gzvcf $unfilter_vcf --site-mean-depth -c > depth_per_site.tsv &

# 3. Calculate missing data per individual
vcftools --gzvcf $unfilter_vcf --missing-indv -c > missing_per_individual.tsv &

# 4. Calculate missing data per site
vcftools --gzvcf $unfilter_vcf --missing-site -c > missing_per_site.tsv &


# check the mane depth 
# remove ind has less average read depth 5 MC1 SR1 
# use mean plus 2 * standard deviation as maxDP
# remove that individual have over 20% missing data
# filter missing data 10% or 20% 


# vcf filteting

# 1 . remove individual and non-biallelic SNPs
vcftools --gzvcf $unfilter_vcf \
    --remove exclude_list.txt \
    --min-alleles 2 --max-alleles 2 \
    --remove-indels \
    --recode --stdout | gzip -c > step1.g.vcf.gz


# 2. filter read depth and missing data 

vcftools --gzvcf step1.g.vcf.gz \
--maf 0.1 --max-missing 0.80 \
--minQ 20 \
--min-meanDP 10 --max-meanDP 30 \
--minDP 10 --maxDP 30 \
--recode --stdout | gzip -c > final_geno_filter.g.vcf.gz


vcftools --gzvcf final_geno_filter.g.vcf.gz --depth --out final & 
vcftools --gzvcf final_geno_filter.g.vcf.gz --missing-indv --out final &

module purge
module load GCC/13.2.0 BCFtools/1.19

echo $filtered_vcf > vcf_stats.txt
echo filtered_variant >> vcf_stats.txt
bcftools view -H $filtered_vcf | wc -l >> vcf_stats.txt &