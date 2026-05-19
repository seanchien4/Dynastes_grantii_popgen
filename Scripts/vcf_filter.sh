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


# output depth frequency
awk 'NR>1 {          
  for(i=3;i<=NF;i++) 
    print $i         
}' DP_distri.DP.FORMAT |   
grep -v '\.' |        
sort -n |             
uniq -c |             
awk '{print $2"\t"$1}' > depth_frequency.txt


# vcf filteting

vcftools --gzvcf $unfilter_vcf \
--remove exclude_list.txt \
--remove-indels --maf 0.1 --max-missing 0.90 --minQ 30 \
--min-meanDP 10 --max-meanDP 30 \
--minDP 10 --maxDP 30 --recode --stdout | bgzip -c > 35cohort_geno_filter.g.vcf.gz

# check filtered vcf file 
# checking read depth for each individual
filtered_vcf=35cohort_geno_filter.g.vcf.gz
OUT=35cohort
vcftools --gzvcf $filtered_vcf --depth --out $OUT & 
awk '{sum += $3; count++} END {if (count >0) print sum / count}' $OUT.idepth
# missing data
vcftools --gzvcf $filtered_vcf --missing-indv --out $OUT &
# heterozygosity
vcftools --gzvcf $filtered_vcf --het --out $OUT &

module purge
module load GCC/13.2.0 BCFtools/1.19

echo $filtered_vcf > vcf_stats.txt
echo filtered_variant >> vcf_stats.txt
bcftools view -H $filtered_vcf | wc -l >> vcf_stats.txt &