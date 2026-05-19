#!/bin/bash
# population structure
vcf=/scratch/user/schien/D_grantii/Population_Genetics/vcf/38cohort_geno.g.vcf.gz

# filter out low quality individual 
# only keep the individuals are used for analysis
ml purge
ml GCC/12.2.0 GCC/12.3.0 VCFtools/0.1.16

# ============================================================
# STEP 1: Filter with MAF
# ============================================================

for pop in Mogollon Lemmon Portal Utah; do
	echo "Filtering $pop..."
	vcftools --gzvcf $vcf \
	--keep ${pop}.list \
    --remove-indels \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing 0.9 \
    --minQ 30 \
    --maf 0.05 \
    --recode --stdout | bgzip > ${pop}_roh_filtered.g.vcf.gz &
done

# check for how many SNPs 
ml purge
ml GCC/13.2.0 BCFtools/1.19
bcftools view -v snps -H Mogollon.g.vcf.gz | wc -l

# ============================================================
# STEP 2: Call ROH 
# ============================================================

for pop in Mogollon Lemmon Portal Utah; do

	echo "RO $pop..."
	# bcftools index ${pop}_roh_filtered.g.vcf.gz
	bcftools roh \
	-e 
    ${pop}_roh_filtered.g.vcf.gz \
    -o roh_${pop}.txt
  
  # Extract just the RG (region) lines
  grep "^RG" roh_${pop}.txt > ${pop}_roh_regions.txt
  
  # Also extract individual summary (^IS lines)
  grep "^IS" roh_${pop}.txt > ${pop}_roh_summary.txt
done

# ============================================================
# STEP 3: Combine for analysis
# ============================================================

cat Mogollon_roh_regions.txt Lemmon_roh_regions.txt \
    Portal_roh_regions.txt Utah_roh_regions.txt \
    > all_roh_regions.txt


ml purge
ml GCC/13.2.0 BCFtools/1.19
# ROH
bcftools roh -e Mogollon.list -o roh_Mogollon.txt Mogollon.g.vcf.gz
grep "^RG" roh_Mogollon.txt > Mogollon_roh_regions.txt

bcftools roh -e Utah.list -o roh_Utah.txt Utah.g.vcf.gz &
grep "^RG" roh_Utah.txt > Utah_roh_regions.txt

bcftools roh -e Lemmon.list -o roh_Lemmon.txt Lemmon.g.vcf.gz &
grep "^RG" roh_Lemmon.txt > Lemmon_roh_regions.txt

bcftools roh -e Portal.list -o roh_Portal.txt Portal.g.vcf.gz &
grep "^RG" roh_Portal.txt > Portal_roh_regions.txt
