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
    --minQ 30 \
    --maf 0.01 \
    --recode --stdout | bgzip > ${pop}_roh_filtered.g.vcf.gz &
done


# Plink
ml GCC/12.3.0 PLINK/2.00a3.7
for pop in Mogollon Lemmon Portal Utah; do
    echo "ROH: $pop..."
    
    # Convert VCF to PLINK format (if not already done)
    plink --vcf ${pop}_roh_filtered.g.vcf.gz \
        --make-bed \
        --out ${pop}_temp \
        --allow-extra-chr
    
    # Run ROH detection
    plink --bfile ${pop}_temp \
        --allow-extra-chr \
        --homozyg \
        --out roh_${pop}
    
    # Extract and format ROH regions (equivalent to RG lines)
    # PLINK outputs .hom and .hom.indiv files
    awk '{print $1, $3, $4, $5}' roh_${pop}.hom > ${pop}_roh_regions.txt
    
    # Clean up temporary files
    rm ${pop}_temp.*
done


cat Mogollon_roh_regions.txt Lemmon_roh_regions.txt \
    Portal_roh_regions.txt Utah_roh_regions.txt \
    > all_roh_regions.txt


# filtering all together
vcftools --gzvcf $vcf \
    --keep keep.list \
    --remove-indels \
    --min-alleles 2 \
    --max-alleles 2 \
    --max-missing 0.9 \
    --minQ 30 \
    --maf 0.01 \
    --recode --stdout | bgzip > all_roh_filtered.g.vcf.gz &


# check for how many SNPs 
ml purge
ml GCC/13.2.0 BCFtools/1.19
bcftools view -v snps -H Mogollon_roh_filtered.g.vcf.gz | wc -l

# ============================================================
# STEP 2: Call ROH 
# ============================================================

# --estimate-AF - \

for pop in Mogollon Lemmon Portal Utah; do

	echo "ROH: $pop..."
	# bcftools index ${pop}_roh_filtered.g.vcf.gz
	bcftools roh \
    -e ${pop}.list \
    -o roh_${pop}.txt \
    --threads 8 \
    ${pop}_roh_filtered.g.vcf.gz
  
  # Extract just the RG (region) lines
  grep "^RG" roh_${pop}.txt > ${pop}_roh_regions.txt
  
done

# ============================================================
# STEP 3: Combine for analysis
# ============================================================

cat Mogollon_roh_regions.txt Lemmon_roh_regions.txt \
    Portal_roh_regions.txt Utah_roh_regions.txt \
    > all_roh_regions.txt
