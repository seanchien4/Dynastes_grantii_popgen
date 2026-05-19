
# Sean Chien
# checking ID & variants 
module load GCC/13.2.0 BCFtools/1.19
# Show all individuals IDs
unfilter_vcf=../vcf/38cohort_geno.g.vcf.gz
module purge
module load GCC/12.2.0 GCC/12.3.0 VCFtools/0.1.16

# check the mane depth 
# remove ind has less average read depth 5 MC1 SR1 
# use mean plus 2 * standard deviation as maxDP
# remove that individual have over 20% missing data
# filter missing data 10% or 20% 


# vcf filteting

# 1 . remove individual and non-biallelic SNPs
vcftools --gzvcf $unfilter_vcf \
    --keep keep.list \
    --min-alleles 2 --max-alleles 2 \
    --remove-indels \
    --recode --stdout | bgzip -c > step1.g.vcf.gz


# 2. filter read depth and missing data 

vcftools --gzvcf step1.g.vcf.gz \
--maf 0.01 --max-missing 0.80 \
--minQ 20 \
--min-meanDP 10 --max-meanDP 30 \
--minDP 10 --maxDP 30 \
--recode --stdout | bgzip -c > final_geno_for_Dxy.g.vcf.gz

tabix final_geno_for_Dxy.g.vcf.gz

# Run pixy. This will calculate all stats in 10kb windows.
pixy --stats pi fst dxy --vcf final_geno_for_Dxy.g.vcf.gz --populations pop.map \
--window_size 5000 --output_folder . --bypass_invariant_check \
--n_cores 8 