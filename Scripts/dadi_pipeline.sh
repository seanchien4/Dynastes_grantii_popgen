#!/bin/bash

# dadi workflow
# filter out the vcf file
# vcf file
# filtering vcf second step 
module purge
module load GCC/12.2.0 GCC/12.3.0 VCFtools/0.1.16
unfilter_vcf=/scratch/user/schien/D_grantii/Population_Genetics/vcf/38cohort_geno.g.vcf.gz

# remove --min-alleles 2 to keep all sites
vcftools --gzvcf $unfilter_vcf \
    --keep keep.list \
    --min-alleles 2 --max-alleles 2 \
    --remove-indels \
    --recode --stdout | bgzip -c > step1.g.vcf.gz

vcftools --gzvcf step1.g.vcf.gz \
--max-missing 0.80 \
--minQ 30 \
--min-meanDP 10 --max-meanDP 30 \
--minDP 10 --maxDP 30 \
--recode --stdout | bgzip -c > dadi_geno_filter.g.vcf.gz

# checking how many SNPs


# LD 
ml purge
ml GCC/12.3.0 PLINK/2.00a3.7
plink --vcf dadi_geno_filter.g.vcf.gz --allow-extra-chr \
--make-bed --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 --out ld_prune

plink --vcf dadi_geno_filter.g.vcf.gz \
--allow-extra-chr --set-missing-var-ids @:# \
--extract ld_prune.prune.in --make-bed --recode vcf --out dadi_LD_prune

bgzip dadi_LD_prune.vcf

vcftools --gzvcf dadi_LD_prune.vcf.gz --depth --out dadi & 
vcftools --gzvcf dadi_LD_prune.vcf.gz --missing-indv --out dadi &


# filter combination
# Mogollon_Utah.list
# Lemmon_Portal

vcftools --gzvcf dadi_LD_prune.vcf.gz \
--max-missing 0.80 \
--keep Mogollon_Portal.list \
--recode --stdout | bgzip -c > Mogollon_Portal_geno_filter.g.vcf.gz &

vcftools --gzvcf Lemmon_Portal_geno_filter.g.vcf.gz --depth --out Lemmon_Portal &
vcftools --gzvcf Lemmon_Portal_geno_filter.g.vcf.gz --missing-indv --out Lemmon_Portal &


# prep the pop.map file 

# choosing projection
easySFS=/scratch/user/schien/Software/easySFS/easySFS.py
conda activate easySFS

$easySFS -i dadi_LD_prune.vcf.gz -p pop.map --preview

# Lemmon
# (2, 43) (3, 65) (4, 80) (5, 91) (6, 100)    (7, 97) (8, 104)    (9, 77) (10, 81)    

# Chiricahua
# (2, 40) (3, 60) (4, 73) (5, 83) (6, 91) (7, 78) (8, 83) (9, 9)  (10, 9) 

# Utah
# (2, 38) (3, 57) (4, 69) (5, 75) (6, 82) (7, 64) (8, 67) (9, 18) (10, 19)    

# Mogollon
# (2, 44) (3, 66) (4, 81) (5, 92) (6, 100)    (7, 108)    (8, 114)    (9, 120)    (10, 125)   (11, 129)   (12, 134)   (13, 137)   (14, 141)   (15, 144)   (16, 147(17, 147)   (18, 150)   (19, 133)   (20, 135)   (21, 83)    (22, 84)    (23, 24)    (24, 24)

$easySFS -i dadi_LD_prune.vcf.gz -p pop.map --proj 8,8,8,20


# using dadi

# activate dadi conda env 
# python3 dadi_models.py 


python3 dadi_models.py Utah-Mogollon.sfs 
# ================================================================================
# MODEL COMPARISON RESULTS
# ================================================================================

# Model                          LL           k      AIC          ΔAIC        
# --------------------------------------------------------------------------------
# Isolation + Size Change        -77.30       3      160.59       0.00        
# Strict Isolation               -77.38       3      160.77       0.18        

# ================================================================================
# BEST MODEL: Isolation + Size Change
# ================================================================================
# Log-Likelihood : -77.2954
# AIC            : 160.5907
# Num Parameters : 3

# Best-fit parameters:
#   nu1_final                 = 0.082176
#   nu2_final                 = 50.600301
#   T_split                   = 0.030239


python3 dadi_models.py Chiricahua-Mogollon.sfs
# ================================================================================
# MODEL COMPARISON RESULTS
# ================================================================================

# Model                          LL           k      AIC          ΔAIC        
# --------------------------------------------------------------------------------
# Isolation + Size Change        -84.89       3      175.79       0.00        
# Strict Isolation               -84.91       3      175.81       0.03        

# ================================================================================
# BEST MODEL: Isolation + Size Change
# ================================================================================
# Log-Likelihood : -84.8942
# AIC            : 175.7884
# Num Parameters : 3

# Best-fit parameters:
#   nu1_final                 = 0.419840
#   nu2_final                 = 3.308846
#   T_split                   = 0.027598


python3 dadi_models.py Lemmon-Mogollon.sfs 
# ================================================================================
# MODEL COMPARISON RESULTS
# ================================================================================

# Model                          LL           k      AIC          ΔAIC        
# --------------------------------------------------------------------------------
# Isolation + Size Change        -83.43       3      172.86       0.00        
# Strict Isolation               -83.44       3      172.88       0.02        

# ================================================================================
# BEST MODEL: Isolation + Size Change
# ================================================================================
# Log-Likelihood : -83.4314
# AIC            : 172.8628
# Num Parameters : 3

# Best-fit parameters:
#   nu1_final                 = 2.901890
#   nu2_final                 = 7.838110
#   T_split                   = 0.050775


python3 dadi_models.py Lemmon-Utah.sfs
# ================================================================================
# MODEL COMPARISON RESULTS
# ================================================================================

# Model                          LL           k      AIC          ΔAIC        
# --------------------------------------------------------------------------------
# Isolation + Size Change        -52.08       3      110.16       0.00        
# Strict Isolation               -52.31       3      110.61       0.45        

# ================================================================================
# BEST MODEL: Isolation + Size Change
# ================================================================================
# Log-Likelihood : -52.0811
# AIC            : 110.1622
# Num Parameters : 3

# Best-fit parameters:
#   nu1_final                 = 12.973941
#   nu2_final                 = 0.177100
#   T_split                   = 0.079343

python3 dadi_models.py Lemmon-Chiricahua.sfs
# ================================================================================
# MODEL COMPARISON RESULTS
# ================================================================================

# Model                          LL           k      AIC          ΔAIC        
# --------------------------------------------------------------------------------
# Isolation + Size Change        -56.95       3      119.91       0.00        
# Strict Isolation               -56.96       3      119.93       0.02        

# ================================================================================
# BEST MODEL: Isolation + Size Change
# ================================================================================
# Log-Likelihood : -56.9529
# AIC            : 119.9057
# Num Parameters : 3

# Best-fit parameters:
#   nu1_final                 = 8.884187
#   nu2_final                 = 0.332549
#   T_split                   = 0.024317


python3 dadi_models.py Chiricahua-Utah.sfs 
# ================================================================================
# MODEL COMPARISON RESULTS
# ================================================================================

# Model                          LL           k      AIC          ΔAIC        
# --------------------------------------------------------------------------------
# Isolation + Size Change        -49.16       3      104.32       0.00        
# Strict Isolation               -49.28       3      104.55       0.23        

# ================================================================================
# BEST MODEL: Isolation + Size Change
# ================================================================================
# Log-Likelihood : -49.1613
# AIC            : 104.3225
# Num Parameters : 3

# Best-fit parameters:
#   nu1_final                 = 10.275965
#   nu2_final                 = 0.279237
#   T_split                   = 0.134266

