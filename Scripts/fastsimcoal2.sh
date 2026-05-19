#!/bin/bash
# fastsimcoal2


awk '{OFS="\t"; print $1, "0", $2}' ../ref/GCA_029619325.2_ASM2961932v2_genomic.fna.fai > all_scaffolds.bed
awk '$2 > 20000000 {print $1}' ../ref/GCA_029619325.2_ASM2961932v2_genomic.fna.fai  > large_scaffolds.txt
grep -wFf large_scaffolds.txt all_scaffolds.bed > large_scaffolds.bed


# filtering VCF file 
ml purge
ml GCC/12.2.0 GCC/12.3.0 VCFtools/0.1.16
unfilter_vcf=/scratch/user/schien/D_grantii/Population_Genetics/vcf/38cohort_geno.g.vcf.gz
# 1 . remove individual and non-biallelic SNPs
vcftools --gzvcf $unfilter_vcf \
    --keep keep.list \
    --min-alleles 2 --max-alleles 2 \
    --remove-indels \
    --recode --stdout | bgzip -c > step1.g.vcf.gz

tabix step1.g.vcf.gz
ml purge 
ml GCC/13.2.0 BCFtools/1.19
bcftools view  --regions-file large_scaffolds.bed step1.g.vcf.gz -O z -o step1_large_scaffolds.vcf.gz
tabix step1_large_scaffolds.vcf.gz 

# 2. filter read depth and missing data 
ml purge
ml GCC/12.2.0 GCC/12.3.0 VCFtools/0.1.16
vcftools --gzvcf step1_large_scaffolds.vcf.gz \
--max-missing 1 \
--minQ 30 \
--min-meanDP 10 --max-meanDP 30 \
--minDP 10 --maxDP 30 \
--minGQ 20 \
--recode --stdout | bgzip -c > final_geno_fastsimcoal.g.vcf.gz

tabix final_geno_fastsimcoal.g.vcf.gz

# 3. Keep only one SNP per locus (if RADseq)
# to avoid linkage
vcftools --gzvcf final_geno_fastsimcoal.g.vcf.gz \
--thin 1000 \
--recode \
--out filtered_thinned

# Check how many SNPs remain
grep -v "^#" filtered_thinned.recode.vcf | wc -l

bgzip filtered_thinned.recode.vcf
tabix filtered_thinned.recode.vcf.gz


# create pop.map 

# Preview projections
easySFS=/scratch/user/schien/Software/easySFS/easySFS.py
conda activate easySFS
$easySFS \
-i filtered_thinned.recode.vcf.gz \
-p pop.map \
--preview \
-a

(easySFS) [schien@grace4 fastsimcoal2]$ $easySFS \
> -i filtered_thinned.recode.vcf.gz \
> -p pop.map \
> --preview \
> -a

# Graham
# (2, 3507)	

# Lemmon
# (2, 3300)	(3, 4951)	(4, 6061)	(5, 6902)	(6, 7580)	(7, 8146)	(8, 8632)	(9, 9057)	(10, 9433)	

# Chiricahua
# (2, 3346)	(3, 5019)	(4, 6149)	(5, 7008)	(6, 7702)	(7, 8286)	(8, 8791)	

# Mogollon
# (2, 3269)	(3, 4904)	(4, 6030)	(5, 6903)	(6, 7620)	(7, 8232)	(8, 8767)	(9, 9243)	(10, 9672)	(11, 10064)	(12, 10424)	(13, 10758)	(14, 11069)	(15, 11360)	(16, 11634)	(17, 11893)	(18, 12139)	(19, 12372)	(20, 12594)	(21, 12807)	(22, 13010)	(23, 13206)	(24, 13393)	

# Outgroup
# (2, 2582)	(3, 3873)	(4, 4701)	(5, 5299)	(6, 5760)	(7, 6133)	(8, 6445)	(9, 6711)	(10, 6942)	

# 2,8,8,20,10

$easySFS \
-i filtered_thinned.recode.vcf.gz \
-p pop.map \
--proj 2,8,8,20,10 \
-o sfs_output \
--prefix D_grantii \
-f \
--dtype int


ls -lh sfs_output/fastsimcoal2/

# Create organized directory for fastsimcoal2
# Copy SFS to each model directory

# ─────────────────────────────────────────────────
# Time
# ▲
# │  Ancestral population
# │       |
# │  ─────────────────────────────
# │  |                           |
# │  Outgroup               Sky Islands ancestor
# │  (Utah)                      |
# │                    ─────────────────────
# │                    |                   |
# │               Mogollon           Lemmon-Graham
# │               lineage             lineage
# │                                     |
# │                              ───────────────
# │                              |             |
# │                           Lemmon         Graham
# │                                           +
# │                                       Chiricahua?
# ▼ Present

# Population order in fastsimcoal2:
# 0=Graham, 1=Lemmon, 2=Chiricahua, 3=Mogollon, 4=Outgroup


fsc28=/scratch/user/schien/Software/fsc28_linux64/fsc28
$fsc28

# TEST 

$FSC \
--tplfile M0_isolation.tpl --estfile M0_isolation.est \
--multiSFS -E 1 --numsims 100 --maxiter 10 --minSFSCount 1 --quiet

$FSC \
--tplfile M1_steppingstone.tpl --estfile M1_steppingstone.est \
--multiSFS -E 1 --numsims 100 --maxiter 10 --minSFSCount 1 --quiet

$FSC \
--tplfile M2_directflow.tpl --estfile M2_directflow.est \
--multiSFS -E 1 --numsims 100 --maxiter 10 --minSFSCount 1 --quiet
