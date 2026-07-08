#!/bin/bash
# Population structure
#######
# PCA #
#######
module purge
module load GCC/12.2.0 GCC/12.3.0 VCFtools/0.1.16
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
--recode --stdout | bgzip -c > final_geno_filter.g.vcf.gz
tabix final_geno_filter.g.vcf.gz 

# One of the assumptions of PCA is the data is independent
# so we need to consider and filter out LD
# find the linkage disequilibrium decays to the genome wide background 
# PopLDdecays
# calculate LD decay
/scratch/user/schien/D_grantii/Population_Genetics/PCA/PopLDdecay/PopLDdecay/bin/PopLDdecay \
-InVCF ../../vcf/final_geno_filter.g.vcf.gz -OutStat Dynastes_LDdecay
module load GCC/12.3.0 R/4.3.2
# draw the figure
perl PopLDdecay/bin/Plot_OnePop.pl -inFile Dynastes_LDdecay.stat.gz  -output Fig

# set LD to 10kb and filter out 
# plink filtering LD
module purge
module load GCC/12.3.0 PLINK/2.00a3.7
# plink pruning
plink --vcf ../vcf/final_geno_filter.g.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.2 --out final

# plink PCA
plink --vcf ../vcf/final_geno_filter.g.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract final.prune.in \
--make-bed --pca --out final


# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' final.bim > final.bim.tmp
mv final.bim.tmp final.bim

admixture --cv final.bed 2 > log2.out &
for i in {3..9}
do
admixture --cv final.bed $i > log${i}.out &
done



awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > admixture.cv.error
awk '{split($1,name,"."); print $1,name[2]}' final.nosex > admixture.list

for i in {2..9}
do
paste admixture.list final.$i.Q > genetic_prop_$i.txt
done
