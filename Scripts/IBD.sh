#!/bin/bash

# IBD test
# Mantel

# First, create coords.txt

ml GCC/12.3.0 PLINK/2.00a3.7

# ALL
plink --vcf ../../vcf/final_geno_filter.g.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --keep all.list \
--make-bed --out pop_all
plink --vcf ../../vcf/final_geno_filter.g.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# --keep all.list \
--extract pop_all.prune.in --distance square \
--out all_dist

# Payson
plink --vcf ../../vcf/final_geno_filter.g.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --keep pop1_Payson.list \
--make-bed --out pop1_payson
plink --vcf ../../vcf/final_geno_filter.g.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# --keep pop1_Payson.list \
--extract pop1_payson.prune.in --distance square \
--out pop1_payson_dist

# Reserve
plink --vcf ../../vcf/final_geno_filter.g.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --keep pop1_Reserve.list \
--make-bed --out pop1_reserve
plink --vcf ../../vcf/final_geno_filter.g.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# --keep pop1_Reserve.list \
--extract pop1_reserve.prune.in --distance square \
--out pop1_reserve_dist

# Utah
plink --vcf ../../vcf/final_geno_filter.g.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --keep utah.list  \
--make-bed --out utah
plink --vcf ../../vcf/final_geno_filter.g.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# --keep utah.list \
--extract utah.prune.in --distance square \
--out utah_dist


# --keep coords_payson.txt
# plot and do analysis in R


# Cryptic relatedness test
plink --bfile popA_only --genome --allow-extra-chr --out relatedness




# Pairwise Fst test
plink --bfile popA_only --allow-extra-chr --recode vcf --out popA_only_for_vcftools

ml purge
ml GCC/12.2.0 GCC/12.3.0 GCC/13.2.0 VCFtools/0.1.16
vcftools --vcf popA_only_for_vcftools.vcf \
         --weir-fst-pop Mt_Lemmon.txt \
         --weir-fst-pop Payson.txt \
         --out Lemmon_vs_Payson



plink --bfile popA_only --allow-extra-chr \
--fst --within populations.clst --out non_utah_pairwise_fst
