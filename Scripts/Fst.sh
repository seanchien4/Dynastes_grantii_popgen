#!/bin/bash
# Fst
vcf=/scratch/user/schien/D_grantii/Population_Genetics/vcf/final_geno_filter.g.vcf.gz
# filter the vcf file only keep the individuals that are going to be used in the analysis
ml purge 
ml  GCC/12.2.0 GCC/12.3.0 GCC/13.2.0 VCFtools/0.1.16
vcftools \
--gzvcf $vcf \
--keep keep.txt \
--recode \
--recode-INFO-all \
--out Fst

vcftools --gzvcf Fst.recode.vcf \
         --weir-fst-pop utah_samples.txt \
         --weir-fst-pop mogollon_samples.txt \
         --out Utah_vs_Mogollon &

vcftools --gzvcf Fst.recode.vcf \
         --weir-fst-pop utah_samples.txt \
         --weir-fst-pop lemmon_samples.txt \
         --out Utah_vs_Lemmon &

vcftools --gzvcf Fst.recode.vcf \
         --weir-fst-pop utah_samples.txt \
         --weir-fst-pop portal_samples.txt \
         --out Utah_vs_Portal &

vcftools --gzvcf Fst.recode.vcf \
         --weir-fst-pop lemmon_samples.txt \
         --weir-fst-pop mogollon_samples.txt \
         --out Lemmon_vs_Mogollon &

vcftools --gzvcf Fst.recode.vcf \
         --weir-fst-pop portal_samples.txt \
         --weir-fst-pop mogollon_samples.txt \
         --out Portal_vs_Mogollon &

vcftools --gzvcf Fst.recode.vcf \
         --weir-fst-pop portal_samples.txt \
         --weir-fst-pop lemmon_samples.txt \
         --out Portal_vs_Lemmon &

awk '$3 ~ /^-?[0-9]*\.?[0-9]+$/ && $3 > 0 { sum += $3; n++ } END { if (n) print sum / n }' Utah_vs_Mogollon.weir.fst 
awk '$3 ~ /^-?[0-9]*\.?[0-9]+$/ && $3 > 0 { sum += $3; n++ } END { if (n) print sum / n }' Lemmon_vs_Mogollon.weir.fst
awk '$3 ~ /^-?[0-9]*\.?[0-9]+$/ && $3 > 0 { sum += $3; n++ } END { if (n) print sum / n }' Portal_vs_Lemmon.weir.fst
awk '$3 ~ /^-?[0-9]*\.?[0-9]+$/ && $3 > 0 { sum += $3; n++ } END { if (n) print sum / n }' Portal_vs_Mogollon.weir.fst
awk '$3 ~ /^-?[0-9]*\.?[0-9]+$/ && $3 > 0 { sum += $3; n++ } END { if (n) print sum / n }' Utah_vs_Lemmon.weir.fst
awk '$3 ~ /^-?[0-9]*\.?[0-9]+$/ && $3 > 0 { sum += $3; n++ } END { if (n) print sum / n }' Utah_vs_Mogollon.weir.fst
awk '$3 ~ /^-?[0-9]*\.?[0-9]+$/ && $3 > 0 { sum += $3; n++ } END { if (n) print sum / n }' Utah_vs_Portal.weir.fst

