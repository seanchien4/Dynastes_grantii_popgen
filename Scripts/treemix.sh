#!/bin/bash
# TreeMix

# create population map
vcf=/scratch/user/schien/D_grantii/Population_Genetics/vcf/final_geno_filter.g.vcf.gz

ml GCC/12.3.0 PLINK/2.00a3.7
# prune LD and create a new vcf file
plink --vcf ../Dsuite/Dsuit_geno.g.vcf.gz --allow-extra-chr \
--keep keep.txt \
--indep-pairwise 50 10 0.2 --out ld_pruned

plink --vcf ../Dsuite/Dsuit_geno.g.vcf.gz --allow-extra-chr \
--keep keep.txt \
-extract ld_pruned.prune.in --recode vcf --out filtered

gzip filtered.vcf

# modify the file so there is ID for each position
ml purge 
ml GCC/13.2.0 BCFtools/1.19
bcftools annotate --set-id '%CHROM:%POS' -o filtered_with_IDs.vcf.gz -Oz filtered.vcf.gz

ml purge
ml GCC/12.3.0 PLINK/2.00a3.7
plink --vcf filtered_with_IDs.vcf.gz --allow-extra-chr --freq --missing --within population.txt --out plink_output

gzip plink_output.frq.strat

# convert vcf to treemix format 
# virrtual env
ml GCCcore/12.3.0 Python/2.7.18
# pip install pandas
python2.7 /scratch/user/schien/Software/plink2treemix.py plink_output.frq.strat.gz treemix_input.gz

# run treemix
# choose migration edge 
treemix -i treemix_input.gz -o test

for i in {0..5}
do
 treemix -i treemix_input.gz -m $i -root Outgroup -o migration.$i -bootstrap -k 500 -noss > treemix_${i}_log &
done
# checking likelihood
grep "ln(likelihood):" treemix_?_log



for i in {0..5}
do
echo -n "m=$i: "
grep "ln(likelihood):" treemix_${i}_log | tail -1
done

# Best m = 1
# - Likelihood improved from m=0 → m=1
# - No further improvement after m=1
# - Adding more migrations = overfitting


