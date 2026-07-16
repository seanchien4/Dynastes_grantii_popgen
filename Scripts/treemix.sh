
treemix_input=/scratch/user/schien/D_grantii/Population_Genetics/vcf/final_geno_filter.g.vcf.gz

# LD pruning
module purge
module load GCC/12.3.0 PLINK/2.00a3.7
# plink pruning
plink --vcf $treemix_input --double-id --allow-extra-chr \
 --set-missing-var-ids @:# \
 --indep-pairwise 50 10 0.2 --out treemix

plink --vcf $treemix_input --double-id --allow-extra-chr \
 --set-missing-var-ids @:# \
 --extract treemix.prune.in \
 --make-bed --out treemix_pruned

awk '{print $1, $1, $2}' ../Dsuite/pop.map > treemix.clust

plink --bfile treemix_pruned \
  --freq --missing \
  --within treemix.clust \
  --allow-extra-chr \
  --out treemix_freq

gzip treemix_freq.frq.strat

zcat treemix.frq.gz | awk 'NR==1 || $0 !~ /(^| )0,0( |$)/' | gzip > treemix_nomiss.frq.gz
# Verify
zcat treemix_nomiss.frq.gz | wc -l
# should be ~92,808 (92,807 SNPs + 1 header)

zcat treemix_nomiss.frq.gz | head -3

# test 
treemix -i treemix_nomiss.frq.gz -m 0 -root Outgroup -k 500 -o test_m0


source activate treemix
mkdir treemix_out && cd treemix_out

for i in {0..5}
do
  for rep in {1..50}
  do
    treemix -i treemix_nomiss.frq.gz -m $i -root Outgroup \
            -k 500 -bootstrap -seed $RANDOM \
            -o m${i}_rep${rep} > m${i}_rep${rep}.log &
  done
  wait
done
for i in {0..5}
do
  echo -n "m=$i: "
  grep "ln(likelihood):" m${i}_rep*.log | tail -1
done



# Best m = 1
# - Likelihood improved from m=0 → m=1
# - No further improvement after m=1
# - Adding more migrations = overfitting