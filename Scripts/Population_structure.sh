#!/bin/bash
# Sean Chien 

# There are too many SNPs to do analyses for setting filter parameters
# So I am going to randomly pick 10,000 SNPs
# there are 4 parameters need to be considered
# 1. quality 
# 2. depth 
# 3. Minor allele frequency
# 4. missing data
module load GCC/12.3.0 BCFtools/1.18
# unfilter variants
bcftools view -H 25_cohort_genotype.g.vcf.gz | wc -l
# 35346468

# check sample ID
bcftools query -l 25_cohort_genotype.g.vcf.gz
bcftools query -l 25_cohort_genotype.g.vcf.gz | grep DGS01 > exclude.txt
vcftools --gzvcf 25_cohort_genotype_filtered.vcf.gz --remove exclude.txt --recode --stdout | gzip >  25_cohort_genotype_filtered_ExUT.vcf.gz
vcftools --gzvcf 25_cohort_genotype_filtered.vcf.gz --remove exclude.txt --recode --stdout | gzip >  25_cohort_genotype_filtered_ExMexico.vcf.gz


# randomly pick 1000,000
#!/bin/bash
#SBATCH --job-name=bwa              # job name
#SBATCH --time=0-20:00:00           # max job run time dd-hh:mm:ssll
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=20          # CPUs (threads) per command
#SBATCH --mem=360G                  # total memory per node
#SBATCH --output=bwa_stdout         # save stdout to file
#SBATCH --error=bwa_stderr          # save stderr to file

# header
bcftools view -h 25_cohort.g.vcf.gz > header.txt
bcftools view -H 25_cohort.g.vcf.gz | shuf -n 1000000 > random_loci.vcf
# combine
cat header.txt random_loci.vcf > 100000_loci.vcf
rm header.txt random_loci.vcf
gzip 100000_loci.vcf

module purge
module load GCC/11.2.0 GCC/9.3.0 iccifort/2019.5.281 VCFtools/0.1.16
# random loci statistic
# allele frequency
vcftools --gzvcf 100000_loci.vcf.gz --freq2 --out 25_freq --max-alleles 2 &

# mean dp for each individual
vcftools --gzvcf 100000_loci.vcf.gz --depth --out 25_dp &

# mean depth per site
vcftools --gzvcf 100000_loci.vcf.gz --site-mean-depth --out 25_dp_site &

# site quality
vcftools --gzvcf 100000_loci.vcf.gz --site-quality --out 25_qc_site &

# missing data
# Generates a file reporting the missingness on a per-individual basis
vcftools --gzvcf 100000_loci.vcf.gz --missing-indv --out 25_miss &

# missing data site
vcftools --gzvcf 100000_loci.vcf.gz --missing-site --out 25_miss_site &

# heterozygosity
vcftools --gzvcf 100000_loci.vcf.gz --het --out 25_heter &

# go to R to do statistic
# filtering vcf files

vcftools --gzvcf 25_cohort_genotype.g.vcf.gz --missing-indv

# DgSR1 has 40% missing so remove this one
vcftools --gzvcf 25_cohort_genotype.g.vcf.gz \
--remove-indels --maf 0.1 --max-missing 0.9 --minQ 20 \
--min-meanDP 10 --max-meanDP 30 \
--remove-indv DgSR1 \
--minDP 10 --maxDP 30 --recode --stdout | gzip > 25_cohort_genotype_filtered_exSR1.vcf.gz


# PCA 
# One of the assumptions of PCA is the data is independent
# so we need to consider and filter out LD
# find the linkage disequilibrium decays to the genome wide background 
# PopLDdecays
# calculate LD decay
./PopLDdecay/bin/PopLDdecay -InVCF ../vcf_handle/25_cohort_genotype_filtered.vcf.gz  -OutStat Dynastes_LDdecay 

module load GCC/12.3.0 R/4.3.2
# draw the figure
perl PopLDdecay/bin/Plot_OnePop.pl -inFile Dynastes_LDdecay.stat.gz  -output Fig

# set LD to 10kb and filter out 
# plink filtering LD
module purge
module load GCC/12.3.0 PLINK/2.00a3.7

# plink pruning
plink --vcf 36cohort_geno_filter_PCA.g.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 10 10 0.1 --out 36cohort_LD

# plink PCA
plink --vcf 36cohort_geno_filter_PCA.g.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract 36cohort_LD.prune.in \
--make-bed --pca --out 36cohort_PCA

# ADMIXTURE

cp ../PCA/final.b*
cp ../PCA/final.fam .
cp ../PCA/final.nosex .

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the fi$
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




