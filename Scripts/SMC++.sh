#!/bin/bash
# SMC++
# filtering VCF file 
module purge
module load GCC/12.2.0 GCC/12.3.0 VCFtools/0.1.16
unfilter_vcf=/scratch/user/schien/D_grantii/Population_Genetics/vcf/38cohort_geno.g.vcf.gz
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
--recode --stdout | bgzip -c > final_geno_smc.g.vcf.gz

tabix final_geno_smc.g.vcf.gz


# 

# knots=20

echo "#*********************#"
echo "#      Mt. Lemmon     #"
echo "#*********************#"

for ind in $(cat lemmon.list);
do 
	for contig in $(cat chrom.names);
do
	echo $contig
        singularity exec $smcpp smc++ vcf2smc --cores $SLURM_CPUS_PER_TASK \
        $vcf out/ML_${ind}_$contig.smc.gz $contig \
        ML:${ind}
done
done

sbatch smc++_Lemmon_estimate.sh


### Estimate ###

# knots=20
t=1
t1=1000000
k=${SLURM_ARRAY_TASK_ID}
singularity exec $smcpp smc++ estimate --cores $SLURM_CPUS_PER_TASK --timepoints $t $t1  -o analysis/ML_ 5.8e-9 out/ML_*.smc.gz
singularity exec $smcpp smc++ plot ML_$k.pdf -c analysis/ML_/model.final.json
rm ML_$k.pdf
mv ML_$k.csv csv/



