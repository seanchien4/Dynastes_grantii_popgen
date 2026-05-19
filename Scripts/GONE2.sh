#!/bin/bash
# GONE2

# filter out small scaffold
vcf="../smc/final_geno_smc.g.vcf.gz"


# bed file
awk '{OFS="\t"; print $1, "0", $2}' ../ref/GCA_029619325.2_ASM2961932v2_genomic.fna.fai > all_scaffolds.bed
awk '$2 > 20000000 {print $1}' ../ref/GCA_029619325.2_ASM2961932v2_genomic.fna.fai  > large_scaffolds.txt
grep -wFf large_scaffolds.txt all_scaffolds.bed > large_scaffolds.bed

# filter vcf file
# population 
ml purge 
ml GCC/12.2.0 GCC/12.3.0 GCC/13.2.0 VCFtools/0.1.16
vcftools --gzvcf $vcf \
--maf 0.01 --keep mogollon.list --max-missing 1 --minQ 20 \
--min-meanDP 10 --max-meanDP 30 --minDP 10 --maxDP 30 --recode --stdout | bgzip -c > Mogollon_for_GONE.g.vcf.gz &

vcftools --gzvcf $vcf \
--maf 0.01 --keep lemmon.list --max-missing 1 --minQ 20 \
--min-meanDP 10 --max-meanDP 30 --minDP 10 --maxDP 30 --recode --stdout | bgzip -c > Lemmon_for_GONE.g.vcf.gz &

vcftools --gzvcf $vcf \
--maf 0.01 --keep portal.list --max-missing 1 --minQ 20 \
--min-meanDP 10 --max-meanDP 30 --minDP 10 --maxDP 30 --recode --stdout | bgzip -c > Portal_for_GONE.g.vcf.gz &

vcftools --gzvcf $vcf \
--maf 0.01 --keep utah.list --max-missing 1 --minQ 20 \
--min-meanDP 10 --max-meanDP 30 --minDP 10 --maxDP 30 --recode --stdout | bgzip -c > Utah_for_GONE.g.vcf.gz &

tabix -p vcf Mogollon_for_GONE.g.vcf.gz &
tabix -p vcf Lemmon_for_GONE.g.vcf.gz &
tabix -p vcf Portal_for_GONE.g.vcf.gz &
tabix -p vcf Utah_for_GONE.g.vcf.gz &

ml purge 
ml GCC/13.2.0 BCFtools/1.19
bcftools view  --regions-file large_scaffolds.bed Mogollon_for_GONE.g.vcf.gz -o Mogollon_geno_large_scaffolds.vcf.gz &
bcftools view  --regions-file large_scaffolds.bed Lemmon_for_GONE.g.vcf.gz -o Lemmon_geno_large_scaffolds.vcf.gz &
bcftools view  --regions-file large_scaffolds.bed Portal_for_GONE.g.vcf.gz -o Portal_geno_large_scaffolds.vcf.gz &
bcftools view  --regions-file large_scaffolds.bed Utah_for_GONE.g.vcf.gz -o Utah_geno_large_scaffolds.vcf.gz &

# Plink for input files for GONE2
ml purge
ml GCC/12.3.0 PLINK/2.00a3.7
plink --vcf Mogollon_geno_large_scaffolds.vcf.gz --recode --allow-extra-chr --out GONE2_Mogollon --geno 0.0 
plink --vcf Lemmon_geno_large_scaffolds.vcf.gz --recode --allow-extra-chr --out GONE2_Lemmon --geno 0.0
plink --vcf Portal_geno_large_scaffolds.vcf.gz --recode --allow-extra-chr --out GONE2_Portal --geno 0.0
plink --vcf Utah_geno_large_scaffolds.vcf.gz --recode --allow-extra-chr --out GONE2_Utah --geno 0.0

ml GCC/14.2.0 OpenMPI/5.0.7
make MAXLOCI=5000000 gone

GONE2=/scratch/user/schien/Software/GONE2/gone2
r=2.48
cpu=8
$GONE2 GONE2_Mogollon.ped -r $r -t $cpu
$GONE2 GONE2_Lemmon.ped -r $r -t $cpu
$GONE2 GONE2_Portal.ped -r $r -t $cpu
$GONE2 GONE2_Utah.ped -r $r -t $cpu

# BS 100 for each population
#!/bin/bash
#SBATCH --job-name=GONE2_Dynastes_BS
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=schien@tamu.edu
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --array=1-100
#SBATCH --account=132704493814

ml purge 
ml GCC/14.2.0 OpenMPI/5.0.7
GONE2=/scratch/user/schien/Software/GONE2/gone2
currentNe2=/scratch/user/schien/Software/currentNe2/currentne2
i=${SLURM_ARRAY_TASK_ID}
r=2.48
$GONE2 GONE2_Mogollon.ped -r $r -t $SLURM_CPUS_PER_TASK -o Mogollon_$i
$currentNe2 GONE2_Mogollon.ped -r $r -t 8

$GONE2 GONE2_Lemmon.ped -r $r -t $SLURM_CPUS_PER_TASK -o Lemmon_$i
$currentNe2 GONE2_Lemmon.ped -r $r -t 8

$GONE2 GONE2_Portal.ped -r $r -t $SLURM_CPUS_PER_TASK -o Portal_$i
$currentNe2 GONE2_Portal.ped -r $r -t 8

$GONE2 GONE2_Utah.ped -r $r -t $SLURM_CPUS_PER_TASK -o Utah_$i
$currentNe2 GONE2_Utah.ped -r $r -t 8

mv Mogollon_*_GONE2_* GONE_Ne/Mogollon &
mv Lemmon_*_GONE2_* GONE_Ne/Lemmon &
mv Portal_*_GONE2_* GONE_Ne/Portal &
mv Utah_*_GONE2_* GONE_Ne/Utah &


currentNe2=/scratch/user/schien/Software/currentNe2/currentne2
$currentNe2 GONE2_Mogollon.ped -r $r -t 8 -o test
i=${SLURM_ARRAY_TASK_ID}
r=2.48
$GONE2 GONE2_Mogollon.ped -r $r -t $SLURM_CPUS_PER_TASK -o Mogollon_$i
$currentNe2 GONE2_Mogollon.ped -r $r -t $SLURM_CPUS_PER_TASK -o Mogollon_$i

$GONE2 GONE2_Lemmon.ped -r $r -t $SLURM_CPUS_PER_TASK -o Lemmon_$i
$currentNe2 GONE2_Lemmon.ped -r $r -t $SLURM_CPUS_PER_TASK -o Lemmon_$i

$GONE2 GONE2_Portal.ped -r $r -t $SLURM_CPUS_PER_TASK -o Portal_$i
$currentNe2 GONE2_Portal.ped -r $r -t $SLURM_CPUS_PER_TASK -o Portal_$i

$GONE2 GONE2_Utah.ped -r $r -t $SLURM_CPUS_PER_TASK -o Utah_$i
$currentNe2 GONE2_Utah.ped -r $r -t $SLURM_CPUS_PER_TASK -o Utah_$i

