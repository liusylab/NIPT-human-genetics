#!/bin/sh
#SBATCH -J Baoan_PCA10
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=30
#SBATCH --time=10-05:05:30
#SBATCH --mem=100G
#SBATCH --partition=bigmem
#SBATCH --output=slurm_Baoan_unimputed_PCA10
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
plink2 --maf 0.05 --vcf GLs.merged.chr1-22.0.001.vcf.gz dosage=DS --pca 10 --out NIPT_rm_dup_merged.chr1-22.0.05_pca10 --threads 30

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 10s" 
sleep 10
echo "process end at : "
date
