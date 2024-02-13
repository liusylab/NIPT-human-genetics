#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=24
#SBATCH --time=06-00:00:00
#SBATCH --mem=50G
#SBATCH --partition=bigmem
#SBATCH --job-name=baoan_kinship
#SBATCH --output=/kinship
workdir=/kinship
plink2 --vcf $workdir/sample_extract_30n/extract.70608_GLs.merged.chr1-22.0.001.vcf.gz --make-king-table --out $workdir/plink2_kinship/70608_GLs.merged.chr1-22.0.001_maf --threads 24
