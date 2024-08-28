#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06-00:00:00
#SBATCH --mem=1G
#SBATCH --job-name=plink_kinship
#SBATCH --output=plink.kinship.slurm.out

plink2=/share/home/lsy_liusiyang/software/plink2.0/plink2
imputed_vcf=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/glimpse_output/imputed_file_merged/glimpse.chr20_imputed.vcf.gz
prefix=`basename $imputed_vcf|sed -e 's/\.vcf\.gz//'`
kinship_outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/kinship_output
mkdir -p $kinship_outdir

$plink2 --vcf $imputed_vcf --make-king-table --out $kinship_outdir/$prefix --threads 24
