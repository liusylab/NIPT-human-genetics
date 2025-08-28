#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06-00:00:00
#SBATCH --mem=1G
#SBATCH --job-name=emu_pca
#SBATCH --output=emu_pca.slurm.out

emu=/share/home/lsy_liusiyang/software/emu
imputed_vcf=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/glimpse_output/imputed_file_merged/glimpse.chr20_imputed.vcf.gz
prefix=`basename $imputed_vcf|sed -e 's/\.vcf\.gz//'`
emu_outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/emu_output
mkdir -p $emu_outdir

$emu -m -p $imputed_vcf -e 10 -t 30 -o $imputed_vcf.emu10

