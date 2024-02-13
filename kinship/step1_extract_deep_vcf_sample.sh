#!/bin/sh
#SBATCH -J baoan_extract_vcf_sample_10K
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=07-00:00:00
#SBATCH --mem=50G
#SBATCH -p bigmem
#SBATCH --output=slurm_baoan_extract_vcf_sample_10K
workdir=/kinship
bcftools view --force-samples -S identical_samples_kinship_SampleID.txt $workdir/unimpute_vcf_file/baoan_70608_GLs.merged.chr1-22.0.001.vcf.gz -Oz -o $workdir/sample_extract_30n/extract.70608_GLs.merged.chr1-22.0.001.vcf.gz
tabix $workdir/sample_extract_30n/extract.70608_GLs.merged.chr1-22.0.001.vcf.gz
echo "OK"
echo "process end at : "
date
