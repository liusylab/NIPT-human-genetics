#!/bin/bash

###
#Users should substitute the following paths with their respective local directories.
###

#Inputfile
bamlist=../outdir/batch1_final_files/bam.list

#Pipelines, softwares and tools
basevarpip=../../BaseVar2/scripts/create_pipeline.py
tabix=/path/to/tabix
bcftools=/path/to/bcftools

#Resources
ref=/path/to/hg38/Homo_sapiens_assembly38.fasta
ref_fai=/path/to/hg38/Homo_sapiens_assembly38.fasta.fai

#Output
outdir=../outdir/basevar_output
outprefix=NIPT_basevar_chr20

mkdir -p $outdir

echo "#!/bin/sh
#SBATCH -J NIPT-human-genetics_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=00-05:05:30
#SBATCH --mem=18432
#SBATCH --partition=cpu" > basevar.sh
python $basevarpip -Q 20 -q 30 -f $ref --ref_fai $ref_fai -c chr20 --delta 5000000 -t 24 -L $bamlist -o $outdir >> basevar.sh

awk '{print $33}' basevar.sh > $outdir/$outprefix.vcf.list

echo "#!/bin/sh
#SBATCH -J NIPT-human-genetics_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=00-05:05:30
#SBATCH --mem=18432
#SBATCH --partition=cpu" > basevar.merge.sh
awk '{s=s" "$1;}END{print "time '$bcftools' concat --threads 12 -a --rm-dups all -O z -o '$outdir'/'$outprefix'.vcf.gz"s" && '$tabix' -f -p vcf '$outdir'/'$outprefix'.vcf.gz && echo \"** merge vcf done **\""}' $outprefix.vcf.list >> basevar.merge.sh



