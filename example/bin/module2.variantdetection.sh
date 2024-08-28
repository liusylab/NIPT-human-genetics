#!/bin/bash

basevarpip=../../../NIPT-human-genetics/basevar/0.6.1.1/scripts/create_pipeline.py
bamlist=../../../NIPT-human-genetics/example/outdir/batch1_final_files/bam.list


ref=/share/home/lsy_liusiyang/20220708_Alignment/hg38/Homo_sapiens_assembly38.fasta
ref_fai=/share/home/lsy_liusiyang/20220708_Alignment/hg38/Homo_sapiens_assembly38.fasta.fai
bgzip=/share/home/lsy_liusiyang/software/tabix-0.2.6/bgzip
tabix=/share/home/lsy_liusiyang/software/tabix-0.2.6/tabix
bcftools=/share/home/lsy_liusiyang/software/bcftools-1.15/bcftools

outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/basevar_output
outprefix=NIPT_basevar_chr20

mkdir -p $outdir

echo "#!/bin/sh
#SBATCH -J NIPT-human-genetics_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=00-05:05:30
#SBATCH --mem=18432
#SBATCH --partition=cpu" > basevar.sh
python $basevarpip -R $ref --ref_fai $ref_fai -c chr20 --delta 1000000 --nCPU 10 -b $bgzip -t $tabix -L $bamlist -o $outdir >> basevar.tmp.sh

awk '{print $33}' basevar.tmp.sh > $outdir/$outprefix.vcf.list

echo "#!/bin/sh
#SBATCH -J NIPT-human-genetics_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=00-05:05:30
#SBATCH --mem=18432
#SBATCH --partition=cpu" > basevar.merge.tmp.sh
awk '{s=s" "$1;}END{print "time '$bcftools' concat --threads 12 -a --rm-dups all -O z -o '$outdir'/'$outprefix'.vcf.gz"s" && '$tabix' -f -p vcf '$outdir'/'$outprefix'.vcf.gz && echo \"** merge vcf done **\""}' $outprefix.vcf.list >> basevar.merge.tmp.sh



