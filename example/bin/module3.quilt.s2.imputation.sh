#!/bin/bash

ref=/share/home/lsy_liusiyang/20220708_Alignment/hg38/Homo_sapiens_assembly38.fasta
bamlist=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/batch1_final_files/bam.list
quilt_outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/quilt_output
map_outdir=$quilt_outdir/genetic_map_file
reference_path=$quilt_outdir/reference_file

QUILT=/share/home/lsy_liuyanhong/software/QUILT/QUILT.R
regionStart=1
regionEnd=64444167

for i in 20 #to enumerate chromosomes, use for i in {1..22} X
do

reference_haplotype_file=${reference_path}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.biallelic.snp.maf0.001.unique_pos.hap.gz
reference_legend_file=${reference_path}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.biallelic.snp.maf0.001.unique_pos.legend.gz
genetic_map_file=${map_outdir}/CHS_chr$i.txt.gz

echo "#!/bin/sh
#SBATCH -J make_quilt.chr20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=07-00:00:00
#SBATCH --mem=10240
#SBATCH -p cpu
#SBATCH --output=quilt.s2.chr$i.slurm.out
echo \"process will start at :\"
date" >quilt.s2.imputation.chr${i}.tmp.sh

echo "
$QUILT \
--outputdir=$quilt_outdir \
--reference=$ref \
--buffer=250000 \
--chr=chr$i \
--regionStart=$regionStart \
--regionEnd=$regionEnd \
--nGibbsSamples=7 \
--n_seek_its=3 \
--bamlist=$bamlist \
--reference_haplotype_file=$reference_haplotype_file \
--reference_legend_file=$reference_legend_file \
--genetic_map_file=$genetic_map_file \
--reference_populations=CHB,CHS,CDX \
--nGen=1240 \
--save_prepared_reference=TRUE
echo \"process end at : \"
date" >>quilt.s2.imputation.chr${i}.tmp.sh

done
