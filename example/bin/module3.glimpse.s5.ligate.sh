#!/bin/bash

reference_path=/share/home/lsy_liuyanhong/20230511_QUILT/quilt_110_1KGP/input/1KGP_legend_haplotype

GLIMPSE_ligate=/share/home/lsy_liusiyang/software/GLIMPSE-1.1.1/ligate/bin/GLIMPSE_ligate
bgzip=/share/home/lsy_liusiyang/software/tabix-0.2.6/bgzip
bcftools=/share/home/lsy_liusiyang/software/bcftools-1.15/bcftools

glimpse_outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/glimpse_output
mkdir -p ${glimpse_outdir}/imputed_file_merged

#prepare tasks
for i in 20 #can enumerate chromosomes by {1..22} X
do
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=07-00:00:00
#SBATCH --mem=20480
#SBATCH --partition=cpu
#SBATCH --job-name=glimpse.s5.ligate.chr$i
#SBATCH --output=glimpse.s5.ligate.chr$i.slurm.out
echo \"process will start at :\"
date" >glimpse.s5.ligate.chr$i.tmp.sh

echo "ls ${glimpse_outdir}/imputed_file/glimpse.chr${i}.*.imputed.vcf.gz > ${glimpse_outdir}/imputed_file/glimpse.chr${i}_imputed_list.txt
/share/home/lsy_guyuqin/software/GLIMPSE-1.1.1/ligate/bin/GLIMPSE_ligate --input ${glimpse_outdir}/imputed_file/glimpse.chr${i}_imputed_list.txt --output ${glimpse_outdir}/imputed_file_merged/glimpse.chr${i}_imputed.vcf
$bgzip ${glimpse_outdir}/imputed_file_merged/glimpse.chr${i}_imputed.vcf
$bcftools index -f ${glimpse_outdir}/imputed_file_merged/glimpse.chr${i}_imputed.vcf.gz
echo \"process end at : \"
date" >> glimpse.s5.ligate.chr$i.tmp.sh
done
