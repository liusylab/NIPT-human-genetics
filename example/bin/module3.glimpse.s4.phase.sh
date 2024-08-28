#!/bin/bash

bgzip=/share/home/lsy_liusiyang/software/tabix-0.2.6/bgzip
bcftools=/share/home/lsy_liusiyang/software/bcftools-1.15/bcftools
GLIMPSE_phase=/share/home/lsy_liusiyang/software/GLIMPSE-1.1.1/phase/bin/GLIMPSE_phase

glimpse_outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/glimpse_output
map_path=/share/home/lsy_liusiyang/software/GLIMPSE-1.1.1/maps/genetic_maps.b38
reference_path=$glimpse_outdir/reference_file

for i in 20 #to enumerate chromosomes, use for i in {1..22} X
do

mkdir -p ${glimpse_outdir}/imputed_file
MAP=${map_path}/chr${i}.b38.gmap.gz
REF=${reference_path}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.biallelic.snp.maf0.001.vcf.gz
CHUNKFILE=$reference_path/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.chunks.txt
VCF=${glimpse_outdir}/GL_file_merged/glimpse.chr${i}.vcf.gz

cat ${CHUNKFILE}|while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT=${glimpse_outdir}/imputed_file/glimpse.chr${i}.${ID}.imputed.vcf

#prepare tasks
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=07-00:00:00
#SBATCH --mem=20480
#SBATCH --partition=cpu
#SBATCH --job-name=glimpse.s4.phase_chr${i}_${ID}
#SBATCH --output=glimpse.s4.phase_chr${i}_${ID}.slurm.out
echo \"process will start at :\"
date" > glimpse.s4.phase.chr${i}_${ID}.tmp.sh

echo "$GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
$bgzip ${OUT}
$bcftools index -f ${OUT}.gz
echo \"process end at : \"
date" >>glimpse.s4.phase.chr${i}_${ID}.tmp.sh
done 

done

