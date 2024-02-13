#!/bin/sh
work_path=glimpse_pipline/ref_1KGP
for i in {1..22} X
do
mkdir -p ${work_path}/bin/tasks
mkdir -p ${work_path}/imputed_file
MAP=/software/GLIMPSE-1.1.1/maps/genetic_maps.b38/chr${i}.b38.gmap.gz
VCF=$work_path/GL_file_merged/high_dep_100.chr${i}.vcf.gz
REF=${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz
while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT=${work_path}/imputed_file/high_dep_100.chr${i}.${ID}.imputed.vcf
#prepare tasks
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=07-00:00:00
#SBATCH --mem=20480
#SBATCH --partition=cpu
#SBATCH --job-name=phase_chr${i}_${ID}
#SBATCH --output=${work_path}/bin/tasks/slurm_phase_chr${i}_${ID}.sh
echo \"process will start at :\"
date
/software/GLIMPSE-1.1.1/phase/bin/GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
/software/bin/bgzip ${OUT}
/software/bcftools-1.15/bcftools index -f ${OUT}.gz
echo \"process end at : \"
date" > ${work_path}/bin/tasks/sbatch_phase_chr${i}_${ID}.sh
sbatch ${work_path}/bin/tasks/sbatch_phase_chr${i}_${ID}.sh
done < $work_path/chunks.G10K.chr${i}.txt
done

