#!/bin/sh
work_path=glimpse_pipline/ref_1KGP
mkdir -p ${work_path}/bin/tasks
mkdir -p ${work_path}/imputed_file_merged
#prepare tasks
for i in {1..22}
do
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=07-00:00:00
#SBATCH --mem=20480
#SBATCH --partition=cpu
#SBATCH --job-name=ligate
#SBATCH --output=${work_path}/bin/tasks/slurm_ligate_chr${i}.sh
echo \"process will start at :\"
date
ls ${work_path}/imputed_file/high_dep_100.chr${i}.*.imputed.vcf.gz > ${work_path}/imputed_file/high_dep_100.chr${i}_imputed_list.txt
/software/GLIMPSE-1.1.1/ligate/bin/GLIMPSE_ligate --input ${work_path}/imputed_file/high_dep_100.chr${i}_imputed_list.txt --output ${work_path}/imputed_file_merged/high_dep_100.chr${i}_imputed.vcf
/software/bin/bgzip ${work_path}/imputed_file_merged/high_dep_100.chr${i}_imputed.vcf
/software/bcftools-1.15/bcftools index -f ${work_path}/imputed_file_merged/high_dep_100.chr${i}_imputed.vcf.gz
echo \"process end at : \"
date" > ${work_path}/bin/tasks/sbatch_ligate_chr${i}.txt
sbatch ${work_path}/bin/tasks/sbatch_ligate_chr${i}.txt
done
