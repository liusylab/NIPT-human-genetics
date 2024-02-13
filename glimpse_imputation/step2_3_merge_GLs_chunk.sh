#!/bin/sh
work_path=glimpse_pipline/ref_1KGP
mkdir -p ${work_path}/bin/tasks
mkdir -p ${work_path}/GL_file_merged
#prepare tasks
for i in {1..22} X
do
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=07-00:00:00
#SBATCH --mem=20480
#SBATCH --partition=cpu
#SBATCH --job-name=GL_merge_chunk
#SBATCH --output=${work_path}/bin/tasks/slurm_GL_merge_chunk.sh
echo \"process will start at :\"
date
ls ${work_path}/GL_file/*.chr${i}.vcf.gz > ${work_path}/GL_file/high_dep_100.chr${i}_GL_list.txt
bcftools merge -m none -r chr${i} -Oz -o ${work_path}/GL_file_merged/high_dep_100.chr${i}.vcf.gz -l ${work_path}/GL_file/high_dep_100.chr${i}_GL_list.txt
bcftools index -f ${work_path}/GL_file_merged/high_dep_100.chr${i}.vcf.gz
/software/GLIMPSE-1.1.1/chunk/bin/GLIMPSE_chunk --input ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz --region chr${i} --window-size 2000000 --buffer-size 200000 --output ${work_path}/chunks.G10K.chr${i}.txt
echo \"process end at : \"
date" > ${work_path}/bin/tasks/sbatch_GL_merge_chunk_chr${i}.sh
sbatch ${work_path}/bin/tasks/sbatch_GL_merge_chunk_chr${i}.sh
done
