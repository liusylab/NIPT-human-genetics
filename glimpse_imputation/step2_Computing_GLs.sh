#!/bin/sh
work_path=glimpse_pipline/ref_1KGP #1KGP/BIGCS/STROMICS
REFGEN=Homo_sapiens_assembly38.fasta
mkdir -p ${work_path}/bin/tasks
mkdir -p ${work_path}/GL_file
while read line
do
filename=${line##*/}
name=${filename%%.*}
for i in {1..22} X
do
#prepare tasks
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=07-00:00:00
#SBATCH --mem=20480
#SBATCH --partition=cpu
#SBATCH --job-name=GL_$name
#SBATCH --output=${work_path}/bin/tasks/slurm_GL_$name.sh
echo \"process will start at :\"
date
export PATH=/share/home/lsy_guyuqin/miniconda3/envs/env_gyq/bin:\$PATH
VCF=${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz
TSV=${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.tsv.gz
bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T \${VCF} -r chr${i} $line -Ou | bcftools call -Aim -C alleles -T \${TSV} -Oz -o ${work_path}/GL_file/${name}.chr${i}.vcf.gz
bcftools index -f ${work_path}/GL_file/${name}.chr${i}.vcf.gz
echo \"process end at : \"
date" > ${work_path}/bin/tasks/sbatch_GL_$name\_chr${i}.sh
sbatch ${work_path}/bin/tasks/sbatch_GL_$name\_chr${i}.sh
done
done < NIPT_highdepth100_bamlist.txt

