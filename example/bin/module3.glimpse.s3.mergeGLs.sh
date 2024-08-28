#!/bin/sh
bcftools=/share/home/lsy_liusiyang/software/bcftools-1.15/bcftools

glimpse_outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/glimpse_output
glpath=$glimpse_outdir/GL_file
glmergepath=$glimpse_outdir/GL_file_merged

mkdir -p $glmergepath

#prepare tasks
for i in 20 #to enumerate chromosomes, use for i in {1..22} X
do
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=07-00:00:00
#SBATCH --mem=20480
#SBATCH --partition=cpu
#SBATCH --job-name=GL_merge
#SBATCH --output=GL_merge.slurm.out
echo \"process will start at :\"
date" >glimpse.s3.GL_merge.chr${i}.tmp.sh

echo "
ls ${glimpse_outdir}/GL_file/*.chr${i}.vcf.gz > ${glimpse_outdir}/GL_file/glimpse.chr${i}_GL_list.txt
$bcftools merge -m none -r chr${i} -Oz -o ${glimpse_outdir}/GL_file_merged/glimpse.chr${i}.vcf.gz -l ${glimpse_outdir}/GL_file/glimpse.chr${i}_GL_list.txt
$bcftools index -f ${glimpse_outdir}/GL_file_merged/glimpse.chr${i}.vcf.gz

echo \"process end at : \"
date" >> glimpse.s3.GL_merge.chr${i}.tmp.sh

done
