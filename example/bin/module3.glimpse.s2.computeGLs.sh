#!/bin/bash
ref=/share/home/lsy_liusiyang/20220708_Alignment/hg38/Homo_sapiens_assembly38.fasta
bamlist=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/batch1_final_files/bam.list

bcftools=/share/home/lsy_liusiyang/software/bcftools-1.15/bcftools
bgzip=/share/home/lsy_liusiyang/software/tabix-0.2.6/bgzip
tabix=/share/home/lsy_liusiyang/software/tabix-0.2.6/tabix

glimpse_outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/glimpse_output
glpath=$glimpse_outdir/GL_file

mkdir -p $glpath

less -S $bamlist|while read line
do
filename=${line##*/}
name=${filename%%.*}
for i in 20 #to enumerate chromosomes, use for i in {1..22} X
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
#SBATCH --output=GL_$name.chr$i.slurm.out
echo \"process will start at :\"
date" >glimpse.s2.gl.$name.chr$i.tmp.sh

echo "
VCF=${glimpse_outdir}/reference_file/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.biallelic.snp.maf0.001.sites.vcf.gz
TSV=${glimpse_outdir}/reference_file/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.biallelic.snp.maf0.001.sites.tsv.gz
$bcftools mpileup -f ${ref} -I -E -a 'FORMAT/DP' -T \${VCF} -r chr${i} $line -Ou | $bcftools call -Aim -C alleles -T \${TSV} -Oz -o ${glimpse_outdir}/GL_file/${name}.chr${i}.vcf.gz
$bcftools index -f ${glimpse_outdir}/GL_file/${name}.chr${i}.vcf.gz
echo \"process end at : \"
date" >> glimpse.s2.gl.$name.chr$i.tmp.sh

done
done 

