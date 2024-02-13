#!/bin/sh
reference_path=/mnt/OS5300/public_data/reference_panel
work_path=glimpse_pipline/ref_1KGP #1KGP/BIGCS/STROMICS
mkdir -p ${work_path}/bin/tasks
mkdir -p ${work_path}/reference_file
for i in {1..22} X
do
#prepare tasks
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=07-00:00:00
#SBATCH --mem=20480
#SBATCH --partition=cpu
#SBATCH --job-name=chr$i
#SBATCH --output=${work_path}/bin/tasks/slurm_chr$i.sh
#SBATCH --exclude=c007
echo \"process will start at :\"
date
#Conduct normalization and filtration of the reference panel 
bcftools norm -m -any ${reference_path}/chr${i}.vcf.gz -Ou --threads 8 | bcftools view -m 2 -M 2 -v snps -i 'MAF>0.001' --threads 8 -Oz -o ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.vcf.gz
bcftools index -f ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.vcf.gz
#Extracting variable positions in the reference panel
bcftools view -G -m 2 -M 2 -v snps ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.vcf.gz -Oz -o ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz --threads 8
bcftools index -f ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz
bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz | bgzip -c > ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.tsv.gz
tabix -s1 -b2 -e2 ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.tsv.gz
echo \"process end at : \"
date" > ${work_path}/bin/tasks/sbatch_chr$i.sh
sbatch ${work_path}/bin/tasks/sbatch_chr$i.sh
#sleep 3s
done

