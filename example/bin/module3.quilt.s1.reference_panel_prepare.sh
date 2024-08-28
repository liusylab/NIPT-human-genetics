#!/bin/bash
bgzip=/share/home/lsy_liusiyang/software/tabix-0.2.6/bgzip
bcftools=/share/home/lsy_liusiyang/software/bcftools-1.15/bcftools
quilt_outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/quilt_output
reference_path=$quilt_outdir/reference_file

for i in 20 #can enumerate chromosomes by {1..22} X
do
#download 1KGP reference panel or prepare own reference panels
#wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz{,.tbi} -P $reference_path

prefix=`basename $reference_path/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz|sed -e 's/\.vcf\.gz$//'`

#prepare tasks
echo "#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=07-00:00:00
#SBATCH --mem=20480
#SBATCH --partition=cpu
#SBATCH --job-name=quilt.s1.chr$i
#SBATCH --output=quilt.s1.chr$i.slurm.out
#SBATCH --exclude=c007
echo \"process will start at :\"
date" > quilt.s1.chr$i.tmp.sh

echo "
#Conduct normalization and filtration of the reference panel
$bcftools norm -m -any ${reference_path}/$prefix.vcf.gz -Ou --threads 8 | $bcftools view -m 2 -M 2 -v snps -i 'MAF>0.001' --threads 8 -Oz -o $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz
$bcftools index -f $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz
zcat $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz | grep '^#'> $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos.vcf && zcat $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz | awk '!/^#/ {print \$0}'|awk 'length(\$4)==1'|awk 'length(\$5)==1'| awk -F '\\t' '!a[\$2]++' >> $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos.vcf 
$bgzip -f $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos.vcf

#Converting formats
$bcftools convert --haplegendsample $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos.vcf.gz 
sed -i 's/sample population group sex/SAMPLE POP GROUP SEX/g' $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos.samples

echo \"process end at : \"
date" >> quilt.s1.chr$i.tmp.sh
done
