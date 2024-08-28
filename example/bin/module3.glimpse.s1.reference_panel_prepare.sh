#!/bin/bash

bcftools=/share/home/lsy_liusiyang/software/bcftools-1.15/bcftools
bgzip=/share/home/lsy_liusiyang/software/tabix-0.2.6/bgzip
tabix=/share/home/lsy_liusiyang/software/tabix-0.2.6/tabix
GLIMPSE_chunk=/share/home/lsy_liusiyang/software/GLIMPSE-1.1.1/chunk/bin/GLIMPSE_chunk

glimpse_outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/glimpse_output
basevar_outdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/basevar_output
reference_path=$glimpse_outdir/reference_file
mkdir -p $reference_path

for i in 20 #can enumerate chromosomes by {1..22} X
do
#download 1KGP reference panel or prepare own reference panels
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz{,.tbi} -P $reference_path

prefix=`basename $reference_path/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz|sed -e 's/\.vcf\.gz$//'`

#prepare tasks
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=07-00:00:00
#SBATCH --mem=20480
#SBATCH --partition=cpu
#SBATCH --job-name=glimpse.s1.chr$i
#SBATCH --output=glimpse.s1.chr$i.slurm.out
#SBATCH --exclude=c007
echo \"process will start at :\"
date" > glimpse.s1.chr$i.tmp.sh

echo "
#Conduct normalization and filtration of the reference panel 
$bcftools norm -m -any ${reference_path}/$prefix.vcf.gz -Ou --threads 8 | $bcftools view -m 2 -M 2 -v snps -i 'MAF>0.001' --threads 8 -Oz -o $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz
$bcftools index -f $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz

#Extracting variable positions in the reference panel
#This process is better replaced with extracting the basevar sites in previous step by replacing $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz with basevar ouput file
$bcftools view -G -m 2 -M 2 -v snps $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz -Oz -o $reference_path/$prefix.biallelic.snp.maf0.001.sites.vcf.gz --threads 8
$bcftools index -f $reference_path/$prefix.biallelic.snp.maf0.001.sites.vcf.gz
$bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' $reference_path/$prefix.biallelic.snp.maf0.001.sites.vcf.gz | $bgzip -c > $reference_path/$prefix.biallelic.snp.maf0.001.sites.tsv.gz
$tabix -s1 -b2 -e2 $reference_path/$prefix.biallelic.snp.maf0.001.sites.tsv.gz

$GLIMPSE_chunk --input $reference_path/$prefix.biallelic.snp.maf0.001.sites.vcf.gz --region chr${i} --window-size 2000000 --buffer-size 200000 --output $reference_path/$prefix.chunks.txt
echo \"process end at : \"
date" >> glimpse.s1.chr$i.tmp.sh

done

