#!/bin/sh
#---
#Perform genome wide association analysis using plink
#Step1 perform inverse rank transformation (quantile transformation) for phenotypic data
#Step2 generate plink shell
#---

outdir=./230914_GWAS/output/
mkdir -p $outdir

script1=./bin/get.pheno.py
script2=./bin/g.worksh_1.py
R=`which Rscript`

allpheno_table=./230914_GWAS/input/coding1.txt  

allsample=./database/pheno.table
vcflist=./230914_GWAS/input/hospital_vcf_list.txt
covar9=./230914_GWAS/input/plink.PC1-5.maternalAge.bmi.gw.xCov.txt
hweinfo=./database/Baoan.all.snp.info.gz
imputevcf=./database/hg38.all.stitch.bed.sorted.dbsnp.gwas.clinvar.snpeff.vcf.gz

for pheno in VA VK
do
    #pheno=$line	
    mkdir -p $outdir/$pheno/input/
    mkdir -p $outdir/$pheno/output/plink/
    mkdir -p $outdir/$pheno/output/plink_merge/
    mkdir -p $outdir/$pheno/bin

#---
#Step1 clean and normalize phenotypes
#---
python $script1 $allpheno_table $allsample $pheno >$outdir/$pheno/input/${pheno}_pheno.table
phenotable=$outdir/$pheno/input/${pheno}_pheno.table

#---
#Step 2 generate plink shell
#---
python $script2 $vcflist $outdir/$pheno/output/plink/ $covar9 $pheno $phenotable $hweinfo linear > $outdir/$pheno/bin/step1.plink.work.sh

done
