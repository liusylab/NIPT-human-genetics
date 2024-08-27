<p align="center">
  <<F2>a href="https://github.com/liusylab/NIPT-human-genetics">
    <img height="300" src="https://github.com/liusylab/NIPT-human-genetics/blob/main/docs/assets/images/NIPT-Human-Genetics.png">
  </a>
  <h1 align="center">NIPT-human-genetics</h1>
</p>


### Introduction
------------

NIPT-human-genetics is a workflow for analysing large-scale NIPT sequencing data for human genetic investigation such as SNP detection, allele frequency estimation, individual genotype imputation, kinship estimation, population structure inference and genome-wide association studies.


### Citation
--------

- Liu et al., [Utilizing Non-Invasive Prenatal Test Sequencing Data for Human Genetic Investigation](https://www.biorxiv.org/content/10.1101/2023.12.11.570976v1). BioRxiv (2023)

- Liu et al., [Genomic analyses from non-invasive prenatal testing reveal genetic associations, patterns of viral infections, and chinese population history](https://doi.org/10.1016/j.cell.2018.08.016). Cell 175.2 (2018): 347-359


### Pre-requistes
-------------

### Install BWA, Samtools, GATK, BaseVar, GLIMPSE, QUILT, PLINK and EMU
- [BWA](https://github.com/lh3/bwa): <https://github.com/lh3/bwa>
- [GATK](https://github.com/broadinstitute/gatk): <https://github.com/broadinstitute/gatk>
- [Samtools](https://github.com/samtools/samtools/blob/develop/INSTALL): <https://github.com/samtools/samtools/blob/develop/INSTALL>
- [BaseVar](https://github.com/ShujiaHuang/basevar/tree/master): <https://github.com/ShujiaHuang/basevar/tree/master>
- [GLIMPSE](https://odelaneau.github.io/GLIMPSE/docs/installation): <https://odelaneau.github.io/GLIMPSE/docs/installation>
- [QUILT](https://github.com/rwdavies/QUILT/blob/master/README_QUILT1.md): <https://github.com/rwdavies/QUILT/blob/master/README_QUILT1.md> (optional)
- [PLINK](https://www.cog-genomics.org/plink/2.0/): <https://www.cog-genomics.org/plink/2.0/>
- [EMU](https://github.com/Rosemeis/emu): <https://github.com/Rosemeis/emu> (optional)

### Download reference datasets

- [Human genome reference](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz): <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz>

- [GATK bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle): <https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle> 


### System requirements

In order to use the set of tools in NIPT-human-genetics, we require a modern Linux operating system, a version of GCC and JAVA.
The enclosed shell [example](https://github.com/liusylab/NIPT-human-genetics/tree/main/example) were run with SLURM workload manager in a Red Hat 8.4.1-1 system. Commands can be adapted to other workload manager and Linux systems.

### Installation
-----------
```bash
$ git clone https://github.com/liusylab/NIPT-human-genetics.git
$ cd NIPT-human-genetics
```

### Quick Start
-----------
Command line users can refer to this example [example](https://github.com/liusylab/NIPT-human-genetics/tree/main/example) for quick start


---------------------------------------------------------------------------------------

### Simulation experiments assessing the performance of BaseVar (optional)

```bash
$ cd basevar_simulation
$ [step1.basevar.simulation.sh](./basevar_simulation/step1.basevar.simulation.sh)
$ sh step1.basevar.simulation.sh
```

---------------------------------------------------------------------------------------

### Step 1: Alignment and statistics

- [step1.alignment.sh](./example/bin/step1.alignment.sh)


**1. shell script for bwa alignment**

```bash
# set parameter

hg38=Homo_sapiens_assembly38.fasta
hg38_index_prefix=Homo_sapiens_assembly38.fasta.gz
gatk_bundle_dir=hg38_gatk_bundle/hg38
lane_id=$1
sample_id=$2
fq=$3

bwa aln -e 10 -t 8 -i 5 -q 0 $hg38_index_prefix $fq > $outdir/${sample_id}.sai && \
bwa samse -r "@RG\tID:${lane_id}\tPL:COMPLETE\tSM:${sample_id}" $hg38_index_prefix $outdir/${sample_id}.sai $fq | $samtools view -h -Sb - > $outdir/${sample_id}.bam && echo "** bwa done **" && \
samtools sort -@ 8 -O bam -o $outdir/${sample_id}.sorted.bam $outdir/${sample_id}.bam && echo "** bam sorted done **" && \
samtools rmdup $outdir/${sample_id}.sorted.bam $outdir/${sample_id}.sorted.rmdup.bam && echo "** rmdup done **" && \
samtools index $outdir/${sample_id}.sorted.rmdup.bam
```

**2. re-alignment with GATK**

```bash
java -Xmx15g -jar $gatk \
    -T RealignerTargetCreator \
    -R $hg38 \
    -I $outdir/${sample_id}.sorted.rmdup.bam \
    -known $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -known $gatk_bundle_dir/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -o $outdir/${sample_id}.indel_target_intervals.list

samtools index $outdir/${sample_id}.sorted.rmdup.realign.bam
```

**3. BQSR base quality score recalibration with GATK**

```bash
java -jar $gatk \
    -T BaseRecalibrator \
    -nct 8 \
    -R $hg38 \
    -I $outdir/${sample_id}.sorted.rmdup.realign.bam \
    --knownSites $gatk_bundle_dir/dbsnp_146.hg38.vcf.gz \
    --knownSites $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --knownSites $gatk_bundle_dir/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -o $outdir/${sample_id}.recal_data.table

java -jar $gatk \
    -T PrintReads \
    -nct 8 \
    -R $hg38 \
    --BQSR $outdir/${sample_id}.recal_data.table \
    -I $outdir/${sample_id}.sorted.rmdup.realign.bam \
    -o $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam
```

**4. bam statistics with samtools and bedtools**

```bash
samtools stats $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam > $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bamstats
bedtools genomecov -ibam $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam -bga -split | bgzip > $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz && tabix -p bed $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz

```

---------------------------------------------------------------------------------------
### Step 2: Maximum likelihood model for SNP discovery and allele frequency estimation with BaseVar
- [step2_basevar.sh](./example/step2.basevar.sh)

```bash
basevar basetype -R $hg38 \
    --regions chr11:5246595-5248428,chr17:41197764-41276135 \
    --batch-count 50 \
    -L bamfile.list \
    --output-vcf test.vcf.gz \
    --output-cvg test.cvg.tsv.gz \
    --nCPU 4
```


---------------------------------------------------------------
### Step 3: Gibbs sampling and hidden markov model for genotype imputation

### Genotype imputation using GLIMPSE (version 1.1.1)

**S1: Preparing the reference panel and chunks**

```bash
## Conduct normalization and filtration of the reference panel
bcftools norm -m -any ${reference_path}/chr${i}.vcf.gz -Ou --threads 8 | $bcftools view -m 2 -M 2 -v snps -i 'MAF>0.001' --threads 8 -Oz -o ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.vcf.gz
bcftools index -f ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.vcf.gz

## Extracting variable positions in the reference panel
bcftools view -G -m 2 -M 2 -v snps ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.vcf.gz -Oz -o ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz --threads 8
bcftools index -f ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz
bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz | bgzip -c > ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.tsv.gz
tabix -s1 -b2 -e2 ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.tsv.gz
```

**S2: Computing Genotype Likelihoods**

```bash
VCF=${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz
TSV=${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.tsv.gz

bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T \${VCF} -r chr${i} $line -Ou | bcftools call -Aim -C alleles -T \${TSV} -Oz -o ${work_path}/GL_file/${name}.chr${i}.vcf.gz
bcftools index -f ${work_path}/GL_file/${name}.chr${i}.vcf.gz
```

**S3: Merge the genotype likelihood by chunk**

```bash
ls ${work_path}/GL_file/*.chr${i}.vcf.gz > ${work_path}/GL_file/high_dep_100.chr${i}_GL_list.txt
bcftools merge -m none -r chr${i} -Oz -o ${work_path}/GL_file_merged/high_dep_100.chr${i}.vcf.gz -l ${work_path}/GL_file/high_dep_100.chr${i}_GL_list.txt
bcftools index -f ${work_path}/GL_file_merged/high_dep_100.chr${i}.vcf.gz
GLIMPSE_chunk --input ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz --region chr${i} --window-size 2000000 --buffer-size 200000 --output ${work_path}/chunks.G10K.chr${i}.txt
```

**S4: Phasing by GLIMPSE**

```bash
GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
bgzip ${OUT}
bcftools index -f ${OUT}.gz
```

**S5: Ligate**

```bash
ls ${work_path}/imputed_file/high_dep_100.chr${i}.*.imputed.vcf.gz > ${work_path}/imputed_file/high_dep_100.chr${i}_imputed_list.txt
GLIMPSE_ligate --input ${work_path}/imputed_file/high_dep_100.chr${i}_imputed_list.txt --output ${work_path}/imputed_file_merged/high_dep_100.chr${i}_imputed.vcf
bgzip ${work_path}/imputed_file_merged/high_dep_100.chr${i}_imputed.vcf
bcftools index -f ${work_path}/imputed_file_merged/high_dep_100.chr${i}_imputed.vcf.gz
```

**S6: Calculating the accuracy (Optinal)**

```bash
bcftools stats $true_set ${work_path}/imputed_file_merged/high_dep_100.chr20_imputed.vcf.gz -s - -t chr20 --af-tag "AF" --af-bins "0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.5, 1" > ${work_path}/accuracy/high_dep_100.chr20_imputed.txt
```

---------------------------------------------------------------
### Step 3: Gibbs sampling and hidden markov model for genotype imputation
### Genotype imputation using QUILT (version 1.0.4)

**Shell scripts for QUILT**


```bash
/software/QUILT/QUILT.R \
--outputdir=quilt_100_1KGP/output/ \
--reference=quilt_100_1KGP/input/hg38/Homo_sapiens_assembly38.fasta \
--buffer=250000 \
--chr=chr20 \
--regionStart=1 \
--regionEnd=5000000 \
--nGibbsSamples=7 \
--n_seek_its=3 \
--bamlist=quilt_100_1KGP/input/baoan_NIPT_110_bamlist.txt \
--reference_haplotype_file=quilt_100_1KGP/input/1KGP_legend_haplotype/uniq_1KGP/legend_haplotype/1KGP.chr20.hap.gz \
--reference_legend_file=quilt_100_1KGP/input/1KGP_legend_haplotype/uniq_1KGP/legend_haplotype/1KGP.chr20.legend.gz \
--genetic_map_file=quilt_100_1KGP/input/genetic_map_file/CHS_chr20.txt.gz \
--reference_populations=CHB,CHS,CDX \
--nGen=1240 \
--save_prepared_reference=TRUE
```


------------------------------------------
### Step 4: kinship estimation using PLINK (v2.00a3LM)

- [step1_extract_deep_vcf_sample.sh](./kinship/step1_extract_deep_vcf_sample.sh)
- [step2_plink_2_kinship.sh](./kinship/step2_plink_2_kinship.sh)
- [step3_MERGE.R](./kinship/step3_MERGE.R)

**Step 1: Extract**

```bash
#!/bin/sh
#SBATCH -J baoan_extract_vcf_sample_10K
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=07-00:00:00
#SBATCH --mem=50G
#SBATCH -p bigmem
#SBATCH --output=slurm_baoan_extract_vcf_sample_10K

workdir=./kinship
bcftools view --force-samples -S identical_samples_kinship_SampleID.txt $workdir/unimpute_vcf_file/baoan_70608_GLs.merged.chr1-22.0.001.vcf.gz -Oz -o $workdir/sample_extract_30n/extract.70608_GLs.merged.chr1-22.0.001.vcf.gz
tabix $workdir/sample_extract_30n/extract.70608_GLs.merged.chr1-22.0.001.vcf.gz

echo "OK"
echo "process end at : "
date

```

**Step 2: Compute kinship with plink2**

```bash
#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=24
#SBATCH --time=06-00:00:00
#SBATCH --mem=50G
#SBATCH --partition=bigmem
#SBATCH --job-name=baoan_kinship
#SBATCH --output=/kinship
workdir=/kinship
plink2 --vcf $workdir/sample_extract_30n/extract.70608_GLs.merged.chr1-22.0.001.vcf.gz --make-king-table --out $workdir/plink2_kinship/70608_GLs.merged.chr1-22.0.001_maf --threads 24
```

**Step 3: merge by R**

```r
#! /usr/local/bin Rscriptlibrary(dplyr)
library(dplyr)
library(data.table)
library(R.utils)
library(tidyr)
workdir<-/kinship
input_file<-fread('$workdir/plink2_kinship/70608_GLs.merged.chr1-22.0.001_maf.kin0',header=T,sep="\t",stringsAsFactors = F)
Sample_list<-fread('kinship_files/kinship_match.txt',header=T,sep="\t",stringsAsFactors = F)
colnames(input_file)<-c("nipt_id.y","nipt_id.x","NSNP","HETHET","IBS0","KINSHIP")
colnames(Sample_list)<-c("nipt_id.x","nipt_id.y","seq_dep.x","seq_dep.y","coef")
merge_file<-merge(input_file,Sample_list,by=c("nipt_id.x","nipt_id.y"))
write.table(merge_file,'$workdir/plink2_kinship_merge/merge.70608_GLs.merged.chr1-22.0.001_maf.kin0',quote=FALSE,row.names=FALSE,col.names=T,sep = "\t")

```


Principal component analyses using PLINK (v2.00a3LM) or EMU (v.0.9)
-------------------------------------------------------------------


**Perform PCA with EMU**


```bash
#!/bin/bash
#SBATCH -J emu10_baoan_all
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=30
#SBATCH --time=50-00:00:00
#SBATCH --mem=204800
#SBATCH -p bigmem
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
emu -m -p GLs.merged.chr1-22.redup.0.05 -e 10 -t 30 -o GLs.merged.chr1-22.redup.0.05.emu10
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
```


**Perform PCA with PLINK**


```bash
#!/bin/sh
#SBATCH -J Baoan_PCA10
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=30
#SBATCH --time=10-05:05:30
#SBATCH --mem=100G
#SBATCH --partition=bigmem
#SBATCH --output=slurm_Baoan_unimputed_PCA10
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
plink2 --maf 0.05 --vcf GLs.merged.chr1-22.0.001.vcf.gz dosage=DS --pca 10 --out NIPT_rm_dup_merged.chr1-22.0.05_pca10 --threads 30

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 10s" 
sleep 10
echo "process end at : "
date
```

Genome-wide association studies by using PLINK(v2.00a3LM)
---------------------------------------------------------

```bash
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
```


