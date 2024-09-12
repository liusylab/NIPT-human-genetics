<p align="center">
  <h1 align="center">NIPT-human-genetics</h1>
</p>
<p align="center">
  <a href="https://github.com/liusylab/NIPT-human-genetics">
    <img height="300" src="https://github.com/liusylab/NIPT-human-genetics/blob/main/docs/assets/images/NIPT-Human-Genetics.png" align="center">
  </a>
</p>


Introduction
------------

**NIPT-human-genetics** is a **semi-automated** workflow designed for the analysis of large-scale ultra-low-pass non-invasive prenatal test (NIPT) sequencing data in human genetic studies. This includes SNP detection, allele frequency estimation, individual genotype imputation, kinship estimation, population structure inference, and genome-wide association studies (GWAS). The workflow has been utilized in several studies, including [Liu et al](<https://doi.org/10.1016/j.cell.2018.08.016>), [Yang et al](<https://doi.org/10.1182/blood.2023021925>), [Zhen et al](<https://doi.org/10.1007/s00125-023-06065-5>), as well as in other works currently under reviews or in press. The term **semi-automated** indicates that while the workflow faciliates the analysis process, users must manually submit jobs according to their specific computational systems settings, allowing for flexibity in execution.

**NIPT-human-genetics** includes the following modules:
- Alignment and statistical analysis of NIPT sequencing data
- SNP discovery and allele frequency estimation using a maximum likelihood model implemented in BASEVAR
- Genotype imputation using a hidden Markov model implemented by GLIMPSE or QUILT
- Kinship estimation
- Principal component analysis (PCA)
- Genome-wide association study (GWAS)

Additional moduldes, such as more accurate genotype imputation, kinship inference and enhanced GWAS analysis, are currently under development. 

Citation
--------

- Liu et al. [Utilizing Non-Invasive Prenatal Test Sequencing Data for Human Genetic Investigation](https://www.biorxiv.org/content/10.1101/2023.12.11.570976v1). BioRxiv (2023)

- Liu, S. et al. [Genomic Analyses from Non-invasive Prenatal Testing Reveal Genetic Associations, Patterns of Viral Infections , and Chinese Population History](https://doi.org/10.1016/j.cell.2018.08.016). Cell 175, 347-359.e14 (2018).



Pre-requistes
-------------

### Install BWA, Samtools, GATK, BaseVar, GLIMPSE, QUILT, PLINK and EMU
- [BWA](https://github.com/lh3/bwa): <https://github.com/lh3/bwa>
- [GATK](https://github.com/broadinstitute/gatk): <https://github.com/broadinstitute/gatk>
- [Samtools](https://github.com/samtools/samtools/blob/develop/INSTALL): <https://github.com/samtools/samtools/blob/develop/INSTALL>
- [BaseVar](https://github.com/ShujiaHuang/basevar): <https://github.com/ShujiaHuang/basevar>
- [GLIMPSE](https://odelaneau.github.io/GLIMPSE/docs/installation): <https://odelaneau.github.io/GLIMPSE/docs/installation>
- [QUILT](https://github.com/rwdavies/QUILT/blob/master/README_QUILT1.md): <https://github.com/rwdavies/QUILT/blob/master/README_QUILT1.md> (optional)
- [PLINK](https://www.cog-genomics.org/plink/2.0/): <https://www.cog-genomics.org/plink/2.0/>
- [EMU](https://github.com/Rosemeis/emu): <https://github.com/Rosemeis/emu> (optional)

### Download human reference datasets

- **Human genome reference**: <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz>
- **GATK bundle**: <https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle>


### System requirements

In order to use the set of tools included in **NIPT-human-genetics**, a modern Linux operating system is required, along with compatible versions of GCC and JAVA.

The accompanying shell [example](https://github.com/liusylab/NIPT-human-genetics/tree/main/example) was executed using the SLURM workload manager on a *Red Hat 8.4.1-1 system*. However, the commands regarding working directory, pathway for the required softwares, and the workload manager should be adapted for use according to the user's system and personal settings.


Installation
------------

```bash
$ git clone --recursive  https://github.com/liusylab/NIPT-human-genetics.git
```

or

```bash
$ git clone --recursive  git@github.com:liusylab/NIPT-human-genetics.git
```

> **WARNING:** Please try several times if fail to clone the data causing by the network problem.


Quick Start
-----------
Users can refer to this example [example](https://github.com/liusylab/NIPT-human-genetics/tree/main/example) for quick start, where we analyze 10 simulated NIPT samples with sequencing depth of around 0.08 to 0.1 fold.
The input data used in this example (2.1G) exceeds the storage capacity available on GitHub, but you can download it from [zenodo](https://zenodo.org/records/13382182) and place it in `example/data` if you wish to try the example byself. Alternatively, you can begin analyzing your own NIPT data using the workflow modules provided in the `example/bin` directory.


Module 1: Alignment and statistics
----------------------------------

**Bash script for this module [module1.alignment.sh](./example/bin/module1.alignment.sh)**

Below is an explanation of the analyses in this module:


**1. read alignment using bwa**

We used the BWA single-end alignment model to map the single-end reads (typically 35 bp) to the latest human genome reference. The alignment reads were then sorted using Samtools, followed by the removal of potential PCR duplicates, and an index was generated for the BAM file.

```bash
# set parameter
hg38=Homo_sapiens_assembly38.fasta
hg38_index_prefix=Homo_sapiens_assembly38.fasta.gz
gatk_bundle_dir=hg38_gatk_bundle/hg38
lane_id=$1
sample_id=$2
fq=$3

$bwa aln -e 10 -t 8 -i 5 -q 0 $hg38_index_prefix $fq > $outdir/${sample_id}.sai && \
$bwa samse -r "@RG\tID:${lane_id}\tPL:COMPLETE\tSM:${sample_id}" $hg38_index_prefix $outdir/${sample_id}.sai $fq | $samtools view -h -Sb - > $outdir/${sample_id}.bam && echo "** bwa done **" && \
$samtools sort -@ 8 -O bam -o $outdir/${sample_id}.sorted.bam $outdir/${sample_id}.bam && echo "** bam sorted done **" && \
$samtools rmdup $outdir/${sample_id}.sorted.bam $outdir/${sample_id}.sorted.rmdup.bam && echo "** rmdup done **" && \
$samtools index $outdir/${sample_id}.sorted.rmdup.bam
```

**2. re-alignment with GATK**

We used GATK to realign indels in the NIPT reads based on known indel information from prior studies.

```bash
$java -Xmx15g -jar $gatk \
    -T RealignerTargetCreator \
    -R $hg38 \
    -I $outdir/${sample_id}.sorted.rmdup.bam \
    -known $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -known $gatk_bundle_dir/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -o $outdir/${sample_id}.indel_target_intervals.list

$samtools index $outdir/${sample_id}.sorted.rmdup.realign.bam
```

**3. BQSR base quality score recalibration with GATK**

Additionally, GATK was employed to recalibrate base quality in the NIPT reads, using known SNPs and indels as references.

```bash
$java -jar $gatk \
    -T BaseRecalibrator \
    -nct 8 \
    -R $hg38 \
    -I $outdir/${sample_id}.sorted.rmdup.realign.bam \
    --knownSites $gatk_bundle_dir/dbsnp_146.hg38.vcf.gz \
    --knownSites $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --knownSites $gatk_bundle_dir/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -o $outdir/${sample_id}.recal_data.table

$java -jar $gatk \
    -T PrintReads \
    -nct 8 \
    -R $hg38 \
    --BQSR $outdir/${sample_id}.recal_data.table \
    -I $outdir/${sample_id}.sorted.rmdup.realign.bam \
    -o $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam
```

**4. bam statistics with samtools and bedtools**

We used Samtools and Bedtools to calculate alignment statistics for the alignment files.

```bash
$samtools stats $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam > $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bamstats
$bedtools genomecov -ibam $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam -bga -split | bgzip > $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz && tabix -p bed $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz
```

Module 2: SNP detection and allele frequency estimation with BaseVar
--------------------------------------------------------------------

You can adhere to the guidelines provided [here](./basevar) to install **BaseVar** before proceeding with the subsequent steps. This ensures that the installation process is completed successfully, enabling you to carry out the following procedures seamlessly.


**Bash script for this module [module2_basevar.sh](./example/bin/module2.basevar.sh)**

Explanation of the analyses in this module:

In Module 2, we begin conducting some analyses in parallel. In the example, we process the data and perform variant detection and allele frequency estimation in 5 million basepair non-overlapping windows. To facilitate this parallelization, we use the pipeline generator [**create_pipeline.py**](./basevar/scripts/create_pipeline.py), which distributes the computational tasks based on the --delta parameter across a specific chromosome defined by the -c parameter.

```bash
$ python create_pipeline.py -R $ref --ref_fai $ref_fai -c chr20 --delta 5000000 -t 20 -L $bamlist -o $outdir > basevar.chr20.sh
```

Within BaseVar, maximum likelihood and likelihood ratio models are employed to determine the polymorphism of a genomic position and estimate the allele frequencies. Detailed matematical documentation can be found [here](https://www.biorxiv.org/content/10.1101/2023.12.11.570976v1). The latest version of BaseVar is available for download via this [link](https://github.com/ShujiaHuang/basevar) and you can install `basevar` locally by following the guideline of this [documents](https://github.com/ShujiaHuang/basevar?tab=readme-ov-file#installation).

To review each of the parameters, you can use `basevar basetype -h` in Linux/MacOS terminal.

```bash
$ /path/to/basevar basetype -h

BaseVar: A software for calling variants efficiently from low-pass whole genome sequencing data.

About: Calling variants by BaseVar.
Usage: basevar basetype [options] <-R Fasta> <--output-vcf> <--output-cvg> [-I input] ...

optional arguments:
  -I, --input=FILE             BAM/CRAM file containing reads.
  -L, --align-file-list=FILE   BAM/CRAM files list, one file per row.
  -R, --reference FILE         Input reference fasta file.

  -m, --min-af=float           Setting prior precision of MAF and skip uneffective caller positions.
                               Usually you can set it to be min(0.001, 100/x), x is the number of input
                               BAM files.[min(0.001, 100/x)]. In generally, you don't have to worry about
                               this parameter.
  -q, --mapq=INT               Only include reads with mapping quality >= INT. [10]
  -B, --batch-count=INT        INT simples per batchfile. [200]
  -t, --thread=INT             Number of threads. [4]

  -G, --pop-group=FILE         Calculating the allele frequency for specific population.
  -r, --regions=chr:start-end  Skip positions which not in these regions. This parameter could be a list
                               of comma deleimited genome regions(e.g.: chr:start-end) or a file contain
                               the list of regions.
  --output-vcf FILE            Output VCF file.
  --output-cvg FILE            Output position coverage file.

  --filename-has-samplename    If the name of bamfile is something like 'SampleID.xxxx.bam', set this
                               argrument could save a lot of time during get the sample id from BAMfile
                               header information.
  --smart-rerun                Rerun process by checking batchfiles.
  -h, --help                   Show this help message and exit.

```


Here is a simple example for running basevar: 

```bash
$basevar basetype -R $hg38 \
    -B 200 \
    -t 20 \
    -L bamfile.list \
    -r chr2,chr11:5246595-5248428,chr17:41197764-41276135 \
    --output-vcf test.vcf.gz \
    --output-cvg test.cvg.tsv.gz
```

>**Plugins**: Simulation experiments for assessing the performance of BaseVar (optional)

The following bash script is available for evaluating the performance of Basevar based on computational simulations: [plugin.basevar.simulation.sh](<./basevar_simulation/plugin.basevar.simulation.sh>)


```bash
$ cd basevar_simulation
$ sh plugin.basevar.simulation.sh
```


Module 3: Genotype imputation
-----------------------------

### Option 1: Genotype imputation using GLIMPSE (version 1.1.1)

**Bash script for this option are divided into five steps:**

- [module3.glimpse.s1.reference_panel_prepare.sh](./example/bin/module3.glimpse.s1.reference_panel_prepare.sh)
- [module3.glimpse.s2.computeGLs.sh](./example/bin/module3.glimpse.s2.computeGLs.sh)
- [module3.glimpse.s3.mergeGLs.sh](./example/bin/module3.glimpse.s3.mergeGLs.sh)
- [module3.glimpse.s4.phase.sh](./example/bin/module3.glimpse.s4.phase.sh)
- [module3.glimpse.s5.ligate.sh](./example/bin/module3.glimpse.s5.ligate.sh)


**Explanation of the analyses in this option:**

**Step1: Preparing the reference panel and imputation chunks**

After downloading the reference panel or computing our own reference panel using some phasing algorithms, the imputation algorithm we use requires normlization and the selection of only bi-allelic SNPs passing certain minor allele frequency threshold from all variants in the reference panel. In addition, prepare chunk files so as to allow glimpse to impute variants by chunks. This process can be accomplished using the following commands.

- Shell scripts: [module3.glimpse.s1.reference_panel_prepare.sh](./example/bin/module3.glimpse.s1.reference_panel_prepare.sh)

```bash
## download 1KGP reference panel
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz{,.tbi} -P $reference_path

## Conduct normalization and filtration of the reference panel
$bcftools norm -m -any ${reference_path}/chr${i}.vcf.gz -Ou --threads 8 | $bcftools view -m 2 -M 2 -v snps -i 'MAF>0.001' --threads 8 -Oz -o ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.vcf.gz
$bcftools index -f ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.vcf.gz

## Extracting variable positions in the reference panel
$bcftools view -G -m 2 -M 2 -v snps ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.vcf.gz -Oz -o ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz --threads 8
$bcftools index -f ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz
$bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz | bgzip -c > ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.tsv.gz
$tabix -s1 -b2 -e2 ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.tsv.gz

##Prepare chunks for genotype imputation
$GLIMPSE_chunk --input ${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz --region chr${i} --window-size 2000000 --buffer-size 200000 --output ${work_path}/chunks.G10K.chr${i}.txt
```

**Step2: Computing Genotype Likelihoods**

In this step, BCFtools is used to compute genotype likelihood. The jobs are processed in parallel based on chromosomes or regions, depending on the available computational resources.

- Shell scripts: [module3.glimpse.s2.computeGLs.sh](./example/bin/module3.glimpse.s2.computeGLs.sh)

```bash
VCF=${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.vcf.gz
TSV=${work_path}/reference_file/chr${i}.biallelic.snp.maf0.001.sites.tsv.gz

$bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T \${VCF} -r chr${i} $line -Ou | bcftools call -Aim -C alleles -T \${TSV} -Oz -o ${work_path}/GL_file/${name}.chr${i}.vcf.gz
$bcftools index -f ${work_path}/GL_file/${name}.chr${i}.vcf.gz
```

**Step3: Merge the genotype likelihood**

Merge the genotype likelihood files and generate chunks to facilitate phasing and genotype imputation.

- Shell scripts: [module3.glimpse.s3.mergeGLs.sh](./example/bin/module3.glimpse.s3.mergeGLs.sh)

```bash
ls ${work_path}/GL_file/*.chr${i}.vcf.gz > ${work_path}/GL_file/high_dep_100.chr${i}_GL_list.txt
$bcftools merge -m none -r chr${i} -Oz -o ${work_path}/GL_file_merged/high_dep_100.chr${i}.vcf.gz -l ${work_path}/GL_file/high_dep_100.chr${i}_GL_list.txt
$bcftools index -f ${work_path}/GL_file_merged/high_dep_100.chr${i}.vcf.gz
```

**Step4: Phasing by GLIMPSE**

Phase the variants using a hidden Markov model based on genotype likelihood information, paralleling the tasks by chunks defined in Step 1. 

- Shell scripts: [module3.glimpse.s4.phase.sh](./example/bin/module3.glimpse.s4.phase.sh)

```bash
$GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
$bgzip ${OUT}
$bcftools index -f ${OUT}.gz
```

**Step5: Ligate**

Merge the imputation results from the chunks.

- Shell scripts: [module3.glimpse.s5.ligate.sh](./example/bin/module3.glimpse.s5.ligate.sh)

```bash
ls ${work_path}/imputed_file/high_dep_100.chr${i}.*.imputed.vcf.gz > ${work_path}/imputed_file/high_dep_100.chr${i}_imputed_list.txt
$GLIMPSE_ligate --input ${work_path}/imputed_file/high_dep_100.chr${i}_imputed_list.txt --output ${work_path}/imputed_file_merged/high_dep_100.chr${i}_imputed.vcf
$bgzip ${work_path}/imputed_file_merged/high_dep_100.chr${i}_imputed.vcf
$bcftools index -f ${work_path}/imputed_file_merged/high_dep_100.chr${i}_imputed.vcf.gz
```

**Plugins**: Calculating the accuracy (Optional)

If high-depth sequencing data is available from some NIPT samples, these high-depth genotypes can be used as the reference set. By comparing the imputed genotype dosage to the high-depth genotypes, we can assess imputation accuracy, specifically by calculating the squared Person correlation coefficient between the true genotypes and the imputed genotype dosages. The following `bcftools` commands can be used to compute imputation accuracy based on allele frequency bins.

```bash
$bcftools stats $true_set ${work_path}/imputed_file_merged/high_dep_100.chr20_imputed.vcf.gz -s - -t chr20 --af-tag "AF" --af-bins "0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.5, 1" > ${work_path}/accuracy/high_dep_100.chr20_imputed.txt
```

### Option 2: Genotype imputation using QUILT (version 1.0.4)

**Bash script for this option are divided into two steps:**

- Step1: [module3.quilt.s1.reference_panel_prepare.sh](./example/bin/module3.quilt.s1.reference_panel_prepare.sh)
- Step2: [module3.quilt.s2.imputation.sh](./example/bin/module3.quilt.s2.imputation.sh)

**Explanation of the analyses in this option**:

**Step1: Preparing the reference panel**

Similar to GLIMPSE, the imputation algorithm of QUILT also involves selecting only bi-allelic SNPs that meet a specified minor allele frequency threshold from all variants in the reference panel. Additionally, QUILT requires the reference panel to be provided in the hap-legend-sample file format.

```bash
#Conduct normalization and filtration of the reference panel
$bcftools norm -m -any ${reference_path}/$prefix.vcf.gz -Ou --threads 8 | $bcftools view -m 2 -M 2 -v snps -i 'MAF>0.001' --threads 8 -Oz -o $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz
$bcftools index -f $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz
zcat $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz | grep '^#'> $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos.vcf && zcat $reference_path/$prefix.biallelic.snp.maf0.001.vcf.gz | awk '!/^#/ {print \$0}'|awk 'length(\$4)==1'|awk 'length(\$5)==1'| awk -F '\\t' '!a[\$2]++' >> $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos.vcf
$bgzip -f $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos.vcf

#Converting formats
$bcftools convert --haplegendsample $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos.vcf.gz
sed -i 's/sample population group sex/SAMPLE POP GROUP SEX/g' $reference_path/$prefix.biallelic.snp.maf0.001.unique_pos.samples
```

**Step2: Conduct genotype imputation**

Conduct the genotype imputation by chromosomes or by defined regions.

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

Module 4: kinship estimation using PLINK (v2.00a3LM)
----------------------------------------------------

**Bash script for module 4 are divided into two steps:**

- Step1: [module4.plink.kinship.s1.sh](./example/bin/module4.plink.kinship.s1.sh)
- Step2: [module4.plink.kinship.s2.sh](./example/bin/module4.plink.kinship.s2.sh)


**Step 1: Compute kinship with plink2**

```bash
$plink2 --vcf $imputed_vcf --make-king-table --out $kinship_outdir/$prefix --threads 24
```

**Step 2: Exclude samples based on kinship**

```bash
$plink2 --vcf $imputed_vcf --king-cutoff-table 0.177 --out $kinship_outdir/$prefix --threads 24
```


Module 5: Principal component analyses using PLINK (v2.00a3LM) or EMU (v.0.9)
-----------------------------------------------------------------------------

**Bash script for module 5 have two options:**

- Step1: [module5.emu.pca.sh](./example/bin/module5.emu.pca.sh)
- Step2: [module5.plink.pca.sh](./example/bin/module5.plink.pca.sh)


**Option1: Perform PCA with EMU**

```bash
$emu -m -p $imputed_vcf -e 10 -t 30 -o $imputed_vcf.emu10
```

**Option2: Perform PCA with PLINK**

```bash
$plink2 --maf 0.05 --vcf $imputed_vcf dosage=DS --pca 10 --out $imputed_vcf.pca10 --threads 30
```

Module 6: Genome-wide association studies by using PLINK (v2.00a3LM)
--------------------------------------------------------------------

**Bash script for module 6 have two options:**

- Step1: [module6.gwas.s1.sh](./example/bin/module6.gwas.s1.sh)
- Step2: [module6.gwas.s2.sh](./example/bin/module6.gwas.s2.sh)

**Step 1: normalize phenotype**

```bash
script1=gwas/get.pheno.py
script2=gwas/generate_plink_tasks.py
allpheno_table=example/data/phenotype.txt  
allsample=./database/pheno.table
vcflist=$imputed_vcf.list
covariates=example/data/covariates.txt

$python  $allpheno_table $allsample $pheno >$outdir/${pheno}_pheno.table
```

**Step2: conduct gwas with plink**

```bash
$plink2 --vcf $imputed_vcf dosage=DS --fam $fam  --covar $covariates --covar-variance-standardize --glm --out .plink --threads 2 --memory 10000 requir
$python $script2 $vcflist $outdir $covar9 $pheno $phenotable $hweinfo linear > $outdir/$pheno/bin/module1.plink.work.sh

```



