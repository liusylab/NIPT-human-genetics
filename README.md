<p align="center">
  <a href="https://github.com/liusylab/NIPT-human-genetics">
    <img height="300" src="https://github.com/liusylab/NIPT-human-genetics/blob/main/docs/assets/images/NIPT-Human-Genetics.png">
  </a>
  <h1 align="center">NIPT-human-genetics</h1>
</p>


### Introduction
NIPT-human-genetics is a set of tools for analysing large-scale NIPT sequencing data for human genetic investigation such as SNP detection and allele frequency estimation, individual genotype imputation, kinship estimation, population structure inference and genome-wide association studies.

### Citation
- Liu et al., [Utilizing Non-Invasive Prenatal Test Sequencing Data Resource for Human Genetic Investigation](https://www.biorxiv.org/content/10.1101/2023.12.11.570976v1).BioRxiv (2023)
- Liu et al., [genomic analyses from non-invasive prenatal testing reveal genetic associations, patterns of viral infections, and chinese population history](https://doi.org/10.1016/j.cell.2018.08.016).Cell 175.2 (2018): 347-359

### Pre-requistes
#### Install BWA, Samtools, GATK, BaseVar, GLIMPSE, QUILT, PLINK and EMU
- [BWA](https://github.com/lh3/bwa): <https://github.com/lh3/bwa>
- [GATK](https://github.com/broadinstitute/gatk): <https://github.com/broadinstitute/gatk>
- [Samtools](https://github.com/samtools/samtools/blob/develop/INSTALL): <https://github.com/samtools/samtools/blob/develop/INSTALL>
- [BaseVar](https://github.com/ShujiaHuang/basevar/tree/master): <https://github.com/ShujiaHuang/basevar/tree/master>
- [GLIMPSE](https://odelaneau.github.io/GLIMPSE/docs/installation): <https://odelaneau.github.io/GLIMPSE/docs/installation>
- [QUILT](https://github.com/rwdavies/QUILT/blob/master/README_QUILT1.md): <https://github.com/rwdavies/QUILT/blob/master/README_QUILT1.md> (optional)
- [PLINK](https://www.cog-genomics.org/plink/2.0/): <https://www.cog-genomics.org/plink/2.0/>
- [EMU](https://github.com/Rosemeis/emu): <https://github.com/Rosemeis/emu> (optional)

#### Download reference datasets
- [Human genome reference](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz): <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz>
- [GATK bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle): <https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle>
- 

### System requirements
In order to use the set of tools in NIPT-human-genetics, we require a modern Linux operating system, a version of GCC and JAVA.
The enclosed shell examples were run with SLURM workload manager in a Red Hat 8.4.1-1 system. Commands can be adapted to other workload manager and Linux systems.

### Installation
- wget https://github.com/liusylab/NIPT-human-genetics/archive/refs/heads/main.zip
- unzip main.zip
- cd NIPT-human-genetics-main

### Step-by-step tutorial
### Simulation experiments assessing the performance of BaseVar (optional)
- cd basevar_simulation
- [step1.basevar.simulation.sh](./basevar_simulation/step1.basevar.simulation.sh)
- sh step1.basevar.simulation.sh

### Maximum likelihood model for SNP discovery and allele frequency estimation with BaseVar
- [step1_workflow_alignment.sh](./basevar/step1_workflow_alignment.sh)
- [step2_step2_basevar.sh](./basevar/step2_basevar.sh)

### Gibbs sampling and hidden markov model for genotype imputation
#### genotype imputation using GLIMPSE (version 1.1.1)
- [step1_reference_panle_prepare.sh](./glimpse_imputation/step1_reference_panle_prepare.sh)
- [step2_Computing_GLs.sh](./glimpse_imputation/step2_Computing_GLs.sh)
- [step2_3_merge_GLs_chunk.sh](./glimpse_imputation/step2_3_merge_GLs_chunk.sh)
- [step3_phase.sh](./glimpse_imputation/step3_phase.sh)
- [step4_ligate.sh](./glimpse_imputation/step4_ligate.sh)
- [step5_accuracy.sh](./glimpse_imputation/step5_accuracy.sh)

#### genotype imputation using QUILT (version 1.0.4)
- [quilt_1KGP.sh](./quilt_imputation/quilt_1KGP.sh)

### kinship estimation using PLINK (v2.00a3LM)
- [step1_extract_deep_vcf_sample.sh](./kinship/step1_extract_deep_vcf_sample.sh)
- [step2_plink_2_kinship.sh](./kinship/step2_plink_2_kinship.sh)
- [step3_MERGE.R](./kinship/step3_MERGE.R)

### Principal component analyses using PLINK (v2.00a3LM) or EMU (v.0.9)
- [emu10.sh](./pca/emu10.sh)
- [plink10.sh](./pca/plink10.sh)

### Genome-wide association studies using PLINK (v2.00a3LM)
- [gwas.plink.sh](./gwas/gwas.plink.sh)

