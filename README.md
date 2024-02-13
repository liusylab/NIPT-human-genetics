# NIPT-human-genetics
## A tutorial to utilize NIPT sequencing data resource for human genetic investigation

### Maximum likelihood model for SNP discovery and allele frequency estimation
BaseVar (v0.8.0): <https://github.com/ShujiaHuang/basevar/releases/tag/v0.8.0>

### Gibbs sampling and hidden markov model for genotype imputation
#### genotype imputation using GLIMPSE (version 1.1.1)
- step1_reference_panle_prepare.sh
- step2_Computing_GLs.sh
- step2_3_merge_GLs_chunk.sh
- step3_phase.sh
- step4_ligate.sh
- step5_accuracy.sh

#### genotype imputation using QUILT (version 1.0.4)
- quilt_1KGP.sh

### kinship estimation using PLINK (v2.00a3LM)
- step1_extract_deep_vcf_sample.sh
- step2_plink_2_kinship.sh
- step3_MERGE.R

### Principal component analyses using PLINK (v2.00a3LM) or EMU (v.0.9)
- emu10.sh
- plink10.sh

### Genome-wide association studies using PLINK (v2.00a3LM)
- gwas.plink.sh

