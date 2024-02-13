#!/bin/bash
#!/bin/sh
#SBATCH -J 1KGP_accuracy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --time=07-00:00:00
#SBATCH --mem=10240
#SBATCH -p cpu
#SBATCH --output=slurm_G10K_accuracy
#SBATCH --exclude=c007,c003,c005,c006
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
#unfiltered
work_path=glimpse_pipline/ref_1KGP
true_set=HC.VQSR.PASS.chr20.with.ChinaMapAF.vcf.gz
mkdir -p ${work_path}/accuracy/filtered
bcftools stats $true_set ${work_path}/imputed_file_merged/high_dep_100.chr20_imputed.vcf.gz -s - -t chr20 --af-tag "AF" --af-bins "0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.5, 1" > ${work_path}/accuracy/high_dep_100.chr20_imputed.txt
#
bcftools stats $true_set ${work_path}/imputed_file_merged/high_dep_100.chr20_imputed.subset.0.1.vcf.gz -s - -t chr20 --af-tag "AF" --af-bins "0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.5, 1" > ${work_path}/accuracy/high_dep_100.chr20_imputed.subset.0.1.txt
#
bcftools stats $true_set ${work_path}/imputed_file_merged/high_dep_100.chr20_imputed.subset.0.2.vcf.gz -s - -t chr20 --af-tag "AF" --af-bins "0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.5, 1" > ${work_path}/accuracy/high_dep_100.chr20_imputed.subset.0.2.txt

#filtered_INFO:0.4
bcftools stats $true_set ${work_path}/imputed_file_merged/filtered/filtered_high_dep_100.chr20_imputed.vcf.gz -s - -t chr20 --af-tag "AF" --af-bins "0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.5, 1" >${work_path}/accuracy/filtered/filtered_high_dep_100.chr20_imputed.txt
#
bcftools stats $true_set ${work_path}/imputed_file_merged/filtered/filtered_high_dep_100.chr20_imputed.subset.0.1.vcf.gz -s - -t chr20 --af-tag "AF" --af-bins "0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.5, 1" > ${work_path}/accuracy/filtered/filtered_high_dep_100.chr20_imputed.subset.0.1.txt
#
bcftools stats $true_set ${work_path}/imputed_file_merged/filtered/filtered_high_dep_100.chr20_imputed.subset.0.2.vcf.gz -s - -t chr20 --af-tag "AF" --af-bins "0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.5, 1" > ${work_path}/accuracy/filtered/filtered_high_dep_100.chr20_imputed.subset.0.2.txt

echo "work done"
echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date
