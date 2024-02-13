#!/bin/bash
#!/bin/sh
#SBATCH -J make_quilt.chr20.1.5000000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --time=07-00:00:00
#SBATCH --mem=10240
#SBATCH -p cpu
#SBATCH --output=slurm_quilt.chr20.1.5000000
#SBATCH --exclude=c007,c003,c005,c006
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
echo "We're doing the job in quilt.chr20.1.5000000"
#echo "We'll save it in ${1}"
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

echo "work done"
echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date
