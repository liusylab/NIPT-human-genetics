ref=/share/home/lsy_liusiyang/20220708_Alignment/hg38/Homo_sapiens_assembly38.fasta
ref_index_prefix=/share/home/lsy_liusiyang/20220708_Alignment/hg38/Homo_sapiens_assembly38.fasta
gatk_bundle_dir=/share/home/lsy_liusiyang/20220708_Alignment/hg38
fqlist=/share/home/lsy_liusiyang/NIPT-human-genetics/example/data/fq.list

workflow=/share/home/lsy_liusiyang/NIPT-human-genetics/bwa_alignment/alingment_workflow.sh
bwa=/share/home/lsy_liusiyang/software/bwa-0.7.17/bwa
samtools=/share/home/lsy_liusiyang/software/samtools-1.15.1/samtools
gatk=/share/home/lsy_liusiyang/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
java=/usr/bin/java
bedtools=/share/home/lsy_liusiyang/software/bedtools2/bin/bedtools
bgzip=/share/home/lsy_liusiyang/software/tabix-0.2.6/bgzip
tabix=/share/home/lsy_liusiyang/software/tabix-0.2.6/tabix


### Output Because the process is conducted in batches and in parallel, a temporary output directory has been set up
tmpoutdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/1tmp_files
finaloutdir=/share/home/lsy_liusiyang/NIPT-human-genetics/example/outdir/batch1_final_files

cat $fqlist|while read -r line 
do
	fq=`echo $line|awk '{print $1}'`
	cid=`echo $line|awk '{print $2}'`
	lid=`echo $line|awk '{print $3}'`
	samid=`echo $line|awk '{print $4}'`
echo "#!/bin/sh
#SBATCH -J NIPT-human-genetics_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=00-05:05:30
#SBATCH --mem=18432
#SBATCH --partition=cpu
$workflow -f $fq  -c $cid -l $lid -s $samid -a $bwa -t $samtools -g $gatk -j $java -e $bedtools -z $bgzip -x $tabix -r $ref -i $ref_index_prefix -b $gatk_bundle_dir -o $tmpoutdir -u $finaloutdir" >$samid.tmp.sh
done
