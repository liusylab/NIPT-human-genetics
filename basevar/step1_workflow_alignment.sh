#!/bin/sh
gatk=GenomeAnalysisTK.jar
java=java

### Reference
hg38=Homo_sapiens_assembly38.fasta
hg38_index_prefix=Homo_sapiens_assembly38.fasta.gz

### GATK bundle
gatk_bundle_dir=hg38_gatk_bundle/hg38

### tools
bwa=bwa
samtools=samtools

### Input
lane_id=$1
sample_id=$2
fq=$3

### Output Because the process is conducted in batches and in parallel, a temporary output directory has been set up.
final_outdir=$workdir/nipt_bwa/1final_file/$lane_id
outdir=$workdir/nipt_bwa/1tmp_file/$lane_id


if [ ! -d $final_outdir ]
then mkdir -p $final_outdir
fi

if [ ! -d $outdir ]
then mkdir -p $outdir
fi

######################################################################################
################################### Pipeline ########################################
######################################################################################
echo "We're doing the job in $lane_id"
echo "We are calculating $fq"
echo "We'll save it in $final_outdir"
### step 1: Adjust the number of threads for bwa according to your cluster. The -e 10 parameter makes indel detection more sensitive, though its impact is minimal.
echo ""
time $bwa aln -e 10 -t 8 -i 5 -q 0 $hg38_index_prefix $fq > $outdir/${sample_id}.sai && \
    time $bwa samse -r "@RG\tID:${lane_id}\tPL:COMPLETE\tSM:${sample_id}" $hg38_index_prefix $outdir/${sample_id}.sai $fq | $samtools view -h -Sb - > $outdir/${sample_id}.bam && echo "** bwa done **" && \
    time $samtools sort -@ 8 -O bam -o $outdir/${sample_id}.sorted.bam $outdir/${sample_id}.bam && echo "** bam sorted done **" && \
    time $samtools rmdup $outdir/${sample_id}.sorted.bam $outdir/${sample_id}.sorted.rmdup.bam && echo "** rmdup done **" && \
    time $samtools index $outdir/${sample_id}.sorted.rmdup.bam && echo "** index done **" && touch ${outdir}/bwa_sort_rmdup.finish

if [ ! -f ${outdir}/bwa_sort_rmdup.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bwa_sort_rmdup not done **" && exit
fi

### step 2: realigne 需根据数据的大小改内存的需要
### Realign
echo ""
time $java -Xmx15g -jar $gatk \
    -T RealignerTargetCreator \
    -R $hg38 \
    -I $outdir/${sample_id}.sorted.rmdup.bam \
    -known $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -known $gatk_bundle_dir/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -o $outdir/${sample_id}.indel_target_intervals.list && echo "** RealignerTargetCreator done " && touch ${outdir}/RealignerTargetCreator.finish

if [ ! -f ${outdir}/RealignerTargetCreator.finish ]
then echo "** [WORKFLOW_ERROR_INFO] RealignerTargetCreator not done **" && exit
fi

time $java -Xmx15g -jar $gatk \
    -T IndelRealigner \
    -R $hg38 \
    -I $outdir/${sample_id}.sorted.rmdup.bam \
    -known $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -known $gatk_bundle_dir/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --targetIntervals $outdir/${sample_id}.indel_target_intervals.list \
    -o $outdir/${sample_id}.sorted.rmdup.realign.bam && echo "** IndelRealigner done **" && touch ${outdir}/IndelRealigner.finish

if [ ! -f ${outdir}/IndelRealigner.finish ]
then echo "** [WORKFLOW_ERROR_INFO] IndelRealigner not done **" && exit
fi

### step3. BQSR base quality score recalibration
time $samtools index $outdir/${sample_id}.sorted.rmdup.realign.bam && echo "** index done **" && touch ${outdir}/index_sorted_rmdup_realign_bam.finish

if [ ! -f ${outdir}/index_sorted_rmdup_realign_bam.finish ]
then echo "** [WORKFLOW_ERROR_INFO] index_sorted_rmdup_realign_bam not done **" && exit
fi

echo ""
time $java -jar $gatk \
    -T BaseRecalibrator \
    -nct 8 \
    -R $hg38 \
    -I $outdir/${sample_id}.sorted.rmdup.realign.bam \
    --knownSites $gatk_bundle_dir/dbsnp_146.hg38.vcf.gz \
    --knownSites $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --knownSites $gatk_bundle_dir/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -o $outdir/${sample_id}.recal_data.table && echo "** BaseRecalibrator done " && touch ${outdir}/baseRecalibrator.finish

if [ ! -f ${outdir}/baseRecalibrator.finish ]
then echo "** [WORKFLOW_ERROR_INFO] baseRecalibrator not done **" && exit
fi

time $java -jar $gatk \
    -T PrintReads \
    -nct 8 \
    -R $hg38 \
    --BQSR $outdir/${sample_id}.recal_data.table \
    -I $outdir/${sample_id}.sorted.rmdup.realign.bam \
    -o $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam && echo "** PrintReads done **" && touch ${outdir}/PrintReads.finish

if [ ! -f ${outdir}/PrintReads.finish ]
then echo "** [WORKFLOW_ERROR_INFO] PrintReads not done **" && exit
fi


### Convert Bam to Cram to save storage (因为变异检测似乎对cram文件的读取特别慢，未查出来原因，所以更改流程做完变异检测后再进行bam到cram的转换)
#echo ""
#time samtools view -C -T $hg38 $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam > $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cram && echo "** cram done **" && touch ${outdir}/cram.finish
#if [ ! -f ${outdir}/cram.finish ]
#then echo "** [WORKFLOW_ERROR_INFO] cram not done **" && exit
#fi

time $samtools index $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam && echo "** bam index done **" && touch ${outdir}/bam_index.finish

if [ ! -f ${outdir}/bam_index.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bam index not done **" && exit
fi

#time samtools index $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cram && echo "** cram index done **" && touch ${outdir}/cram_index.finish

#if [ ! -f ${outdir}/cram_index.finish ]
#then echo "** [WORKFLOW_ERROR_INFO] cram index not done **" && exit
#fi

### step 5. bam stats
time $samtools stats $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam > $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bamstats && echo "** bamstats done **" && touch ${outdir}/bamstats.finish

if [ ! -f ${outdir}/bamstats.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bamstats not done **" && exit
fi

### step 6. bedtools
time bedtools genomecov -ibam $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam -bga -split | bgzip > $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz && tabix -p bed $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz && echo "** sorted.rmdup.realign.BQSR.cvg.bed.gz done **" && touch ${outdir}/sorted_rmdup_realign_BQSR_cvg_bed.finish

if [ ! -f ${outdir}/sorted_rmdup_realign_BQSR_cvg_bed.finish ]
then echo "** [WORKFLOW_ERROR_INFO] sorted.rmdup.realign.BQSR.cvg.bed.gz not done **" && exit
fi

mv -f $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam* $final_outdir 
mv -f $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz* $final_outdir

### clear up temporary directory
rm -rf $outdir

