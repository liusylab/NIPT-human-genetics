#!/bin/sh

while getopts "f:c:l:s:a:t:g:j:r:i:b:o:u:e:z:x:" opt
do
	case $opt in 
		###INPUT
		f)
		fq=$OPTARG
		;;
		c)
		cell_id=$OPTARG
		;;
		l)
		lane_id=$OPTARG
		;;
		s)
		sample_id=$OPTARG
		;;
		###tools
		a)
		bwa=$OPTARG
		;;
		t)
		samtools=$OPTARG
		;;
		g)
		gatk=$OPTARG
		;;
		j)
		java=$OPTARG
		;;
		e)
		bedtools=$OPTARG
		;;
		z)
		bgzip=$OPTARG
		;;
		x)
		tabix=$OPTARG
		;;
		### Reference
		r)
		ref=$OPTARG
		;;
		i)
		ref_index_prefix=$OPTARG
		;;
		### GATK bundle
		b)
		gatk_bundle_dir=$OPTARG
		;;
		### output directory
		o)
		outdir=$OPTARG
		;;
		u)
		final_outdir=$OPTARG
		;;
	esac
done

echo "parameters:-fq $fq -cid $cell_id -lid $sample_id -sid $sample_id -bwa $bwa -stools $samtools -gatk $gatk -java $java -ref $ref -refindex $ref_index_prefix -gatkb $gatk_bundle_dir -tout $outdir -fout $final_outdir"

if [ ! -d $final_outdir ]
then mkdir -p $final_outdir
fi

if [ ! -d $outdir ]
then mkdir -p $outdir
fi

######################################################################################
################################### Pipeline ########################################
######################################################################################
echo "We're doing the job in $cell_id, $lane_id and $sample_id"
echo "We are calculating $fq"
echo "We'll save it in $final_outdir"
### step 1: Adjust the number of threads (-t) for bwa according to your cluster. The -e 10 parameter makes indel detection more sensitive, though its impact is minimal.
echo ""
time $bwa aln -e 10 -t 8 -i 5 -q 0 $ref_index_prefix $fq > $outdir/${sample_id}.sai && \
    time $bwa samse -r "@RG\tID:${cell_id}_${lane_id}\tPL:COMPLETE\tSM:${sample_id}" $ref_index_prefix $outdir/${sample_id}.sai $fq | $samtools view -h -Sb - > $outdir/${sample_id}.bam && echo "** bwa done **" && \
    time $samtools sort -@ 8 -O bam -o $outdir/${sample_id}.sorted.bam $outdir/${sample_id}.bam && echo "** bam sorted done **" && \
    time $samtools rmdup $outdir/${sample_id}.sorted.bam $outdir/${sample_id}.sorted.rmdup.bam && echo "** rmdup done **" && \
    time $samtools index $outdir/${sample_id}.sorted.rmdup.bam && echo "** index done **" && touch ${outdir}/bwa_sort_rmdup.finish

if [ ! -f ${outdir}/bwa_sort_rmdup.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bwa_sort_rmdup not done **" && exit
fi

### step 2: realignment, adjust the memory usage according to data amount
### Realign
echo ""
time $java -Xmx15g -jar $gatk \
    -T RealignerTargetCreator \
    -R $ref \
    -I $outdir/${sample_id}.sorted.rmdup.bam \
    -known $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -known $gatk_bundle_dir/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -o $outdir/${sample_id}.indel_target_intervals.list && echo "** RealignerTargetCreator done " && touch ${outdir}/RealignerTargetCreator.finish

if [ ! -f ${outdir}/RealignerTargetCreator.finish ]
then echo "** [WORKFLOW_ERROR_INFO] RealignerTargetCreator not done **" && exit
fi

time $java -Xmx15g -jar $gatk \
    -T IndelRealigner \
    -R $ref \
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
    -R $ref \
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
    -R $ref \
    --BQSR $outdir/${sample_id}.recal_data.table \
    -I $outdir/${sample_id}.sorted.rmdup.realign.bam \
    -o $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam && echo "** PrintReads done **" && touch ${outdir}/PrintReads.finish

if [ ! -f ${outdir}/PrintReads.finish ]
then echo "** [WORKFLOW_ERROR_INFO] PrintReads not done **" && exit
fi

time $samtools index $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam && echo "** bam index done **" && touch ${outdir}/bam_index.finish

if [ ! -f ${outdir}/bam_index.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bam index not done **" && exit
fi

:<<'BLOCK'
### Convert Bam to Cram to save storage
###Since variant detection appears to be particularly slow when reading cram files, and the cause hasn't been identified, the workflow has been adjusted to perform variant detection first, followed by the conversion from bam to cram.
echo ""
time samtools view -C -T $ref $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam > $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cram && echo "** cram done **" && touch ${outdir}/cram.finish
if [ ! -f ${outdir}/cram.finish ]
then echo "** [WORKFLOW_ERROR_INFO] cram not done **" && exit
fi

#time samtools index $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cram && echo "** cram index done **" && touch ${outdir}/cram_index.finish

#if [ ! -f ${outdir}/cram_index.finish ]
#then echo "** [WORKFLOW_ERROR_INFO] cram index not done **" && exit
#fi
BLOCK

### step 5. bam stats
time $samtools stats $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam > $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bamstats && echo "** bamstats done **" && touch ${outdir}/bamstats.finish

if [ ! -f ${outdir}/bamstats.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bamstats not done **" && exit
fi

### step 6. bedtools
time $bedtools genomecov -ibam $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam -bga -split | $bgzip > $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz && $tabix -p bed $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz && echo "** sorted.rmdup.realign.BQSR.cvg.bed.gz done **" && touch ${outdir}/sorted_rmdup_realign_BQSR_cvg_bed.finish

if [ ! -f ${outdir}/sorted_rmdup_realign_BQSR_cvg_bed.finish ]
then echo "** [WORKFLOW_ERROR_INFO] sorted.rmdup.realign.BQSR.cvg.bed.gz not done **" && exit
fi

mv -f $outdir/${sample_id}.sorted.rmdup.realign.BQSR.bam* $final_outdir 
mv -f $outdir/${sample_id}.sorted.rmdup.realign.BQSR.cvg.bed.gz* $final_outdir


