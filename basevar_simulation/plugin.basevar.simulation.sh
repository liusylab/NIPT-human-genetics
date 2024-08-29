#---
#do simulation
#---
script=`pwd`/bin/multiAllelicSimulation.R
date=20170427
ind=1000
totalrep=140
dep=0.06
error=0.01

s=1
e=150100

prefix=simulation.ind$ind.start$s.end$e
shell=$prefix.ind$totalrep.bc.q.sh
>$shell
outdir=`pwd`/pileup_sim_new/start$s.end$e
mkdir -p $outdir

for rep in {1..140}
do
    outfile=$outdir/$prefix.rep$rep.pileup
    echo "time Rscript $script $ind $dep $error $s $e $outfile && bgzip -c $outfile >$outfile.gz && tabix -f -b 2 -e 2 $outfile.gz && rm -rf $outfile" >>$shell

    rep2=$((rep-1))
    i_start=`expr $rep2 \* $ind`
    i_start=$((i_start+1))
    i_end=`expr $rep \* $ind`
    samplename=$outdir/$prefix.rep$rep.samname

    for i in $(seq $i_start $i_end);
    do
        echo "ind$i"
    done >$samplename
done

echo "time perl /home/share/software/com_extra/qsubtools/8.4.2/bin/qsub-sge --queue bc.q --lab b2c_rd --resource vf=2G --lines 1 -maxjob 500 --convert no $shell && echo \"** qsub $shell done **\"" >qsub-sge.$shell.sh
echo "nohup sh qsub-sge.$shell.sh 2>qsub-sge.$shell.log2 &" >nohup.qsub-sge.$shell.sh

