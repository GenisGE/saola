EVALADMIX=/home/genis/github/evalAdmix/evalAdmix

bgl=/home/genis/saola/analyses/beagle_pcangsd/results_mergedseq/beagle/all.beagle.gz

resdir=/home/genis/saola/analyses/ngsadmix/results/
outdir=/home/genis/saola/analyses/ngsadmix/results/evaladmix
outname=evalAdmixResultsSaola

#for k in `seq 2 3`
#do
#    qfile=`find $resdir/$k | grep qopt_conv`
#    ffile=`find $resdir/$k | grep fopt_conv`
#    out=$outdir/${outname}_K${k}
#    $EVALADMIX -beagle $bgl -qname $qfile -fname $ffile -P 20 -o $out.corres &> $out.log
#done

bgl=/home/genis/saola/analyses/beagle_pcangsd/results_mergedseq/beagle/all.beagle.gz

resdir=/home/genis/saola/analyses/ngsadmix/results_all26
outdir=/home/genis/saola/analyses/ngsadmix/results_all26/evaladmix
outname=evalAdmixResultsSaola_all26
k=2

qfile=`find $resdir/$k | grep qopt_conv`
ffile=`find $resdir/$k | grep fopt_conv`
out=$outdir/${outname}_K${k}
#$EVALADMIX -beagle $bgl -qname $qfile -fname $ffile -P 20 -o $out.corres &> $out.log

for k in 3 4
do
    

    qfile=`find $resdir/$k | grep qopt_conv`
    ffile=`find $resdir/$k | grep fopt_conv`
    out=$outdir/${outname}_K${k}
    $EVALADMIX -beagle $bgl -qname $qfile -fname $ffile -P 20 -o $out.corres &> $out.log
done

