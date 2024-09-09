SATC="/home/genis/github/SATC/satc.R"

outprefix="/home/genis/saola/genome_masks/v2/saola_ref/satc/results/saolaRefN"

idxstats="/home/genis/saola/genome_masks/v2/saola_ref/satc/saolarefN_good_mergedups_idxstats.list"
bams="/home/genis/saola/genome_masks/v2/saola_ref/satc/saola_bams_good.list"


cat $bams | sed 's/$/.idxstats/g' > $idxstats

norm_scaff=/home/genis/saola/genome_masks/v2/saola_ref/satc/satc_norm_chroms.txt 


Rscript $SATC -i $idxstats -o $outprefix --normScaffolds $norm_scaff # --useMedian TRUE

