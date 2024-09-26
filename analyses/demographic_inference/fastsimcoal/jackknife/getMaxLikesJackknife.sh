


model=model1_prev
nsplits=50
outfile=${model}_jacknife${nsplits}.params

ml_col=$((`head -1 $model/split_1/1/1/1.bestlhoods | wc -w` - 1))
head -1 $model/split_1/1/1/1.bestlhoods > $outfile

for i in `seq 1 $nsplits`
do
    wmax=`cat $model/split_${i}/*/*/*.bestlhoods | cut -f${ml_col} | grep -v MaxEstLhood | cat -n | sort -hrk2 | cut -f1 | head -1 | tr -d '\t ' `
    tail -1 $model/split_${i}/$wmax/$wmax/$wmax.bestlhoods >> $outfile
done

