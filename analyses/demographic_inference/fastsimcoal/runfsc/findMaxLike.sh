


model=model1

ml_col=$((`head -1 $model/1/1/1.bestlhoods | wc -w` - 1))
wmax=`seq 100 | xargs -n1 -I {} cat $model/{}/{}/{}.bestlhoods | cut -f${ml_col} | grep -v MaxEstLhood | cat -n | sort -hrk2 | cut -f1 | head -1 | tr -d '\t ' `

fullpath=`realpath $model/$wmax/`
ln -fs $fullpath $model/max_like_run

