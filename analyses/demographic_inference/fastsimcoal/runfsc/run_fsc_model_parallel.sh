

model=$1  ####model name
nrunsparallel=$2 ##### number of runs to do in parallel, each run uses 8 threads so total cpus used will be nrunsparallel * 8

input=../do_sfs/results/2dsfs_fscformat
nruns=100

# prepare folders to run each split
cp $input/Southern-Northern_folded_fscformat.sfs $model/_jointMAFpop1_0.obs

#cp ${model}/${model}_bestlikes.pv $model/split_${split}/${model}_bestlikes.pv


# run splits in parallel
seq 1 $nruns | xargs -n 1 -P $nrunsparallel bash run_fsc.sh $model


