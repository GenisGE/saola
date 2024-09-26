
FSC27=/maps/projects/alab/scratch/genis/software/fsc27_linux64/fsc27093

model=$1
split=$2
nits=10
ncores=8

for i in `seq 1 $nits`
  do
  mkdir $model/split_${split}/$i
  cp ${model}/$model.est $model/split_${split}/$i/$i.est
  cp ${model}/$model.tpl $model/split_${split}/$i/$i.tpl
  cp $model/split_${split}/_jointMAFpop1_0.obs $model/split_${split}/$i/${i}_jointMAFpop1_0.obs
  cp $model/split_${split}/${model}_bestlikes.pv $model/split_${split}/$i/${model}_bestlikes.pv

  export DIR=$model/split_${split}/$i
  cd $DIR
  echo "start running split $split for model $model run $i"
  $FSC27 -t $i.tpl -n500000 --initValues ${model}_bestlikes.pv -e $i.est -M -L100 -m --cores $ncores --numBatches $ncores -C100 --seed ${i}
  rm $i.est
  rm $i.tpl
  rm ${i}_jointMAFpop1_0.obs
  cd ../../..
done
