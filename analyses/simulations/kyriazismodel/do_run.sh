module load SLiM/4.0.1

nrun=$1

slimscript=saola_double_kyriasiz_nwf.slim

mkdir -p $nrun
cd $nrun
cp ../$slimscript .
runslim=`basename $slimscript`

echo "start running $nrun"
slim $runslim > nwf_double_kyriazis_${nrun}.log

