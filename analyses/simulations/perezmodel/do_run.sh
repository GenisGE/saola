module load SLiM/4.0.1

nrun=$1

slimscript=saola_double_perez_nwf.slim

mkdir -p $nrun
cd $nrun
cp ../$slimscript .
runslim=`basename $slimscript`

echo "start running $nrun"
slim $runslim > nwf_double_perez_${nrun}.log

