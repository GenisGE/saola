
NGSA=/home/genis/saola/analyses/ngsadmix/scripts/NGSadmixConv_all26saolas.sh


echo `seq 2 4` | xargs -n1 -P3 $NGSA

# for clean up
#rm /home/genis/impala/analyses/impalaMap/ngsadmix/results/*/admixResultImpalaMap.*fopt.gz
