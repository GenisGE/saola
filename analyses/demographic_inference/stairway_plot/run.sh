


BLUEPRINT=/home/genis/saola/analyses/stairway_plot/saolaNorthern.mu4e-9.foldedsfs.stairwayplot.blueprint
RUNDIR=/home/genis/saola/analyses/stairway_plot/results/saolaNorthern.mu4e-9.folded  # project directory


ln -s $STAIRWAYPLOT $RUNDIR
cd $RUNDIR

java -cp stairway_plot_es Stairbuilder $BLUEPRINT
bash $BLUEPRINT.sh


