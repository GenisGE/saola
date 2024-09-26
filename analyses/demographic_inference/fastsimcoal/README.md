# Fit 2 popualtion demographic model from the 2dSFS with fascimcoal

First infer the sfs with pipeline in `do_sfs` folder. The inference should be the same as in `analyses/sfs`, but in this case we estiamte 50 splits of block sfs for the jackknife and also format it for fastsimcoal.

Then with scripts in `runfsc` folder find the maximum likelihood the parameters of the demographic model fitted to the whole genome 2dsfs.

Finally with scripts in `jackknife` folder fit parameters to 50 leave-one-block-out 2dsfs and then use the estiamte to get standard errors with jackknfie estiamtor.