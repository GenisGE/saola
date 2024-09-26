# Estiamte 1 dimensinal and 2 dimensional Site Frequency Spectrum 

Pipeline for estiamting 1d and 2d SFSs.

Tested different combination of samples with more or less stringent quality filtering of samples. On the end for northern we used northern4.bamlist.

The pipeline does SFS with both winsfs and realSFS, it was done for comparing results (at the beginning of the project winsfs was not yet developed to realSFS had been used) and on the end we used winsfs.

SFS used downstream for demographic inference with stairway plot (1dsfs) and with fastsimcoal (2dsfs) and for estimating FST.

Used ancestral state of cow for building SFS (but for the demographic inference on the end we folded the SFS). Ancestral state estiamted with pipepline in `ancestral_state` folder.