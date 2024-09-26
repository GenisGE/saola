# estiamte population size trajectories through time with the pairwise sequentially markovian coalescnet (psmc)

Does psmc and also estiamtes genome wide heterozygositeis from the called genotypes.

The snakemake has rules to run psmc with and without transitons, because some samples had excess errors in transitions due to DNA damage.

The pipeline starts from the bam files and calls genotyeps, cleans them up and then runs psmc.
