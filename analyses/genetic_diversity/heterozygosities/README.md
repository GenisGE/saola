# Estiamte heterozygosity

Estiamte heterozygosities from genotype likelihoods. First do saf file per sample, then estiamte per sample site frequency spectrum with winsfs.

Two runs, one using all genome in `allgenome` and another excluding regions where each sample has long roh (> 1 Mbp) in `norohs`. Roh used are the ones inferred with the pipleine in `analyses/genetic_diversity/rohs`

Within each of these there is one run using all mutations and another excludign transition mutations for those samples with excess transition mutations due to DNA damage.
