# Estiamte heterozygosity

Estiamte heterozygosities from genotype likelihoods. First do saf file per sample, then estiamte per sample site frequency spectrum with winsfs.

Two runs, one using all genome in `allgenome` and another excluding regions where each sample has long roh (> 1 Mbp) in `norohs`. Roh used are the ones inferred with the pipleine in `analyses/genetic_diversity/rohs`

Within each of these there is one run using all mutations and another excludign transition mutations for those samples with excess transition mutations due to DNA damage.


## Two other heterozygosity pipelines based on different papers, for comparison of saola heterozygosities with other species

- In folder `hets_asrhino` pipeline as in Liu et al. 2021 Ancient and modern genomes unravel the evolutionary history of the rhinoceros family, Cell, 184:19 4874-4885

- In folder `hetas_aszoonomia` piepeline as in Zoonomia Consortium 2020 A comparative genomics multitool for scientific discovery and conservation, Nature 587, 240-245

See manuscript for details.