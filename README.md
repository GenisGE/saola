# Scripts used for the paper

## "Genomes of critically endangered saola are shaped by population structure and purging"

Folders in repository:

`mapping_and_filtering` contains pipeline used for filtering raw reads, mapping the samples, filtering the bam files and generating some quality control statistics.

`genome_masks` contains pipelines and scripts used for identifying regions of the genome likely to contain mapping and genotyping errors (i.e. repetitive regions, regions showing abnromal depth patterns in the mapped samples, regions with low mappability, regions showing excess of heterozygosity and sex-linked scaffolds).

`analyses` contains scripts for different analyses. They are organized by themes and within each pipeline there is a brief documentation.

`figures` contains scripts used for making the plots in the manuscript.
