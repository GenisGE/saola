# Filtering and SAMTools statistics

Run the SnakeMake file to filter the BAMs and to collect SAMTools statistics:

    $ snakemake -np
    $ snakemake -p -j ${N_THREADS}

Filters set in the `Snakefile` file, under the `filter_bam` rule.

