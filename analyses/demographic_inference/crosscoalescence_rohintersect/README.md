# Pipeline for inference of cross coalescence rates using two unphased inbred genomes

Infer population crosscoalescence rates betweent the two popualtions. We have too low sample size to even attempt phasing, but we took advantadge that saolas are qutie inbred and has thus quite a lot of Runs of Homozygosity (ROH). ROH have a single haplotype so it is already phased. Then if two samples are in roh in the same genomic region, we have one perfectly phased haplotype for each sample in that area.

Because we had used psmc for per sample demographic inference, we also used psmc here. We build 'pseudodiploid' sample when two samples from each popualtion are in roh. Then within population coalescence are done using a comparable region of the genome. As shown in the mansucript, rohs are not totally randomly distributed across the genome, so we want to run the per sample psmc in regions where there are roh in other sampels from the popualtion but not in that specific sample. Also the crosscoalescence is estimated in a reduced region of the genome (~ 100 Mbp), so we want the within population coalescence to be estimated in region of the same size.

See mansucript methods and pipeline for details on how it is done. Pipeline might seem a bit messy because it was built by trial and error and it has not been cleaned. But it is what was used for the mansucript.

ROH used are estiamted in `analyses/genetic_diversity/rohs`
