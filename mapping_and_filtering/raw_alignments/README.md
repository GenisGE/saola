# Read mapping

Read mapping is performed using [PALEOMIX](https://github.com/mikkelschubert/paleomix):

    $ paleomix bam run project.yaml --max-threads ${N_THREADS} --destination output/

Some trimming and filtering of reads is performed by AdapterRemoval during the read-merging/adapter trimming step, but no filtering of mapped reads or of PCR duplicates is performed here.

## Software versions

* PALEOMIX: 6c44fa536e81476b839957e4015d0850b34e5e68
* AdapterRemoval: v2.3.2
* BWA: v0.7.17
* Picard tools: v2.24
* samtools: v1.11.0

