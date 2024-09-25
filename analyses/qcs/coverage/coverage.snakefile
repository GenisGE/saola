# snakemake to estiamte coverage after filtering per sample.
# configfile requrieds:
#    bed: path to bedfile with regions to keep
#    outmain: path to main folder
#    samples: dic where keys:values is sample_name:pathtobamfile

SAMTOOLS="/kellyData/home/genis/software/samtools-1.9/samtools"
AWK_STUFF = "a.awk"


SAMPLES = [x for x in config["samples"].keys()]
OUTMAIN=config["outmain"]

rule all:
     input:
         expand(os.path.join(OUTMAIN, "samples", "{s}.depth"), s=SAMPLES),
         os.path.join(OUTMAIN, 'all_depths.txt')


rule do_depth_samtools:
     input:
        bam= lambda wildcards: config["samples"][wildcards.s]
     output:
        os.path.join(OUTMAIN, "samples", "{s}.depth")
     params:
        bed = config["bed"],
        minmapq = config["minmapq"]
     shell: "{SAMTOOLS} depth -a -Q {params.minmapq} -b {params.bed} {input.bam} | awk -f {AWK_STUFF} > {output}"


rule collect_res:
     input:
        expand(os.path.join(OUTMAIN, "samples", "{s}.depth"), s=SAMPLES)
     output:
        os.path.join(OUTMAIN, "all_depths.txt")
     shell: "cat {input} > {output}"
