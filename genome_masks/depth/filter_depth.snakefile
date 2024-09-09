# snakemake to do depth filters: from group/s of bam files generates bed files with list of regions to keep after excluding regions with very low or very high local depth


# config needs:
# groups: group: bamlist (dict where keys are groups (e.g. lowdepth/highdepth) and values list of samples corresponding to that group.
# chroms: list of chromosomes to do depths on
# angsdparams: minQ, minMapQ, maxdepth
# outmain: path to main output folder
# all_bed: path to bed file with all sites (i.e. for each chromosome included "chr start end" where start is 0 and end is lenght of chromosome)

ANGSD="/home/genis/software/angsd/angsd"
BEDTOOLS="/home/genis/software/bedtools2/bin/bedtools"

R="Rscript"
PYTHON="python3"

COMBINE="scripts/combineDepths.R"
PLOT="scripts/plotDepths.R"
FILTER="scripts/doDepthFilter.py"

OUTMAIN = config["outmain"]

wildcard_constraints:
    group="|".join(config["groups"].keys())

with open(config["chroms"]) as fh:
    CHROMS = [x.rstrip() for x in fh.readlines()]

rule all:
    input:
        os.path.join(OUTMAIN, "plots", "depthFilter.png"),
        expand(os.path.join(OUTMAIN, "bed", "all_{w}.bed"), w = ["keep", "remove"])

        
rule do_depths:
    input:
        bamlist = lambda wildcards: config["groups"][wildcards.group]
    output:
        global_depth = os.path.join(OUTMAIN, "depths", "{group}", "{chrr}_{group}.depthGlobal"),
        pos = temp(os.path.join(OUTMAIN, "depths","{group}", "{chrr}_{group}.pos.gz")),
    params:
        outprefix = lambda wildcards, output: output.pos.replace(".pos.gz",""),
        maxdepth = config["angsdparams"]["maxdepth"],
        minQ = config["angsdparams"]["minQ"],
        minMapQ = config["angsdparams"]["minMapQ"],
        r = "{chrr}"
    log: os.path.join(OUTMAIN, "depths","{group}", "{chrr}_{group}.arg"),
    shell: "{ANGSD} -doCounts 1 -doDepth 1 -dumpCounts 1 -maxdepth {params.maxdepth} -minQ {params.minQ} -minMapQ {params.minMapQ} -r {params.r} -bam {input.bamlist} -out {params.outprefix}"



rule combine_depths:
    input:
        global_depths = expand(os.path.join(OUTMAIN, "depths", "{{group}}", "{chrr}_{{group}}.depthGlobal"), chrr=CHROMS)
    output:
        dist = os.path.join(OUTMAIN, "depths", "all_{group}.depthGlobal"),
        median = os.path.join(OUTMAIN, "depths", "{group}.depthMedian")
    params:
        group = "{group}",
        chroms = config["chroms"],
        maxdepth = config["angsdparams"]["maxdepth"],
        outmain = OUTMAIN
    shell: "{R} {COMBINE} {params.group} {params.chroms} {params.maxdepth} {params.outmain}"


rule plot_depths:
    input:
        dists = expand(os.path.join(OUTMAIN, "depths", "all_{group}.depthGlobal"), group = config["groups"].keys()),
        medians = expand(os.path.join(OUTMAIN, "depths", "{group}.depthMedian"), group = config["groups"].keys())
    output:
        os.path.join(OUTMAIN, "plots", "depthFilter.png")
    params:
        outmain = OUTMAIN,
        groups = " ".join([x for x in config["groups"].keys()])
    shell: "{R} {PLOT} {params.outmain} {params.groups}"

           
rule do_filter:
    input:
        pos = expand(os.path.join(OUTMAIN, "depths", "{{group}}", "{chrr}_{{group}}.pos.gz"), chrr=CHROMS),
        median = os.path.join(OUTMAIN, "depths", "{group}.depthMedian")
    output:
        bed = os.path.join(OUTMAIN, "bed", "{group}_remove.bed")
    log:
        os.path.join(OUTMAIN, "bed", "{group}_remove.log")
    params:
        chroms = config["chroms"],
        indir = os.path.join(OUTMAIN, "depths", "{group}"),
        outprefix = os.path.join(OUTMAIN, "bed", "{group}_remove")
    shell: "{PYTHON} {FILTER} {params.indir} {input.median} {params.chroms} {params.outprefix}"


rule combine_beds:
    input:
        beds = expand(os.path.join(OUTMAIN, "bed", "{group}_remove.bed"), group=config["groups"].keys())
    output:
        r = os.path.join(OUTMAIN, "bed", "all_remove.bed"),
        k = os.path.join(OUTMAIN, "bed", "all_keep.bed")
    params:
        allbed = config["allbed"]
    shell: """
    cat {input.beds} | {BEDTOOLS} sort | {BEDTOOLS} merge > {output.r}
    {BEDTOOLS} subtract -a {params.allbed} -b {output.r} > {output.k}
"""

