# snakemake to estimate 2dsfs, sfs and fst between popualtions pairs.
# version 2 intention is to do a more efficient version which will use downsampled saf
# reminder to explain what config needs. also put angsd parameters to config would be nice


import pandas as pd
import itertools as it

ANGSDDIR="/home/genis/github/angsd"
ANGSD=os.path.join(ANGSDDIR, "angsd")
REALSFS=os.path.join(ANGSDDIR, "misc/realSFS")
WINSFS="/newHome/genis/.cargo/bin/winsfs"
RSCRIPT="Rscript"
PLOT1DSFS="scripts/plotsfs.R"
PLOT2DSFS="scripts/plot2dsfs.R"

OUTMAIN=config["outmain"]
OUTBIG=config["outbig"]


wildcard_constraints:
    p = "|".join(config["bamlists"].keys()),
    p1 = "|".join(config["bamlists"].keys()),
    p2 = "|".join(config["bamlists"].keys()),



rule all:
    input:
        expand(os.path.join(OUTMAIN, "sfs", "{p}_{f}.png"),
               p=config["bamlists"].keys(),
               f=["folded", "unfolded"]),
        expand(os.path.join(OUTMAIN, "sfs", "{p}.sfs"), p=config["bamlists"].keys()),
        expand(os.path.join(OUTMAIN, "2dsfs", "{p}{f}"),
                p= ["Central_Northern", "Central_Northern2", "Central_Northern3", "9176_9264"],
                f=["_marginalized1dSFS_unfolded.png", "_2DSFS_folded.png", "_2DSFS_unfolded.png"]),
        os.path.join(OUTMAIN, "2dsfs", "Central_Northern.sfs"),
        os.path.join(OUTMAIN, "2dsfs", "Central_Northern2.sfs"),
        os.path.join(OUTMAIN, "2dsfs", "Central_Northern3.sfs"),
        os.path.join(OUTMAIN, "2dsfs", "9176_9264.sfs")



rule do_saf:
    input:
        bamlist = lambda wildcards: config["bamlists"][wildcards.p]
    output:
        saf = os.path.join(OUTBIG, "safs", "{p}.saf.gz"),
        saf_idx = os.path.join(OUTBIG, "safs", "{p}.saf.idx"),
        saf_pos = os.path.join(OUTBIG, "safs", "{p}.saf.pos.gz")
    params:
        anc = config['anc'],
        outprefix = os.path.join(OUTBIG, "safs", "{p}"),
        sites = config['sites'],
        rf = config['chroms'],
        minq = config["filters"]["angsd"]["minq"],
        minmapq = config["filters"]["angsd"]["minmapq"],
        notrans = config["filters"]["angsd"]["notrans"],
    log: os.path.join(OUTBIG, "safs", "{p}.arg")
    threads: 5
    shell: "{ANGSD} -b {input.bamlist} -gl 2 -dosaf 5 -anc {params.anc} -out {params.outprefix} -P {threads} -rf {params.rf} -sites {params.sites} -minQ {params.minq} -minMapQ {params.minmapq} -noTrans {params.notrans} -doMajorMinor 5"



rule do_sfs:
    input:
        idx = os.path.join(OUTBIG, "safs", "{p}.saf.idx")
    output:
        sfs = os.path.join(OUTMAIN, "sfs", "{p}.sfs")
    threads: 15
    params:
        w = config["winsfs_params"]["w"]
    log: os.path.join(OUTMAIN, "sfs", "{p}.log")
    shell: "{WINSFS} -vvv -w {params.w} {input.idx} -t {threads} > {output.sfs} 2> {log}"


    
rule do_sfs_realsfs:
    input:
        saf = os.path.join(OUTBIG, "safs", "{p}.saf.idx")
    output:
        sfs = os.path.join(OUTMAIN, "sfs", "{p}_realsfs.sfs")
    log: os.path.join(OUTMAIN, "sfs", "{p}_realsfs.log")
    threads: 15
    shell: """
    {REALSFS} {input.saf} -P {threads} > {output.sfs} 2> {log}
"""



rule do_2dsfs:
    input:
        idx1 = os.path.join(OUTBIG, "safs", "{p1}.saf.idx"),
        idx2 = os.path.join(OUTBIG, "safs", "{p2}.saf.idx")
    output:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{p1}_{p2}.sfs")
    log: os.path.join(OUTMAIN, "2dsfs", "{p1}_{p2}.log")
    threads: 15
    params:
        w = config["winsfs_params"]["w"]
    shell: "{WINSFS} -vvv -w {params.w}  {input.idx1} {input.idx2} -t {threads} > {output.sfs} 2> {log}"



rule plot_sfs:
    input:
        sfs =  os.path.join(OUTMAIN, "sfs", "{p}.sfs")
    output:
        expand(os.path.join(OUTMAIN, "sfs", "{{p}}_{f}.png"),
               f=["folded", "unfolded"])
    params:
        outprefix = os.path.join(OUTMAIN, "sfs", "{p}")
    shell: """
    {RSCRIPT} {PLOT1DSFS} {input.sfs} {params.outprefix}
"""


rule plot_2dsfs:
    input:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{p1}_{p2}.sfs")
    output:
        expand(os.path.join(OUTMAIN, "2dsfs", "{{p1}}_{{p2}}{f}"),
               f = ["_marginalized1dSFS_unfolded.png",
                    "_2DSFS_folded.png",
                    "_2DSFS_unfolded.png"]
        )
    params:
        outprefix = os.path.join(OUTMAIN, "2dsfs", "{p1}_{p2}"),
    shell: """
    {RSCRIPT} {PLOT2DSFS} {input.sfs} {params.outprefix}
"""


