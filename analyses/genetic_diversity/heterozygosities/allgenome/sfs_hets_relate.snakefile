# stolen and modified from kristian snakemake in /home/leopard/users/krishang/old/relatedness/run_2dsfs_cleanref.snakefile

import itertools as it
import pandas as pd


ANGSDDIR="/home/genis/software/angsd"
ANGSD=os.path.join(ANGSDDIR, "angsd")
REALSFS=os.path.join(ANGSDDIR, "misc", "realSFS")
WINSFS="/newHome/genis/.cargo/bin/winsfs"

OUTMAIN=config["outmain"]
OUTBIG=config["outbig"]
BAMLIST=config["bamlist"]
IDLIST=config["idlist"]

with open(BAMLIST,"r") as fh:
    BAMS = [b.rstrip() for b in fh.readlines()]
with open(IDLIST, "r") as fh:
    IDS = [b.rstrip() for b in fh.readlines()]

SAMPLES = {s:b for (s,b) in zip(IDS, BAMS)}



    
rule all:
    input:
        #os.path.join(OUTMAIN, "sfs_2d", "collected.txt"),
        os.path.join(OUTMAIN, "sfs", "collected.txt"),
    #    os.path.join(OUTMAIN,"finished2DSts.txt")


rule all_rel:
    input:
        os.path.join(OUTMAIN, "sfs_2d", "collected.txt")
    
rule all_saf_hets:
     input:
        os.path.join(OUTMAIN, "sfs", "collected.txt")

        
rule per_sample_saf:
    input:
        bam = lambda wildcards: SAMPLES[wildcards.s]
    output:
        saf_idx = os.path.join(OUTBIG, "safs", "{s}.saf.idx"),
        saf = os.path.join(OUTBIG, "safs", "{s}.saf.gz"),
        saf_pos = os.path.join(OUTBIG, "safs", "{s}.saf.pos.gz"),
    log: os.path.join(OUTBIG, "safs", "{s}.arg")
    params:
        outprefix = lambda wildcards, output: output.saf_idx.replace(".saf.idx", ""),
        minq = config["filters"]["minq"],
        minmapq = config["filters"]["minmapq"],
        notrans = config["filters"]["notrans"],
        sites = config["mask"],
        rf = config["chromlist"],
        anc = config["ref"]
    threads: 1
    shell:
        "{ANGSD} -i {input.bam} -out {params.outprefix} -minQ {params.minq} -minMapQ {params.minmapq} -dosaf 1 -sites {params.sites} -rf {params.rf} -anc {params.anc} -GL 2 -noTrans {params.notrans}"


        
rule sfs:
    input:
        saf_idx1 = os.path.join(OUTBIG, "safs", "{s}.saf.idx"),
    output:
        os.path.join(OUTMAIN, "sfs", "{s}.sfs")
    threads: 15
    log:
        os.path.join(OUTMAIN, "sfs", "{s}.log")
    shell:
        "{WINSFS} -t {threads} {input.saf_idx1} > {output} 2> {log}"



rule collect_sfs:
    input:
        expand(os.path.join(OUTMAIN, "sfs", "{s}.sfs"), s=IDS)
    output:
        f=os.path.join(OUTMAIN, "sfs", "collected.txt")
    run:
        import os
        import pandas as pd
        data = []
        names = []
        for x in input:
            name = os.path.basename(x).replace(".sfs", "")
            names.append(name)
            with open(x, 'r') as fh:
                t = fh.readlines()[1]
                data.append([float(x) for x in t.rstrip().split()])
        a = pd.DataFrame(data, index=names, columns=["aa","ad","dd"])
        a["het"] = a["ad"] / a.sum(1)
        a.to_csv(output.f, index=True, header=True, index_label="id", sep=" ")



rule sfs_2d:
    input:
        saf_idx1 = os.path.join(OUTBIG, "safs", "{s1}.saf.idx"),
        saf_idx2 = os.path.join(OUTBIG, "safs", "{s2}.saf.idx"),
    output:
        os.path.join(OUTMAIN, "sfs_2d", "{s1}_{s2}.sfs")
    threads: 15
    log:
        os.path.join(OUTMAIN, "sfs_2d", "{s1}_{s2}.pre_log")
    shell:
        "{WINSFS} -t {threads} {input.saf_idx1} {input.saf_idx2} > {output} 2> {log}"

        

rule collect_2dsfs:
    input:
        expand(os.path.join(OUTMAIN, "sfs_2d", "{s[0]}_{s[1]}.sfs"), s=it.combinations(IDS, 2)),
    output:
        f=os.path.join(OUTMAIN, "sfs_2d", "collected.txt")
    run:
        data = []
        names = []
        for x in input:
            name = os.path.basename(x).replace(".sfs", "")
            names.append(name)
            with open(x, 'r') as fh:
                t = fh.readlines()[1]
                data.append([float(x) for x in t.rstrip().split()])
        a = pd.DataFrame(data, index=names, columns=["aaAA","aaAD","aaDD","adAA","adAD","adDD","ddAA","ddAD","ddDD"])

        a["r0"] = (a["aaDD"]+a["ddAA"])/a["adAD"]
        a["r1"] = a["adAD"] / (a.iloc[:,[1,2,3,5,6,7]].sum(1))
        a["king"] = (a["adAD"] - 2*a[["aaDD", "ddAA"]].sum(1)) / (a.iloc[:,[1,3,5,7]].sum(1) + 2*a["adAD"])
        a.to_csv(output.f, index=True, header=True, index_label="id", sep=" ")


