

import itertools as it
import pandas as pd

OUTMAIN=config["outmain"]

ANGSDDIR="/kellyData/home/genis/software/angsd"
ANGSD=os.path.join(ANGSDDIR, "angsd")
REALSFS=os.path.join(ANGSDDIR, "misc", "realSFS")
WINSFS="/kellyData/home/genis/.cargo/bin/winsfs"

IDS=list(config["samples"].keys())


rule do_all:
    input:
        os.path.join(OUTMAIN, "sfs", "collected.txt")



rule per_sample_saf:
    input:
        bam = lambda wildcards: config["samples"][wildcards.s]
    output:
        saf_idx = os.path.join(OUTMAIN, "safs", "{s}.saf.idx"),
        saf = os.path.join(OUTMAIN, "safs", "{s}.saf.gz"),
        saf_pos = os.path.join(OUTMAIN, "safs", "{s}.saf.pos.gz"),
    log: os.path.join(OUTMAIN, "safs", "{s}.arg")
    params:
        outprefix = lambda wildcards, output: output.saf_idx.replace(".saf.idx", ""),
        minq = 20,
        minmapq = 20,
        notrans = config["notrans"],
        mindep = lambda wildcards: int(config["depths"][wildcards.s] * 1/3),
        maxdep = lambda wildcards: int(config["depths"][wildcards.s] * 2),
        rf = config["chromlist"],
        anc = config["ref"]
    threads: 1
    shell: """
            {ANGSD} -i {input.bam} -out {params.outprefix} \
             -minQ {params.minq} -minMapQ {params.minmapq} \
              -dosaf 1 -rf {params.rf} \
              -anc {params.anc} -GL 2 -noTrans {params.notrans} \
              -setMinDepth {params.mindep} -setMaxDepth {params.maxdep} \
              -uniqueOnly 1 -doCounts 1
              """



rule sfs:
    input:
        saf_idx1 = os.path.join(OUTMAIN, "safs", "{s}.saf.idx"),
    output:
        os.path.join(OUTMAIN, "sfs", "{s}.sfs")
    threads: 15
    log:
        os.path.join(OUTMAIN, "sfs", "{s}.log")
    shell:
        "{REALSFS} -P {threads} {input.saf_idx1} > {output} 2> {log}"



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
                t = fh.readlines()[0]
                data.append([float(x) for x in t.rstrip().split()])
        a = pd.DataFrame(data, index=names, columns=["aa","ad","dd"])
        a["het"] = a["ad"] / a.sum(1)
        a.to_csv(output.f, index=True, header=True, index_label="id", sep=" ")

