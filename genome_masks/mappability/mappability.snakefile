# snakemake to do mappability in reference genome.
# config requries:
# ref: path to reference to do mappability on
# outmain: folder to write results to
# tmp_dir: path to tmp folder (place that can temporarily store big things)
# settings: k: lenght of reads. e: mismatches allowed. m: minimum mappability to keep.


GENMAP="/home/genis/software/genmap-build/bin/genmap"

OUTMAIN=config["outmain"]

K=config["settings"]["k"]
E=config["settings"]["e"]
M=config["settings"]["m"]

rule all:
    input:
        bed = os.path.join(OUTMAIN, "bed", "mappability_m{m}_k{k}_e{e}.bed".format(m=M, k=K, e=E))

rule index_ref:
    input:
        ref = config["ref"]
    output:
        done = os.path.join(OUTMAIN, "index_done.txt")
        #index = os.path.join(OUTMAIN, "index", "index.txt.concat")
    params:
        index = os.path.join(OUTMAIN, "index"),
        tmp_dir = config["tmp_dir"]
    shell: """
    export TMPDIR={params.tmp_dir}
    {GENMAP} index -F {input.ref} -I {params.index} -A skew
    touch {output.done}
"""
    



rule mappability:
    input:
        index = rules.index_ref.output.done
    output:
        bedgraph = os.path.join(OUTMAIN, "mappability", "mappability_K{k}_E{e}.bedgraph".format(k=K, e=E))
    params:
        outprefix = lambda wildacrds, output:  output.bedgraph.replace(".bedgraph", ""),
        index = os.path.join(OUTMAIN, "index"),
        k = K,
        e = E
    threads: 10
    shell:
        "{GENMAP} map -K {params.k} -E {params.e} -I {params.index} -O {params.outprefix} --bedgraph -T {threads}"

        
            
rule make_bed:
    input:
        rules.mappability.output.bedgraph,
    output:
        bed = os.path.join(OUTMAIN, "bed", "mappability_m{m}_k{k}_e{e}.bed".format(m=M, k=K, e=E))
    params:
        m = M
    shell: """
    awk '$4 >= {M} {{print $1"\t"$2"\t"$3}}' {input} > {output.bed} 
"""
