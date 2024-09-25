# snakemake to run ngsrelate. config file must have:
# info: file with sample info, three columsn: 1. id 2. population 3. path to bam file
# outmain: folder to put results in
# sites: list of good regions to use in angsd sites format
# chroms: list of chromosomes to use (single column with chromosome name
# params: dictionary with minQ, minMapQ, rmTrans (wehter remove transitions, 0 no 1 yes)


import pandas as pd

ANGSD="/home/genis/software/angsd/angsd"
NGSRELATE="/home/genis/software/ngsRelate/ngsRelate"


OUTMAIN=config["outmain"]


infofile=config['info']

# this assumes info file has no header and ID is first column,
# population second column, bam file path  third column
info = pd.read_table(infofile, names=["id", "population","bams"])
pops = list(info.population.unique())

sample_size = info.population.value_counts()

rule all:
    input:
        expand(os.path.join(OUTMAIN, "ngsrelate_output", "{p}_ngsrelate.tsv"), p=pops)

        
rule do_bam_lists:
    input:
        config['info']
    output:
        expand(os.path.join(OUTMAIN, "bamlists", "bams_{p}.list"), p = pops),
        expand(os.path.join(OUTMAIN, "bamlists", "samples_{p}.list"), p = pops)
    run:
        for pop in pops:
            bams = info.loc[info.population==pop].bams
            outfile1 = os.path.join(OUTMAIN, "bamlists", "bams_{}.list".format(pop))
            ids = info.loc[info.population==pop].id
            outfile2 = os.path.join(OUTMAIN, "bamlists", "samples_{}.list".format(pop))
            bams.to_csv(outfile1, header=False, index = False)
            ids.to_csv(outfile2, header=False, index = False)
 
            
                            
rule do_gls:
    input:
        bamlist = os.path.join(OUTMAIN, "bamlists", "bams_{p}.list")
    output:
        glfile = os.path.join(OUTMAIN, "ngsrelate_input", "{p}.glf.gz"),
        ffile = os.path.join(OUTMAIN, "ngsrelate_input", "{p}.mafs.gz")
    params:
        outprefix = os.path.join(OUTMAIN, "ngsrelate_input", "{p}"),
        sites = config['sites'],
        rf = config['chroms'],
        minQ = config['params']['minQ'],
        minMapQ = config['params']['minMapQ'],
        rmTrans = config['params']['rmTrans']
    log: os.path.join(OUTMAIN, "ngsrelate_input", "{p}.arg")
    threads: 5
    shell: "{ANGSD} -b {input.bamlist} -rmTrans {params.rmTrans} -gl 2 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -minmaf 0.05 -doGlf 3 -out {params.outprefix} -P {threads} -rf {params.rf} -sites {params.sites} -minQ {params.minQ} -minMapQ {params.minMapQ}"



rule get_f:
    input:
        ffile = os.path.join(OUTMAIN, "ngsrelate_input", "{p}.mafs.gz")
    output:
        f = os.path.join(OUTMAIN, "ngsrelate_input", "{p}.freq")
    shell: "zcat {input.ffile} | cut -f5 | sed 1d > {output.f}"


        
rule run_ngsrelate:
    input:
        gl = rules.do_gls.output.glfile,
        f = rules.get_f.output.f,
        ids = os.path.join(OUTMAIN, "bamlists", "samples_{p}.list")
    output:
        out = os.path.join(OUTMAIN, "ngsrelate_output", "{p}_ngsrelate.tsv"),
    params:    
        n = lambda wildcards: sample_size[wildcards.p]
    threads: 10
    shell: "{NGSRELATE} -g {input.gl} -f {input.f} -n {params.n} -z {input.ids} -O {output.out} -p {threads}"



