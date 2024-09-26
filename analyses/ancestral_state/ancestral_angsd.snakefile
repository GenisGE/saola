

ANGSD="/home/genis/software/angsd/angsd"


OUTMAIN="results"


           
rule all:
    input:
         os.path.join(OUTMAIN, "ancestral_cow.fa.gz")


rule do_fasta_ancestral:
    input:
        bam = config['outgroup']
    output:
        fasta = os.path.join(OUTMAIN, "ancestral_cow.fa.gz")
    params:
        rf=config['chroms'],
        sites=config['sites'],
        outprefix=os.path.join(OUTMAIN, "ancestral_cow")
    log: os.path.join(OUTMAIN, "ancestral_cow.arg")
    shell: "{ANGSD} -i {input.bam} -doCounts 1 -doFasta 2 -rf {params.rf} -sites {params.sites} -out {params.outprefix} -minMapQ 25 -minQ 30"


