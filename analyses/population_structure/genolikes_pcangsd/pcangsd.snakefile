"""
snakemake to estimate genotype likelihoods and run pcangsd
configfile needs:
       - info: file with bamfile in column1 and population assignment in column2 path to bamfile. Population is just for visualization, if unknown can be set to arbitrary value for all
       - chroms: path to file with list of chromosomes to use
       - mask: file in sites angsd format with good regions or sites to use
       - ref: path to reference
       - outmain: path to main output
       - params:
                - angsd: minQ: value
                         minMapQ: value
                         rmtrans: 0/1 (no/yes)
                - pcangsd: e: number principal components to use, if dont' know set to 0 and pcangsd infers it
      - scriptsdir: path to directory with requried scripts
"""

ANGSD="/home/genis/software/angsd/angsd"
PCANGSD="/home/genis/software/pcangsd-v.0.99/pcangsd.py"
PYTHON="python3"
R="Rscript"
BEDTOOLS="/home/genis/software/bedtools2/bin/bedtools"
SAMTOOLS="/home/genis/software/samtools-1.9/samtools"

SCRIPTSDIR = config["scriptsdir"]
PLOTPCA = os.path.join(SCRIPTSDIR,"plotPCA.R")

OUTMAIN=config["outmain"]


CHROMLIST=config["chroms"]
with open(CHROMLIST, "r") as fh:
    CHROMS=[x.rstrip() for x in fh.readlines()]

REF = config["ref"]

# LOAD PARAMS
MINQ=config["params"]["angsd"]["minQ"]
MINMAPQ=config["params"]["angsd"]["minMapQ"]
E=config["params"]["pcangsd"]["e"]


rule all:
    input:
        os.path.join(OUTMAIN, "pcangsd", "pca.png")
        #os.path.join(OUTMAIN, "something")


rule do_bamlist:
    """extract path to bam files from info file"""
    input:
        info = config["info"]
    output:
        bamlist = os.path.join(OUTMAIN, "info", "bams.list")
    shell:"""
    cut -f2 {input.info} > {output.bamlist}
"""


rule do_poplist:
    """extract population assignment from info file"""
    input:
        info = config["info"]
    output:
        poplist = os.path.join(OUTMAIN, "info", "pop.list")
    shell:"""
    cut -f1 {input.info} > {output.poplist}
"""



rule do_beagle_perchrom:
    """estiamte genotype likelihoods per chromosome"""
    input:
        bamlist = os.path.join(OUTMAIN, "info", "bams.list")
    output:
        bgl = temp(os.path.join(OUTMAIN, "beagle", "{chrom}.beagle.gz")),
        maf = temp(os.path.join(OUTMAIN, "beagle", "{chrom}.mafs.gz"))
    log: os.path.join(OUTMAIN, "beagle", "{chrom}.arg")
    params:
        outprefix = os.path.join(OUTMAIN, "beagle", "{chrom}"),
        r = "{chrom}",
        sites=config["mask"],
        rmtrans = config["params"]["angsd"]["rmtrans"]
    shell:"""
    {ANGSD} -GL 2 -out {params.outprefix} -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam {input.bamlist} -minmapQ {MINMAPQ} -minQ {MINQ} -r {params.r} -minmaf 0.05 -rmTrans {params.rmtrans} -sites {params.sites}
"""


rule concat_beagle:
    """put together genotype likelihood chromosome fiels into a single file"""
    input:
        expand(os.path.join(OUTMAIN, "beagle", "{chrom}.beagle.gz"), chrom=CHROMS)
    output:
        bgl = os.path.join(OUTMAIN, "beagle", "all.beagle.gz")
    params:
        indir=os.path.join(OUTMAIN, "beagle"),
        outprefix = os.path.join(OUTMAIN, "beagle", "all")
    log: os.path.join(OUTMAIN, "beagle", "all.arg")
    run:
        c = CHROMS[0]
        shell("zcat {params.indir}/{c}.beagle.gz | awk 'NR==1' > {params.outprefix}.beagle") # need awk shit bc head -1 makes snakemake fail
        for c in CHROMS:
            shell("zcat {params.indir}/{c}.beagle.gz | sed 1d >> {params.outprefix}.beagle")
            shell("cat {params.indir}/{c}.arg >> {params.outprefix}.arg")
        shell("gzip {params.outprefix}.beagle")



rule do_pcangsd:
    """run pcangsd to estimate covariance matrix"""
    input:
        beagle = os.path.join(OUTMAIN, "beagle", "all.beagle.gz")
    output:
        os.path.join(OUTMAIN, "pcangsd", "out.cov"),
    params:
        outprefix = os.path.join(OUTMAIN, "pcangsd", "out"),
        e = E
    log: os.path.join(OUTMAIN, "pcangsd", "out.log")
    threads: 20
    shell:"""
    if [ {params.e} = 0 ]
    then
        {PYTHON} {PCANGSD} -beagle {input.beagle} -o {params.outprefix} -threads {threads} > {log}
    else
        {PYTHON} {PCANGSD} -beagle {input.beagle} -o {params.outprefix} -threads {threads} -e {params.e} > {log}
    fi
"""



rule plot_pca:
    input:
        cov = os.path.join(OUTMAIN, "pcangsd", "out.cov"),
        poplist = os.path.join(OUTMAIN, "info", "pop.list")
    output:
        pcaplot = os.path.join(OUTMAIN, "pcangsd", "pca.png")
    shell:"""
    {R} {PLOTPCA} {input.cov} {input.poplist} {output.pcaplot}
"""
        
    
