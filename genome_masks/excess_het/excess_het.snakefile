"""
snakemake to do excess of heterozygosity filter
configfile needs:
       - info: file with bamfile in column1 and population assignment in column2 path to bamfile. Population is just for visualization, if unknown can be set to arbitrary value for all
       - chroms: path to file with list of chromosomes to use
       - ref: path to reference
       - outmain: path to main output
       - params:
                - angsd: minQ: value
                         minMapQ: value
                - pcangsd: e: number principal components to use, if dont' know set to 0 and pcangsd infers it
                - filter: w: size in kb of window to remoe around "bad" site
                          minF: F below which significantly deviating loci are removed
                          lrt: likelihood ratio test statistic threshold to assign signficiance (e.g. lrt 24 equivales pvalue of 1e-6
      - scriptsdir: path to directory with requried scripts
"""

ANGSD="/home/genis/software/angsd/angsd"
PCANGSD="/home/genis/software/pcangsd-v.0.99/pcangsd.py" # MUST BE BEFORE PCANGSD V1, else -sites_save just prints 0 and 1s and nothing works anymore. maybe asks jonas to change
PYTHON="python3"
R="Rscript"
BEDTOOLS="/home/genis/software/bedtools2/bin/bedtools"
SAMTOOLS="/home/genis/software/samtools-1.9/samtools"

SCRIPTSDIR = config["scriptsdir"]
PLOTPCA = os.path.join(SCRIPTSDIR,"plotPCA.R")
SELECTSITES = os.path.join(SCRIPTSDIR,"select_bad_sites.py")
SUMMARIZE_PLOT = os.path.join(SCRIPTSDIR, "summarizePlotExcessHetFilter.R")

OUTMAIN=config["outmain"]


CHROMLIST=config["chroms"]
with open(CHROMLIST, "r") as fh:
    CHROMS=[x.rstrip() for x in fh.readlines()]

REF = config["ref"]

# LOAD PARAMS
MINQ=config["params"]["angsd"]["minQ"]
MINMAPQ=config["params"]["angsd"]["minMapQ"]
E=config["params"]["pcangsd"]["e"]
MINF=config["params"]["filter"]["minF"]
W=config["params"]["filter"]["w"]
LRT=config["params"]["filter"]["lrt"]


rule all:
    input:
        os.path.join(OUTMAIN, "summary", "excessHet_summary.tsv"),
        os.path.join(OUTMAIN, "summary", "excessHet_plotsBig.pdf"),
        os.path.join(OUTMAIN, "beds", "good.bed"),
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


rule do_fai:
    """do fai file from ref in case there is not"""
    input:
        ref = config["ref"]
    output:
        fai = "{}.fai".format(config["ref"])
    shell:"""
    {SAMTOOLS} faidx {input.ref}
"""


rule do_bedall:
    """do bed file with whole sequence for each included chromosome"""
    input:
        fai = ancient("{}.fai".format(REF)),
        chromlist = CHROMLIST
    output:
        bed = os.path.join(OUTMAIN, "beds", "all.bed")
    shell: """
    cat {input.chromlist} | xargs -n1 -I {{}} grep -P "^{{}}\t" {input.fai} | awk '{{print $1"\t0\t"$2}}' > {output.bed}
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
        r = "{chrom}"
    shell:"""
    {ANGSD} -GL 2 -out {params.outprefix} -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam {input.bamlist} -minmapQ {MINMAPQ} -minQ {MINQ} -r {params.r} -minmaf 0.05 -rmTrans 1 -minInd 10
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



rule do_pcangsd_stuff:
    """run pcangsd to do hwe equilitbirum test"""
    input:
        beagle = os.path.join(OUTMAIN, "beagle", "all.beagle.gz")
    output:
        multiext(os.path.join(OUTMAIN, "pcangsd", "out"),".cov", ".sites", ".lrt.sites.npy", ".inbreed.sites.npy"),
    params:
        outprefix = os.path.join(OUTMAIN, "pcangsd", "out"),
        e = E
    log: os.path.join(OUTMAIN, "pcangsd", "out.log")
    threads: 10
    shell:"""
    if [ {params.e} = 0 ]
    then
        {PYTHON} {PCANGSD} -beagle {input.beagle} -sites_save -inbreedSites -o {params.outprefix} -threads {threads} > {log}
    else
        {PYTHON} {PCANGSD} -beagle {input.beagle} -sites_save -inbreedSites -o {params.outprefix} -threads {threads} -e {params.e} > {log}
    fi
"""



rule plot_pca:
    """plots a pca to check pcangsd run seemed fine. now only does a single png pc1 vs pc2,
    in future maybe do pdf with multiple pcs? maybe can take e value and plot up to that pc?"""
    input:
        cov = os.path.join(OUTMAIN, "pcangsd", "out.cov"),
        poplist = os.path.join(OUTMAIN, "info", "pop.list")
    output:
        pcaplot = os.path.join(OUTMAIN, "pcangsd", "pca.png")
    shell:"""
    {R} {PLOTPCA} {input.cov} {input.poplist} {output.pcaplot}
"""



rule select_bad_sites:
    """from pacangsd output select 'bad' sites (sites with significant excess of heterozygosity) and save them to bed format"""
    input:
        pcangsd = multiext(os.path.join(OUTMAIN, "pcangsd", "out"), ".sites", ".lrt.sites.npy", ".inbreed.sites.npy"),
    output:
        badsitesbed = os.path.join(OUTMAIN, "beds", "badsites.bed")
    params:
        inprefix = os.path.join(OUTMAIN, "pcangsd", "out"),
        minF = MINF,
        lrt = LRT
    shell:"""
    {PYTHON} {SELECTSITES} {params.inprefix} {output.badsitesbed} {params.minF} {params.lrt}
"""


    
rule do_genome_file:
    """do genome file with length of each included chromosome/scaffold"""
    input:
        fai = ancient("{}.fai".format(REF)),
        chromlist = CHROMLIST
    output:
        genomefile = os.path.join(OUTMAIN, "info", "genome_sizes.file")
    shell: """
    cat {input.chromlist} | xargs -n1 -I {{}} grep -P "^{{}}\t" {input.fai} | awk '{{print $1"\t"$2}}' > {output.genomefile}
"""


rule make_bad_bed:
    """expands bad sites to get regions of size W around each bad sites, and then merge overlapping sites"""
    input:
        badsitesbed = os.path.join(OUTMAIN, "beds", "badsites.bed"),
        genome = os.path.join(OUTMAIN, "info", "genome_sizes.file")
    output:
        temp_bed = temp(os.path.join(OUTMAIN, "beds", "temp.bed")),
        badbed = os.path.join(OUTMAIN, "beds", "bad.bed")
    params:
        b = int(W)/2,
    shell:"""
    {BEDTOOLS} slop -b {params.b} -i {input.badsitesbed} -g {input.genome} > {output.temp_bed}
    {BEDTOOLS} merge -i {output.temp_bed} > {output.badbed}
"""


rule do_good_bed:
    """produces bed with regions to keep after applying excess of heterozyosity"""
    input:
        badbed = os.path.join(OUTMAIN, "beds", "bad.bed"),
        allbed = os.path.join(OUTMAIN, "beds", "all.bed")
    output:
        goodbed = os.path.join(OUTMAIN, "beds", "good.bed")
    shell:"""
    {BEDTOOLS} subtract -a {input.allbed} -b {input.badbed} > {output.goodbed}
"""


rule plot_summarize_filter:
    input:
        multiext(os.path.join(OUTMAIN, "pcangsd", "out"), ".cov", ".sites", ".lrt.sites.npy", ".inbreed.sites.npy"),
        badbed = os.path.join(OUTMAIN, "beds", "bad.bed"),
        allbed = os.path.join(OUTMAIN, "beds", "all.bed"),
    output:
        summary_table = os.path.join(OUTMAIN, "summary", "excessHet_summary.tsv"),
        bigpdf = os.path.join(OUTMAIN, "summary", "excessHet_plotsBig.pdf"),
    params:
        inprefix = os.path.join(OUTMAIN, "pcangsd", "out"),
        outprefix = os.path.join(OUTMAIN, "summary", "excessHet"),
        minF = MINF,
        minLRT = LRT,
    shell:"""
    {R} {SUMMARIZE_PLOT} {params.inprefix} {params.outprefix} {input.badbed} {input.allbed} {params.minF} {params.minLRT}
"""
        
    
