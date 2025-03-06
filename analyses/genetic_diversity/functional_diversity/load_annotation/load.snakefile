
BCFTOOLS="/home/genis/software/bcftools-1.15/bcftools"
PYTHON3="python3"
RSCRIPT="Rscript"

SCRIPTSDIR=config["scriptsdir"]
VARIANTCOUNTER=os.path.join(SCRIPTSDIR, "count_variants_categories.py")


OUTMAIN=config["outmain"]
OUTBIG=config["outbig"]

chromfile=config["chromfile"]
samplefile = config["samplefile"]

with open(chromfile, "r") as fh:
    CHROMS = [x.rstrip() for x in fh.readlines()]

with open(samplefile, "r") as fh:
    SAMPLES = [x.rstrip() for x in fh.readlines()]

    
BCF = config["bcf"]
BCF_PRE=os.path.basename(BCF).replace(".bcf.gz", "")

wildcard_constraints:
    chrom = "|".join(CHROMS),
    sample = "|".join(SAMPLES),
    s = "|".join(config["subsets"].keys())


include: "rules/filter.smk"
include: "rules/mask.smk"
include: "rules/index.smk"


rule all:
    input:
        expand(os.path.join(OUTMAIN, "annotation_load", "variant_impact_count_{s}{w}.tsv"),
               s = config["subsets"].keys(),
               w = ["", "_nomissing"]),
        expand(os.path.join(OUTMAIN, "annotation_load", "variant_consequence_count_set3a.tsv"),
               s = config["subsets"].keys(),
               w = ["", "_nomissing"]),
        expand(os.path.join(OUTBIG, "annotation_load", "bcf", BCF_PRE + "_annotatedvariants_{s}{w}_ancestralderived_variable_{c}.bcf.gz"),
               s=config["subsets"].keys(),w = ["", "_nomissing"],
               c=["LOW", "MODERATE", "HIGH"]),
        expand(os.path.join(OUTMAIN, "sfs",  BCF_PRE + "_annotatedvariants_{s}_nomissing_ancestralderived_variable_{c}_2dsfs_NorthernCentral.png"),
               s=config["subsets"].keys(),
               c=["LOW", "MODERATE", "HIGH"])
        #expand(os.path.join(OUTMAIN, "gerp_local", BCF_PRE + "_mingerp0_notransitions_10dp_3het.{sample}_{chrom}.local_genomic_gerps.txt"), sample=SAMPLES, chrom=CHROMS), 
       # expand(os.path.join(OUTMAIN, "gerp", BCF_PRE + "_mingerp0{transitions}{dp}{het}_gerpload.tsv"),
        #       transitions=["", "_notransitions"],
        #       dp=["", "_6dp","_8dp","_10dp"],
        #       het=["","_1het", "_2het", "_3het"])


rule keep_annotated_only:
    input:
        bcf = config["bcf"]
    output:
        bcf = os.path.join(OUTBIG, "annotation_load", "bcf", BCF_PRE + "_annotatedvariants.bcf.gz")
    threads: 3
    shell: """
        {BCFTOOLS} view  --threads {threads} -i "INFO/consequence != '.'" -Ob -o {output.bcf} {input.bcf}
"""


rule subset_samples:
    input:
        bcf = os.path.join(OUTBIG, "annotation_load", "bcf", BCF_PRE + "_annotatedvariants.bcf.gz")
    output:
        bcf = os.path.join(OUTBIG, "annotation_load", "bcf", BCF_PRE + "_annotatedvariants_{s}.bcf.gz")
    params:
        samples = lambda wildcards: config["subsets"][wildcards.s]["samples"]
    threads: 3
    shell: """
    {BCFTOOLS} view --threads {threads} -s {params.samples} -Ob -o {output.bcf} {input.bcf}
"""


rule count_variant_categories:
    input:
        bcf = lambda wildcards: os.path.join(OUTBIG, "annotation_load", "bcf", BCF_PRE + "_annotatedvariants_{s}" + config["subsets"][wildcards.s]["filters"] + ".bcf.gz"),
        csi = lambda wildcards: os.path.join(OUTBIG, "annotation_load", "bcf", BCF_PRE + "_annotatedvariants_{s}" + config["subsets"][wildcards.s]["filters"] + ".bcf.gz.csi")
    output:
        tsv = os.path.join(OUTMAIN, "annotation_load", "variant_impact_count_{s}.tsv")
    shell: """
    {PYTHON3} {VARIANTCOUNTER} {input.bcf} {output.tsv}
"""
