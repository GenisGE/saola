
SNPEFF="java -jar /home/genis/software/snpEff/snpEff.jar"
BCFTOOLS="/home/genis/software/bcftools-1.15/bcftools"
VCFANNO="/home/genis/software/vcfanno_linux64"
BGZIP="bgzip"
TABIX="tabix"

PYTHON="python3"
MEANCON="/home/genis/impala/analyses/goatMapV3/load_called_genos/scripts/genload_as_meancon.py"
VEP="/home/genis/software/ensembl-vep/vep"


OUTMAIN=config["outmain"]
OUTBIG=config["outbig"]
CHROMFILE=config["chromfile"]

BCF=config["bcf"]
BCF_PREFIX = os.path.basename(BCF).replace(".bcf.gz", "")

with open(CHROMFILE, "r") as fh:
    CHROMS=[x.rstrip() for x in fh.readlines()]

rule all:
    input:
        os.path.join(OUTBIG, "annotation", BCF_PREFIX + "_VariantEffectPrediction.vep.txt"),
        os.path.join(OUTBIG, "annotation", BCF_PREFIX + "_VariantEffectPrediction.bed.gz.tbi"),
        os.path.join(OUTBIG, "bcf", BCF_PREFIX + "_vep.bcf.gz.csi"),


rule bcf_to_vcf:
    input:
        bcf=config["bcf"]
    output:
        vcf = temp(os.path.join(OUTMAIN, "bcf", "temp.vcf.gz"))
    shell: """
    {BCFTOOLS} view -Oz -o {output.vcf} {input.bcf}
"""
        


rule link_bcf:
    input:
        bcf = os.path.join(OUTBIG, "bcf", BCF_PREFIX + ".bcf.gz")
    output:
        bcf = os.path.join(OUTMAIN, "bcf", BCF_PREFIX + ".bcf.gz")
    shell: """
1;5202;0c    ln -s {input.bcf} {output.bcf}
"""

rule index_bcf:
    input:
        bcf = "{path}.bcf.gz"
    output:
        csi = "{path}.bcf.gz.csi"
    threads: 3
    shell:"""
    {BCFTOOLS} index --threads {threads} {input.bcf}
"""

rule variant_effect_vep_chr:
    input:
        bcf = os.path.join(OUTBIG, "bcf", BCF_PREFIX + ".bcf.gz"),
        gff = config["gff"],
        ref = config["ref"]
    output:
        vep = os.path.join(OUTBIG, "annotation", BCF_PREFIX + "_VariantEffectPrediction_{chrom}.vep.txt")
    threads: 20
    params:
        buffer_size = 500000
    shell: """
    {BCFTOOLS} view -r {wildcards.chrom} -m2 -M2 -Ov --threads {threads} {input.bcf} | {VEP} --format vcf --fork {threads} --buffer_size {params.buffer_size} --gff {input.gff} --fasta {input.ref}  --output_file {output.vep}
"""

    
rule combine_vep:
    input:
        expand(os.path.join(OUTBIG, "annotation", BCF_PREFIX + "_VariantEffectPrediction_{chrom}.vep.txt"), chrom=CHROMS),
    output:
        vep = os.path.join(OUTBIG, "annotation", BCF_PREFIX + "_VariantEffectPrediction.vep.txt")
    shell: """
    cat {input} > {output}
"""

    
rule vep_annotation_to_bed:
    input:
        vep = os.path.join(OUTBIG, "annotation", BCF_PREFIX + "_VariantEffectPrediction.vep.txt")
    output:
        bed = os.path.join(OUTBIG, "annotation", BCF_PREFIX + "_VariantEffectPrediction.bed.gz")
    shell: """
    grep -v "^#" {input.vep} | awk '{{split($2,pos,":"); print pos[1]"\t"pos[2]-1"\t"pos[2]"\t"$0}}' | {BGZIP} -c > {output.bed}
"""

rule tabix_annotation_bed:
    input:
        bed = os.path.join(OUTBIG, "annotation", BCF_PREFIX + "_VariantEffectPrediction.bed.gz")
    output:
        tbi = os.path.join(OUTBIG, "annotation", BCF_PREFIX + "_VariantEffectPrediction.bed.gz.tbi")
    shell: """
    {TABIX} -p bed {input.bed}
"""


rule bcf_to_vcf:
    input:
        bcf =  os.path.join(OUTBIG, "bcf", BCF_PREFIX + ".bcf.gz"),
    output:
        vcf = temp(os.path.join(OUTMAIN, "bcf", "temp.vcf.gz"))
    shell: """
    {BCFTOOLS} view -Oz -o {output.vcf} {input.bcf}
"""

    
rule add_annotation_bcf:
    input:
        vcf = os.path.join(OUTMAIN, "bcf", "temp.vcf.gz")
    output:
        bcf = os.path.join(OUTBIG, "bcf", BCF_PREFIX + "_vep.bcf.gz")
    threads: 20
    params:
        vcfanno_config = config["vcfanno_config_annotation"]
    shell: """
    {VCFANNO} -p {threads} {params.vcfanno_config} {input.vcf} | {BCFTOOLS} view -Ob -o {output.bcf}
"""

    