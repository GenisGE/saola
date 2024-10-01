

rule bcf_to_vcf:
    """Convert BCF to VCF"""
    input:
        bcf = "{path}.bcf.gz"
    output:
        vcf = "{path}.vcf.gz"
    threads: 4
    shell: """
    {BCFTOOLS} view --threads {threads} -Oz -o {output.vcf} {input.bcf}
"""

