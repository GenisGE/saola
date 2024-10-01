
rule index:
    """ Index BCF """
    input:
        vcf="{path}.bcf.gz",
    output:
        csi="{path}.bcf.gz.csi",
    threads: 4
    shell:
        """
        {BCFTOOLS} index --threads {threads} {input.vcf}
        """


rule index_vcf:
    """ Index VCF """
    input:
        vcf="{path}.vcf.gz",
    output:
        csi="{path}.vcf.gz.csi",
    threads: 4
    shell:
        """
        {BCFTOOLS} index --threads {threads} --threads {threads} {input.vcf}
        """

        
rule tabix_index_vcf:
    """ Index VCF with tabix"""
    input:
        vcf = "{path}.vcf.gz"
    output:
        tbi = "{path}.vcf.gz.tbi"
    threads: 4
    shell: """
    {BCFTOOLS} index --threads {threads} --tbi {input.vcf}
"""


rule n_sites:
    """ Get number of sites from index """
    input:
        csi="{path}.bcf.gz.csi",
    output:
        txt="{path}.n_sites.txt",
    threads: 1
    shell:
        """
        {BCFTOOLS} index --nrecords {input.csi} > {output.txt}
        """
