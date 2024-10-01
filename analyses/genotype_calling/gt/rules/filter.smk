rule filter_variable:
    """ Filter to retain only variable sites. """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_variable.bcf.gz",
    log:
        "{path}_variable.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} view \
                --min-alleles 2 \
                --output-type b \
                --threads {threads} \
                {input.vcf} \
                > {output.vcf} \
        ) 2> {log}
        """


rule filter_nomultiallelics:
    """ Filter to remove all indels """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_nomultiallelics.bcf.gz",
    log:
        "{path}_nomultiallelics.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} view \
                --max-alleles 2 \
                --output-type b \
                --threads {threads} \
                {input.vcf} \
                > {output.vcf} \
        ) 2> {log}
        """


rule filter_noindels:
    """ Filter to remove all indels """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_noindels.bcf.gz",
    log:
        "{path}_noindels.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} view \
                --include 'STRLEN(REF)=1 & (STRLEN(ALT)=1 | ALT=".")' \
                --output-type b \
                --threads {threads} \
                {input.vcf} \
                > {output.vcf} \
        ) 2> {log}
        """


rule filter_missing:
    """ Filter to retain only sites with no missing genotypes """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_nomissing.bcf.gz",
    log:
        "{path}_nomissing.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} view \
                --exclude 'GT[*] = "mis"' \
                --output-type b \
                --threads {threads} \
                {input.vcf} \
            > {output.vcf} \
        ) 2> {log}
        """


rule filter_sites:
    """ Filter to retain only sites in BED file """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_sites.bcf.gz",
    params:
        bed=lambda wildcards: config["sites"],
    log:
        "{path}_sites.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} view \
                --targets-file {params.bed} \
                --output-type b \
                --threads {threads} \
                {input.vcf} \
            > {output.vcf} \
        ) 2> {log}
        """


rule filter_transitions:
    """ Filter to remove sites with transition mutations """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_notransitions.bcf.gz",
    log:
        "{path}_notransitions.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} view \
                --exclude 'REF =="A" & ALT[0] == "G" | REF =="G" & ALT[0] == "A" | REF =="T" & ALT[0] == "C" | REF =="C" & ALT[0] == "T"' \
                --output-type b \
                --threads {threads} \
                {input.vcf} \
            > {output.vcf} \
        ) 2> {log}
        """


