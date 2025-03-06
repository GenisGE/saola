rule mask_low_depth:
    """ Set low-depth genotypes as missing and recompute tags """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_{dp}dp.bcf.gz",
    log:
        "{path}_{dp}dp.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} +setGT \
                --output-type u \
                --threads {threads} \
                {input.vcf} \
                -- \
                --include "FMT/DP<{wildcards.dp}" \
                --target-gt q \
                --new-gt . | \
            {BCFTOOLS} +fill-tags \
                --output-type b \
                --threads {threads} \
                -- \
                --tags AN,AC,MAF \
                > {output.vcf} \
        ) 2> {log}
        """


rule mask_het_allele_support:
    """ Set heterozygous genotypes with low support as missing and recompute tags """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_{het}het.bcf.gz",
    log:
        "{path}_{het}het.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} +setGT \
                --output-type u \
                --threads {threads} \
                {input.vcf} \
                -- \
                --include "FMT/GT=='het' & (FMT/AD[*:0]<{wildcards.het} | FMT/AD[*:1]<{wildcards.het})" \
                --target-gt q \
                --new-gt . | \
            {BCFTOOLS} +fill-tags \
                --output-type b \
                --threads {threads} \
                -- \
                --tags AN,AC,MAF \
                > {output.vcf} \
        ) 2> {log}
        """
