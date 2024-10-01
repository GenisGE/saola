rule depth:
    """ Count depths """
    input:
        vcf="{path}.bcf.gz",
    output:
        txt="{path}.depth.txt",
    params:
        awk="{a[$1+$2+$3+$4]++} END {for(val in a){print val, a[val]}}",
    shell:
        """
        {BCFTOOLS} query -f '%INFO/DP4\\n' {input.vcf} |
            awk -F"," '{params.awk}' |
            sort -gk1,1 \
        > {output.txt}
        """


rule plot_depth:
    """ Plot depth bar plot """
    input:
        txt=rules.depth.output.txt,
    output:
        pdf="{path}.depth.pdf",
    params:
        width=10,
        height=6,
    shell:
        """
        Rscript scripts/plot_depth.R {input.txt} {output.pdf} {params.width} {params.height}
        """
