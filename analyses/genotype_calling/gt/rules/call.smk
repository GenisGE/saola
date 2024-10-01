###################################################################################################
####################                     HELPER FUNCTIONS                      ####################
###################################################################################################


def get_baq():
    """Get bcftools baq argument based on config"""
    baq = config["call_parameters"]["baq"]
    if baq == "with":
        baq = "--full-BAQ"
    elif baq == "without":
        baq = "--no-BAQ"
    elif baq == "partial":
        baq = ""
    else:
        print("BAQ option must be one of 'with', 'without', or 'partial'", file=sys.stderr)
        sys.exit(1)
    return baq


def get_call_vcfs(wildcards):
    """Get all call VCFs for regions with size more than cutoff"""
    vcfs = expand(str(rules.call.output.vcf), region=get_large_regions(), allow_missing=True)
    return vcfs


def get_large_regions():
    """Get list of large regions in FAI blocked into scatter/gather blocks"""
    fai = config["ref"]["fai"]
    block_size = config["call_parameters"]["block_size"]
    cutoff = config["call_parameters"]["size_cutoff"]
    regions = []
    with open(fai) as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"')
        for row in reader:
            name = row[0]
            size = int(row[1])
            if size >= cutoff:
                for start in range(1, size, block_size):
                    end = start + block_size - 1
                    end = min(end, size + 1)
                    region = f"{name}:{start}-{end}"
                    regions.append(region)
    return regions


###################################################################################################
####################                         RULES                             ####################
###################################################################################################


rule call:
    """ Call per region """
    input:
        bamlist=rules.bamlist.output.bamlist,
    output:
        vcf=temp(TMP_DIR / "split" / REF_NAME / (REF_NAME + "_{group}_{region}.bcf.gz")),
    params:
        ref=config["ref"]["fa"],
        baq=get_baq(),
        bq=config["call_parameters"]["bq"],
        mq=config["call_parameters"]["mq"],
    log:
        TMP_DIR / "split" / REF_NAME / (REF_NAME + "_{group}_{region}.log"),
    threads: 1
    shell:
        """
        ( \
            {BCFTOOLS} mpileup \
                {params.baq} \
                --annotate FORMAT/DP,FORMAT/SP,FORMAT/AD \
                --bam-list {input.bamlist} \
                --fasta-ref {params.ref} \
                --min-BQ {params.bq} \
                --min-MQ {params.mq} \
                --output-type u \
                --per-sample-mF \
                --regions '{wildcards.region}' \
                --threads {threads} \
            | {BCFTOOLS} call \
                --multiallelic-caller \
                --output-type b \
                --output '{output.vcf}' \
                --threads {threads} \
        ) 2> '{log}'
        """


rule vcflist:
    """ Put VCFs in list to avoid escapes and file number limit """
    input:
        vcfs=get_call_vcfs,
    output:
        vcflist=VCF_DIR / (REF_NAME + "_{group}.vcflist"),
    run:
        with open(output.vcflist, "w") as outfile:
            print("\n".join(input.vcfs), file=outfile)


rule concat_calls:
    """ Concatenate VCFs """
    input:
        vcfs=get_call_vcfs,
        vcflist=rules.vcflist.output.vcflist,
    output:
        vcf=protected(VCF_DIR / (REF_NAME + "_{group}.bcf.gz")),
    log:
        VCF_DIR / (REF_NAME + "_{group}.log"),
    threads: 4
    shell:
        """
        ( \
            {BCFTOOLS} concat \
                --file-list {input.vcflist} \
                --naive \
                --output-type b \
                --threads {threads} \
                > {output.vcf} \
        ) 2> {log}
        """
