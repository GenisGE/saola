###################################################################################################
####################                     HELPER FUNCTIONS                      ####################
###################################################################################################


def get_ingroup_files(wildcards):
    """Helper to get ingroup-related files for merging"""
    vcf = rules.merge_prepare_ingroup.output.vcf
    index = str(vcf) + ".csi"
    return {"ingroup_vcf": vcf, "ingroup_index": index}


def get_outgroup_files(wildcards):
    """Helper to get outgroup-related files for merging"""
    outgroups = wildcards.outgroups.split("-")
    vcfs = expand(rules.merge_prepare_outgroup.output.vcf, outgroup=outgroups)
    indexes = [f"{vcf}.csi" for vcf in vcfs]
    return {"outgroup_vcfs": vcfs, "outgroup_indexes": indexes}


###################################################################################################
####################                         RULES                             ####################
###################################################################################################


rule merge_prepare_ingroup:
    """
    Prepare ingroup data for merging outgroups by splitting multiallelics and normalising indels
    """
    input:
        vcf=VCF_DIR / (REF_NAME + f"_{INGROUP}_variable_sites.bcf.gz"),
    output:
        vcf=temp(
            RESULTS_DIR
            / "tmp"
            / "merge"
            / (REF_NAME + f"_{INGROUP}_variable_sites_normsplit.bcf.gz")
        ),
    params:
        ref=config["ref"]["fa"],
    log:
        RESULTS_DIR / "tmp" / "merge" / (REF_NAME + f"_{INGROUP}_variable_sites_normsplit.log"),
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} norm \
                --atomize \
                --atom-overlaps . \
                --multiallelics - \
                --fasta-ref {params.ref} \
                --output-type b \
                --threads {threads} \
                {input.vcf} \
                > {output.vcf} \
        ) 2> {log}
        """


rule merge_prepare_ids:
    """ Prepare IDs to add to outgroups from the ingroup VCF """
    input:
        vcf=rules.merge_prepare_ingroup.output.vcf,
    output:
        ids=temp(RESULTS_DIR / "tmp" / "merge" / (REF_NAME + f"_{INGROUP}.ids")),
    params:
        query='"%CHROM:%POS:%REF:%ALT\\n%CHROM:%POS:%REF:.\\n"',
    log:
        RESULTS_DIR / "tmp" / "merge" / (REF_NAME + f"_{INGROUP}.ids.log"),
    threads: 1
    shell:
        """
        ( \
            {BCFTOOLS} query \
                --format {params.query} \
                {input.vcf} \
                > {output.ids} \
        ) 2> {log}
        """


rule merge_prepare_outgroup:
    """
    Prepare ingroup data for merging outgroups by splitting multiallelics and normalising indels,
    as well as filtering on IDs from ingroup
    """
    input:
        vcf=VCF_DIR / (REF_NAME + "_{outgroup}_sites.bcf.gz"),
        ids=rules.merge_prepare_ids.output.ids,
    output:
        vcf=temp(
            RESULTS_DIR / "tmp" / "merge" / (REF_NAME + "_{outgroup}_sites_normsplitfilter.bcf.gz")
        ),
    params:
        id='"%CHROM:%POS:%REF:%ALT"',
        ref=config["ref"]["fa"],
    log:
        RESULTS_DIR / "tmp" / "merge" / (REF_NAME + "_{outgroup}_sites_normsplitfilter.log"),
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} norm \
                --atomize \
                --atom-overlaps . \
                --multiallelics - \
                --fasta-ref {params.ref} \
                --output-type u \
                --threads {threads} \
                {input.vcf} | \
            {BCFTOOLS} annotate \
                --set-id {params.id} \
                --output-type u \
                --threads {threads} | \
            {BCFTOOLS} view \
                --include "ID=@{input.ids}" \
                --output-type b \
                --threads {threads} \
            > {output.vcf} \
        ) 2> {log}
        """


rule merge:
    """
    Merge outgroups into VCF, gather multiallelics,
    set any remaining half-genotypes missing, and recompute AN,AC,MAF tags
    """
    input:
        unpack(get_ingroup_files),
        unpack(get_outgroup_files),
    output:
        vcf=VCF_DIR / (REF_NAME + f"_{INGROUP}_variable_sites_merge{{outgroups}}.bcf.gz"),
    params:
        ref=config["ref"]["fa"],
    log:
        VCF_DIR / (REF_NAME + f"_{INGROUP}_variable_sites_merge{{outgroups}}.log"),
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} merge \
                --output-type u \
                --threads {threads} \
                {input.ingroup_vcf} \
                {input.outgroup_vcfs} | \
            {BCFTOOLS} norm \
                --fasta-ref {params.ref} \
                --multiallelics + \
                --output-type u \
                --threads {threads} | \
            {BCFTOOLS} +setGT \
                --output-type u \
                --threads {threads} \
                -- \
                --target-gt ./x \
                --new-gt . | \
            {BCFTOOLS} +fill-tags \
                --output-type b \
                --threads {threads} \
                -- \
                --tags AN,AC,MAF \
            > {output.vcf} \
        ) 2> {log}
        """
