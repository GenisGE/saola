###################################################################################################
####################                     HELPER FUNCTIONS                      ####################
###################################################################################################


def samplename_from_bam(bam):
    """Helper to get the VCF sample name from the BAM path"""
    name = Path(bam).name
    samplename = name.split(".")[0]
    return samplename


###################################################################################################
####################                         RULES                             ####################
###################################################################################################


rule bamlist:
    """ Create BAM list for a particular grouping """
    output:
        bamlist=RESULTS_DIR / "bamlist" / (REF_NAME + "_{group}.bamlist"),
    run:
        with open(output.bamlist, "w") as f:
            for group, bam in zip(GROUP_ASSIGNMENTS, BAMS):
                if group == wildcards.group:
                    print(bam, file=f)


rule grouplist:
    """ Create samples list for a particular grouping """
    output:
        sampleslist=RESULTS_DIR / "sampleslist" / (REF_NAME + "_{group}.sampleslist"),
    run:
        with open(output.sampleslist, "w") as f:
            for group, bam in zip(GROUP_ASSIGNMENTS, BAMS):
                if group == wildcards.group:
                    samplename = samplename_from_bam(bam)
                    print(samplename, file=f)
