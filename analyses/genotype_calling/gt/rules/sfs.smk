###################################################################################################
####################                          SETUP                            ####################
###################################################################################################

SFS_SAMPLE_NAMES = list()
SFS_POPULATION_ASSIGNMENTS = list()
with open(config["populations"]) as f:
    reader = csv.reader(f, delimiter="\t", quotechar='"')
    for sample_name, population in reader:
        SFS_SAMPLE_NAMES.append(sample_name)
        SFS_POPULATION_ASSIGNMENTS.append(population)


POPULATIONS = set(SFS_POPULATION_ASSIGNMENTS)
POPULATION_CONSTRAINT = "|".join(POPULATIONS)

POPULATION_PAIRS = [f"{x}-{y}" for x, y in it.combinations(sorted(list(POPULATIONS)), 2)]
POPULATION_PAIRS_CONSTRAINT = "|".join(POPULATION_PAIRS)

wildcard_constraints:
    population=POPULATION_CONSTRAINT,
    population_pair=POPULATION_PAIRS_CONSTRAINT,
#    populations=f"{POPULATION_CONSTRAINT}|{POPULATION_PAIRS_CONSTRAINT}",


###################################################################################################
####################                     HELPER FUNCTIONS                      ####################
###################################################################################################


def get_names(wildcards):
    return wildcards.populations.split("-")


def get_counts(wildcards):
    """Helper to get counts for multiple populations"""
    populations = get_names(wildcards)
    counts = expand(rules.count_population.output.txt, population=populations)
    return counts


def get_sampleslists(wildcards):
    """Helper to get population list for multiple populations"""
    populations = get_names(wildcards)
    sampleslists = expand(rules.population_list.output.sampleslist, population=populations)
    return sampleslists


###################################################################################################
####################                         RULES                             ####################
###################################################################################################


rule symlink_vcf:
    """ Symlink main VCF into SFS temp dir to avoid polluting main results folder """
    input:
        vcf=VCF_DIR / f"{REF_NAME}_{INGROUP}.bcf.gz",
    output:
        vcf=SFS_DIR / "temp" / f"{REF_NAME}_{INGROUP}.bcf.gz",
    shell:
        """ ln -s $(pwd)/{input.vcf} {output.vcf} """


rule population_list:
    """ Create samples list for a single populations """
    output:
        sampleslist=SFS_DIR / (REF_NAME + "_{population}.sampleslist"),
    run:
        with open(output.sampleslist, "w") as f:
            for samplename, population in zip(SFS_SAMPLE_NAMES, SFS_POPULATION_ASSIGNMENTS):
                if population == wildcards.population:
                    print(samplename, file=f)


rule count_population:
    """ Get AC,AN counts from population for all sites """
    input:
        vcf=SFS_DIR
        / "temp"
        / f"{REF_NAME}_{INGROUP}_sites_nomultiallelics_noindels_10dp_2het.bcf.gz",
        sampleslist=SFS_DIR / (REF_NAME + "_{population}.sampleslist"),
    output:
        txt=SFS_DIR / f"{REF_NAME}_{{population}}.counts.txt.gz",
    log:
        SFS_DIR / f"{REF_NAME}_{{population}}.counts.log",
    threads: 4
    shell:
        """
        ( \
            {BCFTOOLS} view \
                --samples-file {input.sampleslist} \
                --output-type u \
                --threads {threads} \
                {input.vcf} | \
            {BCFTOOLS} +fill-tags \
                --output-type u \
                --threads {threads} \
                -- \
                --tags AN,AC,MAF | \
            {BCFTOOLS} query \
                --format '%INFO/AC,%INFO/AN\\n' | \
            tr '.' '0' | \
            gzip > {output.txt} \
        ) 2> {log}
        """


rule create_sfs:
    """ Create SFS in ANGSD format from AC,AN counts """
    input:
        counts=get_counts,
        sampleslists=get_sampleslists,
    output:
        sfs=SFS_DIR / f"{REF_NAME}_{{populations}}.sfs",
    shell:
        """
        {PYTHON} scripts/create_sfs.py {input.sampleslists} {input.counts} > {output.sfs}
        """


rule fold_sfs:
    """ Fold SFS using δaδi's logic """
    input:
        sfs=rules.create_sfs.output.sfs,
    output:
        sfs=SFS_DIR / f"{REF_NAME}_{{populations}}.folded.sfs",
    shell:
        """
        {PYTHON} scripts/fold_sfs.py {input.sfs} > {output.sfs}
        """


rule stats_1d_sfs:
    """ Stats for 1D SFS """
    input:
        sfs=SFS_DIR / f"{REF_NAME}_{{population}}.sfs",
    output:
        txt=SFS_DIR / f"{REF_NAME}_{{population}}.stats.txt",
    shell:
        """
        {RSCRIPT} scripts/stats_1d_sfs.R {input.sfs} > {output.txt}
        """


rule plot_1d_sfs:
    """ Plot 1D SFS """
    input:
        sfs=SFS_DIR / f"{REF_NAME}_{{population}}.sfs",
        samples=rules.population_list.output.sampleslist,
    output:
        pdf=SFS_DIR / f"{REF_NAME}_{{population}}.pdf",
    shell:
        """
        {RSCRIPT} scripts/plot_1d_sfs.R {input.sfs} {input.samples} {output.pdf}
        """


rule stats_2d_sfs:
    """ Stats for 2D SFS """
    input:
        sfs=SFS_DIR / f"{REF_NAME}_{{population_pair}}.sfs",
    output:
        txt=SFS_DIR / f"{REF_NAME}_{{population_pair}}.stats.txt",
    shell:
        """
        {RSCRIPT} scripts/stats_2d_sfs.R {input.sfs} > {output.txt}
        """


rule plot_2d_sfs:
    """ Plot 2D SFS """
    input:
        sfs=SFS_DIR / f"{REF_NAME}_{{population_pair}}.sfs",
    output:
        pdf=SFS_DIR / f"{REF_NAME}_{{population_pair}}.pdf",
    shell:
        """
        {RSCRIPT} scripts/plot_2d_sfs.R {input.sfs} {output.pdf}
        """


rule plot_fst:
    """ Plot matrix of Fst values """
    input:
        stats=expand(rules.stats_2d_sfs.output.txt, population_pair=POPULATION_PAIRS),
    output:
        pdf=SFS_DIR / f"{REF_NAME}.fst.pdf",
        csv=SFS_DIR / f"{REF_NAME}.fst.csv",
    shell:
        """
        {RSCRIPT} scripts/plot_fst.R {output.pdf} {output.csv} {input.stats}
        """


rule all_sfs:
    """ Fake rule to create all SFS output """
    input:
        expand(rules.plot_1d_sfs.output.pdf, population=POPULATIONS),
        expand(rules.stats_1d_sfs.output.txt, population=POPULATIONS),
        expand(rules.plot_2d_sfs.output.pdf, population_pair=POPULATION_PAIRS),
        expand(rules.stats_2d_sfs.output.txt, population_pair=POPULATION_PAIRS),
        expand(rules.fold_sfs.output.sfs, populations=list(POPULATIONS)+list(POPULATION_PAIRS)),
        SFS_DIR / f"{REF_NAME}.fst.pdf",
