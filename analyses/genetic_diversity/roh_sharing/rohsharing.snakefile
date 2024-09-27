# before running do 
# module load gcc/11.2.0
# module load R

BEDTOOLS="/maps/projects/alab/people/lpq293/scratch_data/software/bedtools2/bin/bedtools"
RSCRIPT="Rscript"
ROHSHARING="scripts/checkOverlappingRohs.R"
OVERLAPROHSHARCDS="scripts/overlapWithCDS.R"
PLOTROHSHARING="scripts/plotRohSharingAndCDS.R"

OUTMAIN=config["outmain"]

wildcard_constraints:
    s = "|".join(config["samples"]),
    g = "|".join(config["groups"].keys())

rule all:
    input:
        expand(os.path.join(OUTMAIN, "shuffle", "{i}", "rohsharing_tables", "rohsharing_summary_withcdsprop.tsv"),
                    i = range(1, config["n_shuffles"] + 1)),
        os.path.join(OUTMAIN, "beds", "intersects_beds.bed"),
        os.path.join(OUTMAIN, "rohsharing_tables", "rohsharing_summary.tsv"),
        os.path.join(OUTMAIN, "beds", "intersects_beds_withcdsseq.bed"),
        os.path.join(OUTMAIN, "rohsharing_tables", "rohsharing_summary_withcdsprop.tsv"),
        os.path.join(OUTMAIN, "plots", "barplots_rohSharing.png"),
        os.path.join(OUTMAIN, "plots", "barplots_rohCDSOverlap.png")



rule make_roh_beds:
    input:
        rohfile = config["rohfile"]
    output:
        rohbed = os.path.join(OUTMAIN, "beds", "{s}_rohs.bed")
    params:
        fai = config["fai"]
    shell: """
    sed 1d {input.rohfile} | awk '$1 == "{wildcards.s}"' | awk '{{print $4"\t"$7-1"\t"$8}}' | {BEDTOOLS} sort -faidx {params.fai} > {output.rohbed}
"""


rule intersect_rohs:
    input:
        rohbeds = expand(os.path.join(OUTMAIN, "beds", "{s}_rohs.bed"), s=config["samples"])
    output:
        rohintersect = os.path.join(OUTMAIN, "beds", "intersects_beds.bed")
    params:
        g = config["genomebed"],
        names = " ".join(config["samples"])
    shell: """
    {BEDTOOLS} multiinter -g {params.g} -empty -names {params.names} -header -i {input} > {output.rohintersect}
"""


rule do_rohsharing:
    input:
        rohs =  os.path.join(OUTMAIN, "beds", "intersects_beds.bed"),
        pl = multiext(config["inplink"], ".bed", ".fam", ".bim")
    output:
        outsv = os.path.join(OUTMAIN, "rohsharing_tables", "rohsharing_summary.tsv"),
    params:
        plprefix = config["inplink"],
        samples=",".join(config["samples"])
    shell: """
    {RSCRIPT} {ROHSHARING} {params.plprefix} {input.rohs} {output.outsv} {params.samples}
"""


rule keep_gff_cds_chroms_autosomeorder:
    input:
        gff = config["gff"]
    output:
        bed = os.path.join(OUTMAIN, "beds", "cds_chromsfilt_autosomesort.bed")
    params:
        gbed=config["genomebed"],
    shell: """
    zcat {input.gff} | awk '$3 == "CDS"' | {BEDTOOLS} intersect -a stdin -b {params.gbed} | {BEDTOOLS} merge | {BEDTOOLS} sort -faidx {params.gbed} > {output.bed}
"""


rule intersect_wins_cds:
    input:
        winbed = os.path.join(OUTMAIN, "beds", "intersects_beds.bed"),
        gff = os.path.join(OUTMAIN, "beds", "cds_chromsfilt_autosomesort.bed")
    output:
        intbed = os.path.join(OUTMAIN, "beds", "intersects_beds_withcdsseq.bed")
    shell: """
    {BEDTOOLS} intersect -wao -a {input.winbed} -b {input.gff} > {output.intbed}
"""


rule combine_rohsharing_cds:
    input:
        rohsharing = os.path.join(OUTMAIN, "rohsharing_tables", "rohsharing_summary.tsv"),
        cdsbed = os.path.join(OUTMAIN, "beds", "intersects_beds_withcdsseq.bed")
    output:
        rohsharingcds = os.path.join(OUTMAIN, "rohsharing_tables", "rohsharing_summary_withcdsprop.tsv"),
    shell: """
    {RSCRIPT} {OVERLAPROHSHARCDS} {input.cdsbed} {input.rohsharing} {output.rohsharingcds}
"""

rule plot_rohsharing_cds:
    input:
        rohsharingcds = os.path.join(OUTMAIN, "rohsharing_tables", "rohsharing_summary_withcdsprop.tsv")
    output:
        png1 = os.path.join(OUTMAIN, "plots", "barplots_rohSharing.png"),
        png2 = os.path.join(OUTMAIN, "plots", "barplots_rohCDSOverlap.png")
    params:
        outpre = os.path.join(OUTMAIN, "plots", "barplots")
    shell: """
    {RSCRIPT} {PLOTROHSHARING} {input.rohsharingcds} {params.outpre}
"""






### estimate confidence in roh sharing by shuffling roh bed files for each sample, then redo the whole pipeline.
### that should give us random result

rule shuffle_roh_beds:
    input:
        rohbed = os.path.join(OUTMAIN, "beds", "{s}_rohs.bed")
    output:
        rohbed = temp(os.path.join(OUTMAIN, "shuffle", "{i}", "beds", "{s}_rohs.bed"))
    params:
        fai = config["fai"]
    shell: """
    {BEDTOOLS} shuffle -i {input.rohbed} -g {params.fai} -seed {wildcards.i} | {BEDTOOLS} sort -faidx {params.fai} > {output.rohbed}
"""


rule intersect_shuffled_rohs:
    input:
        rohbeds = expand(os.path.join(OUTMAIN, "shuffle", "{{i}}", "beds", "{s}_rohs.bed"), s=config["samples"])
    output:
        rohintersect = temp(os.path.join(OUTMAIN, "shuffle", "{i}", "beds", "intersects_beds.bed"))
    params:
        g = config["genomebed"],
        names = " ".join(config["samples"])
    shell: """
    {BEDTOOLS} multiinter -g {params.g} -empty -names {params.names} -header -i {input} > {output.rohintersect}
"""


rule do_rohsharing_shuffled:
    input:
        rohs =  os.path.join(OUTMAIN, "shuffle", "{i}", "beds", "intersects_beds.bed"),
        pl = multiext(config["inplink"], ".bed", ".fam", ".bim")
    output:
        outsv = temp(os.path.join(OUTMAIN, "shuffle", "{i}", "rohsharing_tables", "rohsharing_summary.tsv")),
    params:
        plprefix = config["inplink"],
        samples=",".join(config["samples"])
    shell: """
    {RSCRIPT} {ROHSHARING} {params.plprefix} {input.rohs} {output.outsv} {params.samples}
"""



rule intersect_wins_cds_shuffled:
    input:
        winbed = os.path.join(OUTMAIN, "shuffle", "{i}", "beds", "intersects_beds.bed"),
        gff = os.path.join(OUTMAIN, "beds", "cds_chromsfilt_autosomesort.bed")
    output:
        intbed = os.path.join(OUTMAIN, "shuffle", "{i}",  "beds", "intersects_beds_withcdsseq.bed")
    shell: """
    {BEDTOOLS} intersect -wao -a {input.winbed} -b {input.gff} > {output.intbed}
"""


rule combine_rohsharing_cds_shuffled:
    input:
        rohsharing = os.path.join(OUTMAIN, "shuffle", "{i}", "rohsharing_tables", "rohsharing_summary.tsv"),
        cdsbed = os.path.join(OUTMAIN,"shuffle", "{i}",  "beds", "intersects_beds_withcdsseq.bed")
    output:
        rohsharingcds = os.path.join(OUTMAIN, "shuffle", "{i}", "rohsharing_tables", "rohsharing_summary_withcdsprop.tsv"),
    shell: """
    {RSCRIPT} {OVERLAPROHSHARCDS} {input.cdsbed} {input.rohsharing} {output.rohsharingcds}
"""

