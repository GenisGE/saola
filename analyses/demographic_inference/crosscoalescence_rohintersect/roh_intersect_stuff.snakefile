import itertools as it

BEDTOOLS="/home/genis/software/bedtools2/bin/bedtools"
BCFTOOLS="/home/genis/software/bcftools-1.15/bcftools"

# https://github.com/lh3/psmc
VCFUTILS="/home/genis/software/bcftools-1.15/misc/vcfutils.pl"
SAMTOOLS="/home/genis/software/samtools-1.9/samtools"

PSMC_DIR="/home/genis/software/psmc"
FQ2PSMCFA=os.path.join(PSMC_DIR,"utils/fq2psmcfa")
PSMC=os.path.join(PSMC_DIR,"psmc")

PYTHON3="python3"
FILLFQ="scripts/fill_bcf_tags.py"

RSCRIPT="Rscript"
REDUCEBED="scripts/reduceRohBed.R"

OUTMAIN=config["outmain"]
OUTBIG=config["outbig"]

BCF=config["bcf"]
BCF_PRE=os.path.basename(BCF).replace(".bcf.gz", "")


wildcard_constraints:
    s = "|".join(config["samples"] + config["othersamples"]),
    s1 = "|".join(config["samples"] + config["othersamples"]),
    s2 = "|".join(config["samples"] + config["othersamples"])



rule all:
    input:
        psmc = expand(
            os.path.join(OUTMAIN, "psmc_output",  BCF_PRE + "_pseudodiploid_fake_{s[0]}_{s[1]}.psmc"),
            s=it.combinations(config["samples"], 2)),
        hets = expand(
            os.path.join(OUTMAIN, "hets", BCF_PRE + "_pseudodiploid_{s[0]}_{s[1]}_rohsoverlap.het.txt"),
            s=it.combinations(config["samples"], 2)),
        hets_subtracted_roh = expand( os.path.join(OUTMAIN, "hets_sample_subtractroh", BCF_PRE + "_sample{s[0]}_rohs{s[1]}subtracted{s[0]}.het.txt"),
                                      s=it.permutations(config["samples"], 2)),
        psmc_subtractedrohs = expand(os.path.join(OUTMAIN, "psmc_output",  BCF_PRE + "_sample{s[0]}_rohs{s[1]}subtracted{s[0]}_fqtag.psmc"),
                                     s=it.permutations(config["samples"], 2)),
        sample_hets = expand(os.path.join(OUTMAIN, "hets_sample", BCF_PRE + "_{s}.het.txt"),
                             s = config["samples"]),
        psmc_sample =  expand(os.path.join(OUTMAIN, "psmc_output",  BCF_PRE + "_sample{s}_commonrohin_otheroh_fqtag.psmc"), s=config["samples"]),
        sample_hets_otherohs = expand(os.path.join(OUTMAIN, "hets_sample_otheroh", BCF_PRE + "_sample{s}_commonrohin_otheroh.het.txt"), s=config["samples"])


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



rule make_roh_beds_toremove:
    input:
        rohfile = config["rohfileremove"]
    output:
        rohbed = os.path.join(OUTMAIN, "beds", "{s}_rohs_toremove.bed")
    params:
        fai = config["fai"]
    shell: """
    sed 1d {input.rohfile} | awk '$1 == "{wildcards.s}"' | awk '{{print $4"\t"$7-1"\t"$8}}' | {BEDTOOLS} sort -faidx {params.fai} > {output.rohbed}
"""



rule make_noroh_beds:
    input:
        gbed = config["genomebed"],
        rohbed = os.path.join(OUTMAIN, "beds", "{s}_rohs.bed")
    output:
        norohbed = os.path.join(OUTMAIN, "beds", "{s}_norohs.bed")
    shell: """
        {BEDTOOLS} subtract -a {input.gbed} -b {input.rohbed}  > {output.sites}
"""



rule intersect_roh_beds:
    input:
        roha = os.path.join(OUTMAIN, "beds", "{s1}_rohs.bed"),
        rohb = os.path.join(OUTMAIN, "beds", "{s2}_rohs.bed"),
    output:
        roh = os.path.join(OUTMAIN, "beds", "{s1}_{s2}_common_rohs.bed"),
    shell: """
    {BEDTOOLS} intersect -a {input.roha} -b {input.rohb} > {output.roh}
"""



rule make_pseudodiploid_roh_prebcf:
    input:
        bcf = config["bcf"],
        rohsfile = os.path.join(OUTMAIN, "beds", "{s1}_{s2}_common_rohs.bed"),
    output:
        prebcf = os.path.join(OUTBIG, "bcf", "pre_" + BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap.noheader.vcf")
    shell: """
    {BCFTOOLS} view -R {input.rohsfile} -s {wildcards.s1},{wildcards.s2} {input.bcf} | {BCFTOOLS} view -e 'F_MISSING > 0 || GT == "het"' | {BCFTOOLS} query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT\t[%GT-]\n' | sed 's|0/0|0|g;s|1/1|1|g' | sed 's|-|/|g;s|/$||g;s|1/0|0/1|g' > {output.prebcf}
"""



rule get_bcf_header:
    input:
        bcf = config["bcf"]
    output:
        h = os.path.join(OUTBIG, "bcf", "new_header_{s1}{s2}merged.txt")
    params:
        rmstring = "##FORMAT=<ID=PL|##FORMAT=<ID=DP|##FORMAT=<ID=SP|##FORMAT=<ID=AD",
        newfields = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfake{s1}{s2}"
    shell: """
    cat <({BCFTOOLS} view -h {input.bcf} | egrep -v "{params.rmstring}" | grep -v "^#CHROM") <(echo "{params.newfields}")  > {output.h}
"""



rule do_proper_bcf:
    input:
        h = os.path.join(OUTBIG, "bcf", "new_header_{s1}{s2}merged.txt"),
        prebcf = os.path.join(OUTBIG, "bcf", "pre_" + BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap.noheader.vcf"),
    output:
        bcf = os.path.join(OUTBIG, "bcf",  BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap.bcf.gz"),
        tmpvcf = temp(os.path.join(OUTBIG, "bcf", "tmp_bcf_{s1}_{s2}.vcf"))
    params:
        fai = config["fai"]
    shell: """
    cat {input.h} {input.prebcf} > {output.tmpvcf}
    {BCFTOOLS} reheader --fai {params.fai} {output.tmpvcf} | {BCFTOOLS} view -Ob -o {output.bcf}
"""



rule index_bcf:
    input:
        bcf = "{path}.bcf.gz"
    output:
        csi =  "{path}.bcf.gz.csi"
    shell: """
    {BCFTOOLS} index {input.bcf}
    """



rule add_info_fq_to_bcf:
    input:
        bcf = os.path.join(OUTBIG, "bcf", BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap.bcf.gz"),
        csi = os.path.join(OUTBIG, "bcf", BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap.bcf.gz.csi")
    output:
        bcf = os.path.join(OUTBIG, "bcf",  BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap_fqfilled.bcf.gz")
    shell: """
    {PYTHON3} {FILLFQ} {input.bcf} {output.bcf}
"""



rule gen_fq:
    input:
         bcf = os.path.join(OUTBIG, "bcf", BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap_fqfilled.bcf.gz")
    output:
        fq = os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_pseudodiploid_fake_{s1}_{s2}.fq.gz")
    threads: 3
    shell:
        """
        {BCFTOOLS} view --threads {threads} {input.bcf} | {VCFUTILS} vcf2fq -d 0 -D 10000 -Q 0 | gzip > {output.fq}
        """



rule gen_psmcfa:
    input:
        fq = os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_pseudodiploid_fake_{s1}_{s2}.fq.gz")
    output:
        psmcfa = os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_pseudodiploid_fake_{s1}_{s2}.psmcfa"),
    shell:
        """
        {FQ2PSMCFA} -q 0 -g 500 {input.fq} > {output.psmcfa}
        """

rule run_psmc:
    input:
        psmcfa = os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_pseudodiploid_fake_{s1}_{s2}.psmcfa"),
    output:
        psmc = os.path.join(OUTMAIN, "psmc_output",  BCF_PRE + "_pseudodiploid_fake_{s1}_{s2}.psmc")
    shell:
        """{PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o {output.psmc} {input.psmcfa}"""




#############################
##### do heterozygosity #####
#############################
rule do_bcf_stats:
    input:
        bcf = os.path.join(OUTBIG, "bcf", BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap_fqfilled.bcf.gz")
    output:
        stats = os.path.join(OUTMAIN, "hets", BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap.stats.txt")
    shell: """
    {BCFTOOLS} stats -s - {input.bcf} > {output.stats}
"""


rule get_het:
    input:
       stats = os.path.join(OUTMAIN, "hets", BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap.stats.txt")
    output:
       het = os.path.join(OUTMAIN, "hets", BCF_PRE + "_pseudodiploid_{s1}_{s2}_rohsoverlap.het.txt")
    shell: """
    echo fake_{wildcards.s1}_{wildcards.s2} `grep '^PSC' {input.stats} | awk '{{print $6/($4+$5+$6)}}'` > {output.het}
"""




rule do_sample_stats:
    input:
        bcf = config["bcf"]
    output:
        stats = os.path.join(OUTMAIN, "hets_sample", BCF_PRE + "_{s}.stats")
    shell: """
    {BCFTOOLS} view -s {wildcards.s} {input.bcf} | {BCFTOOLS} view -i 'FMT/DP >= 10' | {BCFTOOLS} stats -s {wildcards.s} > {output.stats} 
"""


rule get_sample_het:
    input:
        stats = os.path.join(OUTMAIN, "hets_sample", BCF_PRE + "_{s}.stats")
    output:
        het = os.path.join(OUTMAIN, "hets_sample", BCF_PRE + "_{s}.het.txt")
    shell: """
        echo {wildcards.s} `grep '^PSC' {input.stats} | awk '{{print $6/($4+$5+$6)}}'` > {output.het}
"""





###########################################################################
#### do in single sample, using subtracted rohs and matched in length #####
###########################################################################


rule subtract_rohs:
    input:
        rohs1 = os.path.join(OUTMAIN, "beds", "{s1}_rohs_toremove.bed"),
        rohs2 = os.path.join(OUTMAIN, "beds", "{s2}_rohs.bed"),
    output:
        rohs = os.path.join(OUTMAIN, "beds", "{s2}_rohs_minus_{s1}_rohs.bed"),
    shell: """
    {BEDTOOLS} subtract -a {input.rohs2} -b {input.rohs1} > {output.rohs}
"""

rule reduce_subracted_rohs:
    input:
        rohs1 = os.path.join(OUTMAIN, "beds", "{s2}_rohs_minus_{s1}_rohs.bed"),
        rohsref = os.path.join(OUTMAIN, "beds", "{s1}_{s2}_common_rohs.bed")
    output:
        rohs = os.path.join(OUTMAIN, "beds", "{s2}_rohs_minus_{s1}_rohs_reduced.bed"),
    shell: """
    {RSCRIPT} {REDUCEBED} {input.rohs1} {input.rohsref} {output.rohs}
"""



rule extract_sample_subtractedroh:
    input:
        bcf = config["bcf"],
        rohs = os.path.join(OUTMAIN, "beds", "{s2}_rohs_minus_{s1}_rohs_reduced.bed"),
    output:
        bcf =  os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_nofqtag.bcf.gz")
    shell: """
    {BCFTOOLS} view -s {wildcards.s1} -R {input.rohs} {input.bcf} | {BCFTOOLS} view -e 'F_MISSING > 0' -Ob -o {output.bcf}
"""



rule add_fq_tag_sample_subtractedroh:
    input:
        bcf =  os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_nofqtag.bcf.gz"),
        csi = os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_nofqtag.bcf.gz.csi"),
    output:
        bcf = os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_fqtag.bcf.gz"),
    shell: """
    {PYTHON3} {FILLFQ} {input.bcf} {output.bcf}
"""




rule gen_fq_sample_subtractedroh:
    input:
        bcf = os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_fqtag.bcf.gz"),
    output:
        fq = os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_fqtag.fq.gz")
    threads: 3
    shell:
        """
        {BCFTOOLS} view --threads {threads} {input.bcf} | {VCFUTILS} vcf2fq | gzip > {output.fq}
        """



rule gen_psmcfa_sample_subtractedroh:
    input:
        fq = os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_fqtag.fq.gz")
    output:
        psmcfa =  os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_fqtag.psmcfa")
    shell:
        """
        {FQ2PSMCFA} -q 0 -g 500 {input.fq} > {output.psmcfa}
        """

    
rule run_psmc_sample_subtractedroh:
    input:
        psmcfa =  os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_fqtag.psmcfa")
    output:
        psmc = os.path.join(OUTMAIN, "psmc_output",  BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_fqtag.psmc"),
    shell:
        """{PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o {output.psmc} {input.psmcfa}"""




rule do_sample_stats_subtractedroh:
    input:
        bcf = os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}_fqtag.bcf.gz"),
    output:
        stats = os.path.join(OUTMAIN, "hets_sample_subtractroh", BCF_PRE +  "_sample{s1}_rohs{s2}subtracted{s1}.stats")
    shell: """
    {BCFTOOLS} view -s {wildcards.s1} {input.bcf} | {BCFTOOLS} view -i 'FMT/DP >= 10' | {BCFTOOLS} stats -s {wildcards.s1} > {output.stats} 
"""


rule get_sample_het_subtractedroh:
    input:
        stats = os.path.join(OUTMAIN, "hets_sample_subtractroh", BCF_PRE +  "_sample{s1}_rohs{s2}subtracted{s1}.stats")
    output:
        het =  os.path.join(OUTMAIN, "hets_sample_subtractroh", BCF_PRE + "_sample{s1}_rohs{s2}subtracted{s1}.het.txt")
    shell: """
        echo {wildcards.s1}_subtractedrohs{wildcards.s2} `grep '^PSC' {input.stats} | awk '{{print $6/($4+$5+$6)}}'` > {output.het}
"""


    



    
#############################################################
#### do in single sample, using rohs from other samples #####
#############################################################
rule extract_sample_otherohs:
    input:
        bcf = config["bcf"],
        rohfile = lambda wildcards: os.path.join(OUTMAIN, "beds", config["otherroh"][wildcards.s]  + "_common_rohs.bed"),
    output:
        bcf =  os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s}_commonrohin_otheroh_nofqtag.bcf.gz")
    shell: """
    {BCFTOOLS} view -s {wildcards.s} -R {input.rohfile} {input.bcf} | {BCFTOOLS} view -e 'F_MISSING > 0' -Ob -o {output.bcf}
"""



rule add_fq_tag_sample_otherrohs:
    input:
        bcf =  os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s}_commonrohin_otheroh_nofqtag.bcf.gz"),
        csi = os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s}_commonrohin_otheroh_nofqtag.bcf.gz.csi")
    output:
        bcf =  os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s}_commonrohin_otheroh_fqtag.bcf.gz")
    shell: """
    {PYTHON3} {FILLFQ} {input.bcf} {output.bcf}
"""




rule gen_fq_sample_otherrohs:
    input:
         bcf = os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s}_commonrohin_otheroh_fqtag.bcf.gz")
    output:
        fq = os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_sample{s}_commonrohin_otheroh_fqtag.fq.gz")
    threads: 3
    shell:
        """
        {BCFTOOLS} view --threads {threads} {input.bcf} | {VCFUTILS} vcf2fq | gzip > {output.fq}
        """



rule gen_psmcfa_sample_otherrohs:
    input:
        fq =  os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_sample{s}_commonrohin_otheroh_fqtag.fq.gz")
    output:
        psmcfa =  os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_sample{s}_commonrohin_otheroh_fqtag.psmcfa")
    shell:
        """
        {FQ2PSMCFA} -q 0 -g 500 {input.fq} > {output.psmcfa}
        """

    
rule run_psmc_sample_otherrohs:
    input:
        psmcfa =  os.path.join(OUTBIG, "psmc_input", BCF_PRE + "_sample{s}_commonrohin_otheroh_fqtag.psmcfa")
    output:
        psmc = os.path.join(OUTMAIN, "psmc_output",  BCF_PRE + "_sample{s}_commonrohin_otheroh_fqtag.psmc"),
    shell:
        """{PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o {output.psmc} {input.psmcfa}"""




rule do_sample_stats_otherrohs:
    input:
        bcf = os.path.join(OUTBIG, "bcf", BCF_PRE + "_sample{s}_commonrohin_otheroh_fqtag.bcf.gz")
    output:
        stats = os.path.join(OUTMAIN, "hets_sample_otheroh", BCF_PRE + "_sample{s}_commonrohin_otheroh.stats")
    shell: """
    {BCFTOOLS} view -s {wildcards.s} {input.bcf} | {BCFTOOLS} view -i 'FMT/DP >= 10' | {BCFTOOLS} stats -s {wildcards.s} > {output.stats} 
"""


rule get_sample_het_otherrohs:
    input:
        stats = os.path.join(OUTMAIN, "hets_sample_otheroh", BCF_PRE + "_sample{s}_commonrohin_otheroh.stats")
    output:
        het =  os.path.join(OUTMAIN, "hets_sample_otheroh", BCF_PRE + "_sample{s}_commonrohin_otheroh.het.txt")
    shell: """
        echo {wildcards.s} `grep '^PSC' {input.stats} | awk '{{print $6/($4+$5+$6)}}'` > {output.het}
"""
