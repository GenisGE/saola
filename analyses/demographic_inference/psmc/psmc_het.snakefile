# snakemake to run psmc. can do it using different filtering in the form of bedfiles as wildcards, and differnt minimum allele support to call hets.
# config needs:
#    samples: dic with pairs of sample:pathtobamfile
#    depths: dic with depth thresholds specific to each sample sample:[mindepth, maxdepth]
#    beds: dic with filtername:path_to_bed_with_filters
#    allele_support: list of allele support values to call heterozygotes
#    plot_params: dic with g:generation_time mu:mutation_rate
#    outmain: path to main output directory

import os

# https://github.com/lh3/psmc
VCFUTILS="/home/genis/software/bcftools/misc/vcfutils.pl"
SAMTOOLS="/home/genis/software/samtools-1.9/samtools"
BCFTOOLS="/home/genis/software/bcftools/bcftools"

PSMC_DIR="/home/genis/software/psmc"
FQ2PSMCFA=os.path.join(PSMC_DIR,"utils/fq2psmcfa")
PSMC=os.path.join(PSMC_DIR,"psmc")
#PSMC_PLOT=os.path.join(PSMC_DIR,"utils/psmc_plot.pl")
PSMC_PLOT="scripts/plotPsmc.R"

OUTMAIN = config["outmain"]
OUTBIG=config["outbig"]
BEDS = config["beds"]
REF=config["ref"]
allele_support = config["allele_support"]


MIN_BQ=30
MIN_MQ=25

wildcard_constraints:
    t = "|".join(allele_support),
    sample = "|".join(config["samples"].keys()),
    bed = "|".join(BEDS.keys()),



    
rule all:
    input:
        expand(os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "{sample}.psmc"),
               sample=config["samples"].keys(),
               bed=BEDS.keys(),
               t=allele_support
        ),
#        expand(os.path.join(OUTMAIN, "psmc_output_figures", "{bed}_{t}.png"),
#               bed=BEDS.keys(),
#               t=allele_support
#        ),
        expand(os.path.join(OUTMAIN, "hets", "{sample}_{bed}_{t}.het"),
               sample=config["samples"].keys(),
               bed=BEDS.keys(),
               t=allele_support
        ),
        expand(os.path.join(OUTMAIN, "hets", "{sample}_{bed}_{t}_notransitions.het"),
               sample=config["samples"].keys(),
               bed=BEDS.keys(),
               t=allele_support
        ),
        expand(os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "{sample}_notransitions.psmc"),
               sample=config["samples"].keys(),
               bed=BEDS.keys(),
               t=allele_support
        ),



rule gen_bcftools_genome_wide:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        vcf = os.path.join(OUTBIG, "vcf", "{sample}.bcf.gz")
    threads: 2
    shell:
        """{BCFTOOLS} mpileup  -B -Q {MIN_BQ}  -q {MIN_MQ} --threads {threads} -O u --fasta-ref {REF} --per-sample-mF -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR {input}  | {BCFTOOLS} call -Ob -o {output} --threads {threads} -c
#{BCFTOOLS} index {output}
"""

        
        

# rules to run psmc

rule gen_fq_mindepthx:
    input:
         os.path.join(OUTBIG, "vcf", "{sample}.bcf.gz")
    output:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.fq.gz"),
    params:
        mindepth=lambda wildcards: config["depths"][wildcards.sample][0],
        maxdepth=lambda wildcards: config["depths"][wildcards.sample][1],
        B = lambda wildcards: BEDS[wildcards.bed]
    threads: 3
    shell:
        """
        {BCFTOOLS} view --threads {threads} -i 'sum(INFO/DP4)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view --threads {threads} -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t} ) || GT=="hom"' | awk -f scripts/rm_indels.awk |  {VCFUTILS} vcf2fq -d {params.mindepth} -D {params.maxdepth} | gzip > {output}
        """


rule gen_psmcfa:
    input:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.fq.gz"),
    output:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.psmcfa"),
    shell:
        """
        {FQ2PSMCFA} -q 20 {input} > {output}
        """

rule run_psmc:
    input:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.psmcfa"),
    output:
        os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "{sample}.psmc")
    shell:
        """{PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o {output} {input}"""


rule plot:
    input:
        expand(os.path.join(OUTMAIN, "psmc_output", "{{bed}}", "{{t}}", "{sample}.psmc"),
               sample=config["samples"].keys())
    output:
        png=os.path.join(OUTMAIN, "psmc_output_figures", "{bed}_{t}.png")
    params:
        dir = os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}"),
        g=config["plot_params"]["g"],
	mu=config["plot_params"]["mu"]
    shell: "Rscript {PSMC_PLOT} -d {params.dir} -g {params.g} -m {params.mu} -o {output.png}"






# rules to estimate genome_wide heterozygoisties with bcftools, using same sites as in psmc

rule get_het_mindepthx:
    input:
        os.path.join(OUTBIG, "vcf", "{sample}.bcf.gz")
    output:
        os.path.join(OUTMAIN, "bcf_stats", "{bed}", "{t}", "{sample}.bcf.stats"),
    params:
        mindepth=lambda wildcards: config["depths"][wildcards.sample][0], # 30/3
        maxdepth=lambda wildcards: config["depths"][wildcards.sample][1],  # 30*2
        B = lambda wildcards: BEDS[wildcards.bed]
    shell:
        """
        {BCFTOOLS} view -i 'sum(INFO/DP4)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t} ) || GT=="hom"' | awk -f scripts/rm_indels.awk | {BCFTOOLS} stats -s - > {output}
        """

rule estimate_het:
    input:
       os.path.join(OUTMAIN, "bcf_stats", "{bed}", "{t}", "{sample}.bcf.stats")
    output:
       os.path.join(OUTMAIN, "hets", "{sample}_{bed}_{t}.het")
    shell: "grep '^PSC' {input} | awk '{{print $6/($4+$5+$6)}}' > {output}"






           

### RULES TO RUN PSMC WIHTOUT TRANSITONS


rule gen_fq_mindepthx_notransitions:
    input:
         os.path.join(OUTBIG, "vcf", "{sample}.bcf.gz")
    output:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}_notransitions.fq.gz"),
    params:
        mindepth=lambda wildcards: config["depths"][wildcards.sample][0],
        maxdepth=lambda wildcards: config["depths"][wildcards.sample][1],
        B = lambda wildcards: BEDS[wildcards.bed]
    shell:
        """
        {BCFTOOLS} view -i 'sum(INFO/DP4)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} | {BCFTOOLS} view -e 'REF =="A" & ALT[0] == "G" | REF =="G" & ALT[0] == "A" | REF =="T" & ALT[0] == "C" | REF =="C" & ALT[0] == "T"' |  {BCFTOOLS} view -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t} ) || GT=="hom"' | awk -f scripts/rm_indels.awk |  {VCFUTILS} vcf2fq -d {params.mindepth} -D {params.maxdepth} | gzip > {output}
        """


rule gen_psmcfa_notransitions:
    input:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}_notransitions.fq.gz"),
    output:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}_notransitions.psmcfa"),
    shell:
        """
        {FQ2PSMCFA} -q 20 {input} > {output}
        """

rule run_psmc_notransitions:
    input:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}_notransitions.psmcfa"),
    output:
        os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "{sample}_notransitions.psmc")
    shell:
        """{PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o {output} {input}"""


        


# rules to estimate genome_wide heterozygoisties with bcftools, using same sites as in psmc WITHOUT TRANSITIONS

rule get_het_mindepthx_notransitions:
    input:
        os.path.join(OUTBIG, "vcf", "{sample}.bcf.gz")
    output:
        os.path.join(OUTMAIN, "bcf_stats", "{bed}", "{t}", "{sample}_notransitions.bcf.stats"),
    params:
        mindepth=lambda wildcards: config["depths"][wildcards.sample][0], # 30/3
        maxdepth=lambda wildcards: config["depths"][wildcards.sample][1],  # 30*2
        B = lambda wildcards: BEDS[wildcards.bed]
    shell:
        """
        {BCFTOOLS} view -i 'sum(INFO/DP4)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t} ) || GT=="hom"' | awk -f scripts/rm_indels.awk |  {BCFTOOLS} view -e 'REF =="A" & ALT[0] == "G" | REF =="G" & ALT[0] == "A" | REF =="T" & ALT[0] == "C" | REF =="C" & ALT[0] == "T"' | {BCFTOOLS} stats -s - > {output}
        """

rule estimate_het_notransitions:
    input:
       os.path.join(OUTMAIN, "bcf_stats", "{bed}", "{t}", "{sample}_notransitions.bcf.stats")
    output:
       os.path.join(OUTMAIN, "hets", "{sample}_{bed}_{t}_notransitions.het")
    shell: "grep '^PSC' {input} | awk '{{print $6/($4+$5+$6)}}' > {output}"

