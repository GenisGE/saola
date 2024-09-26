

BCFTOOLS="/home/genis/software/bcftools-1.15/bcftools"
BEDTOOLS="/home/genis/software/bedtools2/bin/bedtools"
PYTHON="PYTHONPATH=/home/genis/.local/lib/python3.6/site-packages/ python3" # NEED PYTHONPATH BECAUSE IN TMUX PYTHON DOES NOT POINT TO THE CORRECT LOCAL FOLDER. THIS ADDS THE FODLER TO sys.path VARIABLE. FIX SOME DAY
PIXY="/home/genis/software/pixy/pixy/__main__.py"
TABIX="tabix"
RSCRIPT="Rscript"

PLOTPIXY="scripts/summaryPixyOut.R"
DEPLETION="scripts/getDepletionBeds.R"
MAKETABLE="scripts/makeTableOut.R"
CDSOVERLAPWIN="scripts/sumOverlapCdsWin.R"
PLOTCDSOVERLAPWIN="scripts/plotCDSOverlapThetaWins.R"

INTERVENE="/home/genis/.local/bin/intervene"

OUTMAIN=config["outmain"]
OUTBIG=config["outbig"]

BCF=config["bcf"]
BCF_PRE=os.path.basename(BCF).replace(".bcf.gz", "")

wildcard_constraints:
    win = "[0-9]*",#"|".join(config["winsize"]),
    miss = "[0-9\.]*"#"|".join(config["maxmiss"])

rule all:
    input:
        expand(os.path.join(OUTMAIN, "pixy_out", "output_maxmiss{miss}_win{win}_combinedPops_pi.txt"),
               miss=config["maxmiss"], win = config["winsize"]),
        expand(os.path.join(OUTMAIN, "depletion_intersect", "maxmiss{miss}_win{win}", "out_maxmiss{miss}_win{win}_{f}.png"),
               miss=config["maxmiss"], win = config["winsize"], f=["venn", "upset"]),
        expand(os.path.join(OUTMAIN, "pixy_out", "output_maxmiss{miss}_win{win}_pi.txt"),
               miss=config["maxmiss"], win = config["winsize"]),
        expand(os.path.join(OUTMAIN, "pixy_plots", "output_maxmiss{miss}_win{win}", "output_maxmiss{miss}_win{win}_{f}.png"),
               miss = config["maxmiss"],
               win = config["winsize"],
               f=["pidxyNsitesHist", "pidxyPiDxyHist", "manhattans", "minusLog10manhattans"]),
        expand(os.path.join(OUTMAIN, "pixy_processed_out", "output_maxmiss{miss}_win{win}_allwins_cds_overlap.bed"),
               miss = config["maxmiss"],
               win = config["winsize"],
        ),
        expand(os.path.join(OUTMAIN , "pixy_processed_out", "output_maxmiss{miss}_win{win}_cdsOverlapPiDxy.png"),
               miss = config["maxmiss"],
               win = config["winsize"],
        )


rule bcf_filtermissing_vcf:
    input:
        bcf = config["bcf"]
    output:
        vcf = os.path.join(OUTBIG, "vcf", BCF_PRE + "_maxmiss{miss}.vcf.gz")
    params:
        s = ",".join(config["keepsamples"])
    threads: 3
    shell: """
    {BCFTOOLS} view -s {params.s} --threads {threads} {input.bcf} | {BCFTOOLS} view --threads {threads} -i 'F_MISSING < {wildcards.miss}' -Oz -o {output.vcf}
"""


rule tabix_vcf:
    input:
        vcf = os.path.join(OUTBIG, "vcf", BCF_PRE + "_maxmiss{miss}.vcf.gz")
    output:
        tbi = os.path.join(OUTBIG, "vcf", BCF_PRE + "_maxmiss{miss}.vcf.gz.tbi")
    shell: """
    {TABIX} {input.vcf}
"""

    

rule run_pixy:
    input:
        vcf = os.path.join(OUTBIG, "vcf", BCF_PRE + "_maxmiss{miss}.vcf.gz"),
        tbi = os.path.join(OUTBIG, "vcf", BCF_PRE + "_maxmiss{miss}.vcf.gz.tbi"),
        popfile = config["popfile"]
    output:
        out = expand(os.path.join(OUTMAIN, "pixy_out", "output_maxmiss{{miss}}_win{{win}}_{f}.txt"),
                      f= ["pi", "dxy"])
    params:
        outdir = os.path.join(OUTMAIN, "pixy_out"),
        outprefix = "output_maxmiss{miss}_win{win}",
        chunk_size = 500000
    threads: 20
    shell: """
    {PYTHON} {PIXY} --stats pi fst dxy --vcf {input.vcf} --populations {input.popfile} --n_cores {threads} --window_size {wildcards.win} --fst_type hudson --output_folder {params.outdir} --output_prefix {params.outprefix} --chunk_size {params.chunk_size}
"""



rule do_combined_pop_file:
    input:
        popfile = config["popfile"]
    output:
        popfile = os.path.join(OUTMAIN, "info", "combined_pops.list")
    shell: """
    awk '{{print $1"\tCombined"}}' {input.popfile} > {output.popfile}
"""

    

rule run_pixy_combinedPi:
    input:
        vcf = os.path.join(OUTBIG, "vcf", BCF_PRE + "_maxmiss{miss}.vcf.gz"),
        tbi = os.path.join(OUTBIG, "vcf", BCF_PRE + "_maxmiss{miss}.vcf.gz.tbi"),
        popfile = os.path.join(OUTMAIN, "info", "combined_pops.list")
    output:
        out = os.path.join(OUTMAIN, "pixy_out", "output_maxmiss{miss}_win{win}_combinedPops_pi.txt"),
    params:
        outdir = os.path.join(OUTMAIN, "pixy_out"),
        outprefix = "output_maxmiss{miss}_win{win}_combinedPops",
        chunk_size = 500000
    threads: 20
    shell: """
    {PYTHON} {PIXY} --stats pi --vcf {input.vcf} --populations {input.popfile} --n_cores {threads} --window_size {wildcards.win} --fst_type hudson --output_folder {params.outdir} --output_prefix {params.outprefix} --chunk_size {params.chunk_size}
"""


rule plot_pixy:
    input:
        pixy = expand(os.path.join(OUTMAIN, "pixy_out", "output_maxmiss{{miss}}_win{{win}}_{f}.txt"),
                      f= ["pi", "dxy"])
    output:
        png = expand(os.path.join(OUTMAIN, "pixy_plots", "output_maxmiss{{miss}}_win{{win}}", "output_maxmiss{{miss}}_win{{win}}_{f}.png"),
                     f = ["pidxyNsitesHist", "pidxyPiDxyHist", "manhattans", "minusLog10manhattans"])
    params:
        inpre = os.path.join(OUTMAIN, "pixy_out", "output_maxmiss{miss}_win{win}"),
        outpre = os.path.join(OUTMAIN, "pixy_plots", "output_maxmiss{miss}_win{win}", "output_maxmiss{miss}_win{win}"),
        g = config["genomefile"]
    shell: """
    {RSCRIPT} {PLOTPIXY} {params.inpre} {params.outpre} {params.g}
"""



rule save_combined_table:
    input:
        pixy = expand(os.path.join(OUTMAIN, "pixy_out", "output_maxmiss{{miss}}_win{{win}}_{f}.txt"),
                    f= ["pi", "dxy"])
    output:
        tsv = os.path.join(OUTMAIN, "pixy_processed_out", "output_maxmiss{miss}_win{win}_allCombinedWindowRes.tsv"),
        bed = os.path.join(OUTMAIN, "pixy_processed_out", "output_maxmiss{miss}_win{win}_allCombinedWindows.bed"),
    params:
        inprefix = os.path.join(OUTMAIN, "pixy_out", "output_maxmiss{miss}_win{win}"),
        outprefix = os.path.join(OUTMAIN, "pixy_processed_out", "output_maxmiss{miss}_win{win}"),
        g = config["genomefile"]
    shell: """
    {RSCRIPT} {MAKETABLE} {params.inprefix} {params.outprefix} {params.g}
"""


rule keep_gff_cds_chroms_autosomeorder:
    input:
        gff = config["gff"]
    output:
        bed = os.path.join(OUTMAIN, "beds", "cds_chromsfilt_autosomesort.bed")
    params:
        gbed=config["genomebed"],
        gfile=config["genomefile"]
    shell: """
    zcat {input.gff} | awk '$3 == "CDS"' | {BEDTOOLS} intersect -a stdin -b {params.gbed} | {BEDTOOLS} merge | {BEDTOOLS} sort -faidx {params.gfile} > {output.bed}
"""


rule intersect_wins_cds:
    input:
        winbed = os.path.join(OUTMAIN, "pixy_processed_out","output_maxmiss{miss}_win{win}_allCombinedWindows.bed"),
        gff = os.path.join(OUTMAIN, "beds", "cds_chromsfilt_autosomesort.bed")
    output:
        intbed = os.path.join(OUTMAIN, "pixy_processed_out", "output_maxmiss{miss}_win{win}_allwins_allcds_overlap.bed")
    shell: """
    {BEDTOOLS} intersect -wao -a {input.winbed} -b {input.gff} > {output.intbed}
"""


rule get_win_proprotion_cds_overlap:
    input:
        intbed = os.path.join(OUTMAIN, "pixy_processed_out", "output_maxmiss{miss}_win{win}_allwins_allcds_overlap.bed")
    output:
        outbed = os.path.join(OUTMAIN, "pixy_processed_out", "output_maxmiss{miss}_win{win}_allwins_cds_overlap.bed")
    shell: """
    {RSCRIPT} {CDSOVERLAPWIN} {input.intbed} {output.outbed}
"""

    
rule plot_cds_overlap_theta_wins:
    input:
        tsv= os.path.join(OUTMAIN, "pixy_processed_out", "output_maxmiss{miss}_win{win}_allCombinedWindowRes.tsv"),
        bed = os.path.join(OUTMAIN, "pixy_processed_out", "output_maxmiss{miss}_win{win}_allwins_cds_overlap.bed"),
    output:
        png = os.path.join(OUTMAIN , "pixy_processed_out", "output_maxmiss{miss}_win{win}_cdsOverlapPiDxy.png")
    params:
        outprefix = os.path.join(OUTMAIN , "pixy_processed_out", "output_maxmiss{miss}_win{win}")
    shell: """
    {RSCRIPT} {PLOTCDSOVERLAPWIN} {input.tsv} {input.bed} {params.outprefix}
"""
    
    
rule get_depletion_beds:
    input:
        pixy = expand(os.path.join(OUTMAIN, "pixy_out", "output_maxmiss{{miss}}_win{{win}}_{f}.txt"),
                    f= ["pi", "dxy"])
    output:
        beds = expand(os.path.join(OUTMAIN, "depletion_beds", "out_maxmiss{{miss}}_win{{win}}_{f}.bed"), f=["depletedPiNorthern", "depletedPiCentral", "depletedDxy"])
    params:
        inprefix = os.path.join(OUTMAIN, "pixy_out", "output_maxmiss{miss}_win{win}"),
        outprefix = os.path.join(OUTMAIN, "depletion_beds", "out_maxmiss{miss}_win{win}"),
        g = config["genomefile"]
    shell: """
    {RSCRIPT} {DEPLETION} {params.inprefix} {params.outprefix} {params.g}
"""

    


rule intersect_do_upset:
    input:
        beds = expand(os.path.join(OUTMAIN, "depletion_beds", "out_maxmiss{{miss}}_win{{win}}_{f}.bed"), f=["depletedPiNorthern", "depletedPiCentral", "depletedDxy"])
    output:
        venn = os.path.join(OUTMAIN, "depletion_intersect", "maxmiss{miss}_win{win}", "out_maxmiss{miss}_win{win}_upset.png"),
        sets = directory(os.path.join(OUTMAIN, "depletion_intersect", "maxmiss{miss}_win{win}", "sets"))
    params:
        names="pi_depleted_Northern,pi_depleted_Central,dxy_depleted",
        outdir = os.path.join(OUTMAIN, "depletion_intersect", "maxmiss{miss}_win{win}"),
        project = "out_maxmiss{miss}_win{win}"
    shell: """
    {INTERVENE} upset -i {input.beds} --names {params.names} -o {params.outdir} --save-overlaps --project {params.project} --figtype png  --mbcolor "#000000" --sbcolor "#000000"
"""


rule intersect_do_venn:
    input:
        beds = expand(os.path.join(OUTMAIN, "depletion_beds", "out_maxmiss{{miss}}_win{{win}}_{f}.bed"), f=["depletedPiNorthern", "depletedPiCentral", "depletedDxy"])
    output:
        venn = os.path.join(OUTMAIN, "depletion_intersect", "maxmiss{miss}_win{win}", "out_maxmiss{miss}_win{win}_venn.png"),
        sets = directory(os.path.join(OUTMAIN, "depletion_intersect", "maxmiss{miss}_win{win}", "sets"))
    params:
        names="pi_depleted_Northern,pi_depleted_Central,dxy_depleted",
        outdir = os.path.join(OUTMAIN, "depletion_intersect", "maxmiss{miss}_win{win}"),
        project = "out_maxmiss{miss}_win{win}"
    shell: """
    {INTERVENE} venn -i {input.beds} --names {params.names} -o {params.outdir} --save-overlaps --project {params.project} --figtype png
"""
