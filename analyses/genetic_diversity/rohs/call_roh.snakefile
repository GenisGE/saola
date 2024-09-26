
PLINK="/home/genis/software/plink"

OUTMAIN=config["outmain"]

R="Rscript"

PLOTSCRIPTSDIR = config["plotscriptsdir"]
PLOTROHDIST = os.path.join(PLOTSCRIPTSDIR, "plotROHcalls.R")
PLOTROHGENOME = os.path.join(PLOTSCRIPTSDIR, "plotROHgenomeROHcalls.R")
PLOTROHDISTMERGED = os.path.join(PLOTSCRIPTSDIR, "plotROHcalls_mergeROHs.R")
PLOTROHGENOMEMERGED = os.path.join(PLOTSCRIPTSDIR, "plotROHgenomeROHcalls_mergeROHs.R")
PLOTROHGENOMEDENSITY = os.path.join(PLOTSCRIPTSDIR, "plotROH_fromxiaodong.R")

with open(config["inplink"]+".fam", "r") as fh:
    SAMPLESS = [l.split()[0] for l in fh.readlines()]


with open(config["autosome_bed"], "r") as fh:
    CHROMS = [x.split()[0] for x in fh.readlines()]

    
rule all:
    input:
        expand(os.path.join(OUTMAIN, "roh_plot", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.png"), s =SAMPLESS,
               win_snps = config["roh_pars"]["win_snps"],
               win_het = config["roh_pars"]["win_het"],
               min_kb_roh = config["roh_pars"]["min_kb_roh"],
               min_kb_den = config["roh_pars"]["min_kb_den"]),

        expand(os.path.join(OUTMAIN, "roh_plot_merged", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.png"), s =SAMPLESS,
               win_snps = config["roh_pars"]["win_snps"],
               win_het = config["roh_pars"]["win_het"],
               min_kb_roh = config["roh_pars"]["min_kb_roh"],
               min_kb_den = config["roh_pars"]["min_kb_den"]),

        expand(os.path.join(OUTMAIN, "plink", "{s}.bed"), s=SAMPLESS),
        
        expand(os.path.join(OUTMAIN, "roh", "roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.{outf}"), outf= ["hom", "hom.indiv"],
               win_snps = config["roh_pars"]["win_snps"],
               win_het = config["roh_pars"]["win_het"],
               min_kb_roh = config["roh_pars"]["min_kb_roh"],
               min_kb_den = config["roh_pars"]["min_kb_den"]),
        
        expand(os.path.join(OUTMAIN, "roh_plot", "roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.png"),
               win_snps = config["roh_pars"]["win_snps"],
               win_het = config["roh_pars"]["win_het"],
               min_kb_roh = config["roh_pars"]["min_kb_roh"],
               min_kb_den = config["roh_pars"]["min_kb_den"]),
        
        expand(os.path.join(OUTMAIN, "roh_plot_merged", "roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.png"),
               win_snps = config["roh_pars"]["win_snps"],
               win_het = config["roh_pars"]["win_het"],
               min_kb_roh = config["roh_pars"]["min_kb_roh"],
               min_kb_den = config["roh_pars"]["min_kb_den"]),

        expand(os.path.join(OUTMAIN, "roh_plot_density", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}_rohdensity.png"),
               s = SAMPLESS,
               win_snps = config["roh_pars"]["win_snps"],
               win_het = config["roh_pars"]["win_het"],
               min_kb_roh = config["roh_pars"]["min_kb_roh"],
               min_kb_den = config["roh_pars"]["min_kb_den"]),

              
rule make_keep_files:
    input:
        fam = config["inplink"] + ".fam"
    output:
        expand(os.path.join(OUTMAIN, "plink", "{s}.keep"), s=SAMPLESS)
    params:
        outdir = os.path.join(OUTMAIN, "plink")
    run:
        with open(input.fam, "r") as fh:
            for l in fh.readlines():
                ids = l.split()[0:2]
                s = ids[0]
                out = "{d}/{s}.keep".format(d=params.outdir, s=s)
                with open(out, "w") as fh2:
                    fh2.write(" ".join(ids) + "\t")


rule subset_plink_indivs:
    input:
        plink = multiext(config["inplink"], ".bed", ".bim", ".fam"),
        keepid = os.path.join(OUTMAIN, "plink", "{s}.keep")
    output:
        plink = multiext(os.path.join(OUTMAIN, "plink", "{s}"), ".bed", ".bim", ".fam")
    params:
        inplink=config["inplink"],
        outplink=os.path.join(OUTMAIN, "plink", "{s}"),
        chroms = " ".join(CHROMS)
    shell: """
    {PLINK} --bfile {params.inplink} --keep {input.keepid} --geno 0 --make-bed --out {params.outplink} --chr {params.chroms} --allow-extra-chr --chr-set 29 --set-missing-var-ids @:#:\$1:\$2
"""


rule call_rohs:
    input:
        plink =  multiext(os.path.join(OUTMAIN, "plink", "{s}"), ".bed", ".bim", ".fam")
    output:
        rohs =  multiext(os.path.join(OUTMAIN, "roh", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}"), ".hom", ".hom.indiv"),
        rohsummary=temp(os.path.join(OUTMAIN, "roh", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.hom.summary"))
    log: os.path.join(OUTMAIN, "roh", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.log")
    params:
        inprefix = os.path.join(OUTMAIN, "plink", "{s}"),
        outprefix = os.path.join(OUTMAIN, "roh", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}"),
    shell: """
    {PLINK} --bfile {params.inprefix} --out {params.outprefix} \
    --homozyg 'extend' --homozyg-kb {wildcards.min_kb_roh}  \
    --homozyg-density {wildcards.min_kb_den}  --homozyg-window-snp {wildcards.win_snps} \
    --homozyg-window-het {wildcards.win_het} \
    --allow-extra-chr  --chr-set 29
"""



rule combine_rohs_indivs:
    input:
        expand(os.path.join(OUTMAIN, "roh", "{s}_roh_min{{min_kb_roh}}kb_den{{min_kb_den}}_win{{win_snps}}_het{{win_het}}.hom.indiv"), s=SAMPLESS),
    output:
        f = os.path.join(OUTMAIN, "roh", "roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.hom.indiv")
    run:
        shell("cat {i} | sed -n 1p > {o} ".format(i=input[0], o=output.f))
        for f in input:
            shell("sed 1d {i} >> {o}".format(i=f, o=output.f))
        

rule combine_rohs_hom:
    input:
        expand(os.path.join(OUTMAIN, "roh", "{s}_roh_min{{min_kb_roh}}kb_den{{min_kb_den}}_win{{win_snps}}_het{{win_het}}.hom"), s=SAMPLESS),
    output:
        f = os.path.join(OUTMAIN, "roh", "roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.hom")
    run:
        shell("cat {i} | sed -n 1p > {o} ".format(i=input[0], o=output.f))
        for f in input:
            shell("sed 1d {i} >> {o}".format(i=f, o=output.f))
        


rule plot_rohs_combined:
    """ this rule uses rscript hardoced for impala samples, change for other data"""
    input:
        roh = os.path.join(OUTMAIN, "roh", "roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.hom")
    output:
        png = os.path.join(OUTMAIN, "roh_plot", "roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.png")
    params:
        g = config["autosome_bed"],
        fam =  "{prefix}.fam".format(prefix=config["inplink"])
    shell: """
    {R} {PLOTROHDIST} {input.roh} {output.png} {params.g} {params.fam}
"""


rule plot_rohs_individual:
    input:
        roh = os.path.join(OUTMAIN, "roh", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.hom"),
        plink = multiext(os.path.join(OUTMAIN, "plink", "{s}"), ".bed", ".bim", ".fam")
    output:
        png = os.path.join(OUTMAIN, "roh_plot", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.png"),
    params:
        inplink = os.path.join(OUTMAIN, "plink", "{s}"),
        fai = config["fai"]
    shell: """
    {R} {PLOTROHGENOME} {params.inplink} {input.roh} {output.png} {wildcards.s} {params.fai}
"""


rule plot_roh_density:
    """ plots roh with snp and heterzogysosity density information. Script made by Xiaodong"""
    input:
        roh = os.path.join(OUTMAIN, "roh", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.hom"),
        plink = multiext(os.path.join(OUTMAIN, "plink", "{s}"), ".bed", ".bim", ".fam")
    output:
        png = os.path.join(OUTMAIN, "roh_plot_density", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}_rohdensity.png")
    params:
        inplink = os.path.join(OUTMAIN, "plink", "{s}"),
    shell: """
    {R} {PLOTROHGENOMEDENSITY} -p {params.inplink} -s {wildcards.s} --homfile {input.roh} -o {output.png}
"""

        
rule plot_rohs_combined_merged:
    """ this rule uses rscript hardoced for impala samples, change for other data"""
    input:
        roh = os.path.join(OUTMAIN, "roh", "roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.hom")
    output:
        png = os.path.join(OUTMAIN, "roh_plot_merged", "roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.png")
    params:
        g = config["autosome_bed"],
        fam =  "{prefix}.fam".format(prefix=config["inplink"])
    shell: """
    {R} {PLOTROHDISTMERGED} {input.roh} {output.png} {params.g} {params.fam}
"""


rule plot_rohs_individual_merged:
    input:
        roh = os.path.join(OUTMAIN, "roh", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.hom"),
        plink = multiext(os.path.join(OUTMAIN, "plink", "{s}"), ".bed", ".bim", ".fam")
    output:
        png = os.path.join(OUTMAIN, "roh_plot_merged", "{s}_roh_min{min_kb_roh}kb_den{min_kb_den}_win{win_snps}_het{win_het}.png"),
    params:
        inplink = os.path.join(OUTMAIN, "plink", "{s}"),
        fai = config["fai"]
    shell: """
    {R} {PLOTROHGENOMEMERGED} {params.inplink} {input.roh} {output.png} {wildcards.s} {params.fai}
"""
