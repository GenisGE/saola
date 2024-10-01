# a1kg genotype calling pipeline

An overview of how to use the genotype calling pipeline developed for the a1kg project.

1. [Install](#Install)
2. [Quickstart](#Quickstart)
3. [Input TSV](#Input-TSV)  
4. [Configuration](#Configuration)
5. [Targets](#Targets)
	1. [Filters](#Filters) 
	2. [Masks](#Masks) 
	3. [Merging](#Merging) 
	4. [Plink](#Plink) 
6. [SFS](#SFS)
	1. [Folding](#Folding) 
	2. [Plotting](#Plotting) 
	3. [Statistics](#Statistics) 
		1. [Fst](#Fst) 


## Install

The easiest thing is to clone the entire project and change to the `gt/` subdir

```sh
git clone git@github.com:popgenDK/QCSeq.git && cd gt
```

If you're copying the code from elsewhere, ensure that you have the entire folder `gt/` and that this is the current working directory. The `snakefile` itself is *not* sufficient, since the logic is split across multiple files in `rules/` files, some rules depend on scripts in the `scripts/` subdirectory, and so on. 

## Quickstart

Create an [input TSV](#Input-TSV) and a [configuration file](#Configuration), and you can run

```sh
snakemake path/to/target --configfile path/to/configfile -n
```

For a dry run, or 

```sh
nice snakemake path/to/target --configfile path/to/configfile -j $threads
```

to run with `$threads` threads using a default nice value. See [Targets](#Targets) for details on available targets.

## Input TSV

To point the pipeline to the input data, please construct a two-column TSV, in which the first column contains group names, and the second column contains paths to BAM files. For instance,

```
groupA	/path/to/individual1.bam
groupA	/path/to/individual2.bam
groupB	/path/to/individual3.bam
groupB	/path/to/individual4.bam
groupB	/path/to/individual5.bam
```

The group names serve two purposes. 

- First, all groups are called separately in the initial genotyping. That is, for the above, individuals 1 and 2 would be called together, and individuals 3, 4, and 5 would be called together.
- Second, one group can be defined as the "ingroup" in the `config.yaml` (see [Configuration](#Configuration)), which will define a special group for the purposes of later merging groups past genotype calling. See also [`Merging`](#Merging) below.

In the typical case, there will be a species under study (e.g. giraffe) be one group (and also the "ingroup"), and then various other groups for your outgroups (e.g. okapi, pronghorn). It is perfectly possible to run with just a single group, if no split genotyping is preferred. See also input TSVs from past runs for examples.

## Configuration

A YAML configuration file is required to run the pipeline. A template is provided under `configs/template.yaml`. Some notes on the keys:

- `ingroup`: Used for merging, must refer to a group name in the input TSV file. See [Input TSV](#Input-TSV) for more.
- `results_subdir`: Name of subdirectory under `results/` where generated files will appear.
- `ref`: Name of reference file (used for file name prefixes), as well as paths to reference and reference index.
- `sites` Path to sites file used for filtering, can be left blank for initial calling and added later for downstream sites filtering. The sites file is a BED file containing *good* sites, i.e. those that should be kept.
- `tsv`: Path to input TSV file. See [Input TSV](#Input-TSV).
- `populations`: A two-column TSV file linking sample names to populations. Only required for SFS, and can left blank. See (SFS)[#SFS] section for details.
- `tmp_dir`: Path to location for `tmp` files. Could be simply `tmp`. For large sample sizes and high depth, significant time can be saved by running on `brenda` using the `tmpSSD`, but this requires permissions. Talk to Anders for details.

*Important*: This run-specific configuration file is used in conjunction with default settings defined in the top-level `config.yaml`. Check that you agree with the cutoffs defined herein. To override these defaults, simply replicate the keys in your specific config file with the desired values. See examples in the existing YAML files.

## Targets

The pipeline uses standard Snakemake logic for determining the target from the CLI, and basic familiarity with Snakemake is recommended. The [docs](https://snakemake.readthedocs.io/en/stable/) are excellent. It may also be necessary to look at the source rules in `rules/*.smk` to find the name pattern of the targets to create.

Some general pointers on targets to get started, however. Suppose we have set `results_subdir: animal` Then all files will be created in the `results/animal` sub-directory. In fact, with a few exceptions (plink files, in particular), most files will be created in `results/animal/vcf`. The main genotype call file is of the form `{refname}_{group}.bcf.gz`, so to call the `A` group (from the first column of the [input TSV](#Input-TSV)), the corresponding target will be `results/animal/vcf/Animal_A.bcf.gz` (where `Animal` is the mapping reference name, again from the input TSV). Then, to call all groups `A`, `B`, `C`, we could use a [brace expand](https://www.gnu.org/software/bash/manual/html_node/Brace-Expansion.html) and do

```sh
nice snakemake results/animal/vcf/Animal_{A,B,C}.bcf.gz --configfile path/to/configfile -j $threads
```

For convenience, common targets are defined the main snakefile using the [pseudo-rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#targets-and-aggregation) pattern. That is, we could have also run,

```sh
nice snakemake call_all_groups --configfile path/to/configfile -j $threads
```

to get the same result.


### Filters

Once you have a variant file, various filters and masks can be applied to your variants simply by adding suffixes to the file name. For instance, having generated `results/animal/vcf/Animal_A.bcf.gz` above, we could target `results/animal/vcf/Animal_A_variable.bcf.gz` to get a BCF file with only the variable sites. A (not necessarily exhaustive) list of current filter suffixes:

- `variable`: Remove monomorphic sites.
- `sites`: Remove sites not contained in sites BED file.
- `nomultiallelics`: Remove multiallelic sites.
- `noindels`: Remove indels.
- `missing`: Remove sites with any missing genotypes.

These can be composed in any order, e.g. `results/animal/vcf/Animal_A_variable_sites_missing.bcf.gz` will contain all variable sites after using the sites BED file with no missing genotypes.

Note that to maintain flexibility, the pipeline has to generate an intermediate file for each filter (and mask) in a chain. Therefore, it is good to always start by applying the `variable` filter where possible, and it may be necessary to clean unneeded files from `results/animal/vcf` occasionally.

###  Masks

The pipeline also has the concept of a mask. Whereas filters remove sites, masks set genotypes to missing. These include (curly braces here are [wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards)):

- `{dp}dp`: Mask genotypes with less than `{dp}` depth.
- `{het}het`: Mask genotypes with less than `{het}` heterozygous support.

These compose with filters in any order you like, e.g. `results/animal/vcf/Animal_A_variable_sites_10dp_3het_missing.bcf.gz`. Note that the order of the suffixes correspond to application order, from left to right. Hence, this would remove missing sites *after* masking depth and heterozygous support.

*Important*: Please note that masking may make various information in the masked variant record outdated. All masks update the `FORMAT/AC`, `FORMAT/AN`, and `FORMAT/MAF` tags after masking, but any other tags that are computed from genotype information will not be up-to-date. The [`fill-tags` plugin](https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html) for `bcftools` may be of use.

### Merging

When calling multiple groups separately, it is often required to merge the outgroups back into the ingroup. The pipeline treats this as somewhat akin to a left join on variants. In overview, all alleles observed in the ingroup will be kept, and any novel alleles introduced by the outgroups will be set to missing.

Note that unlike masks and filters, merging must take place after filtering variable sites and applying the sites filter. The target is then `results/animal/vcf/Animal_A_variable_sites_mergeB-C.bcf.gz` to merge outgroups `B` and `C` into the ingroup. See also the `merge_all_outgroups` pseudo-rule.

### Plink

The pipeline has limited support for conversion to [plink](https://www.cog-genomics.org/plink/1.9/formats) files. The naming follows the same rules as outlined above, but go in the `results/animal/plink` subdirectory and carry the usual plink extensions. To give an example, a valid target might be `results/animal/plink/Animal_A_variable_sites_mergeA-B.bed`.

## SFS

After calling genotypes, it is possible to possible to create site frequency spectra in one through four dimensions based on the genotype calls. In order to do so, the `populations` field of the [configuration file](#Configuration) must be filled with a path to a file containing a two-column TSV file: the first column should contain sample names, and the second column should contain population assignments. Sample names should match those in the variant files, `bcftools query -l $vcf` may be helpful to see which sample names have been used in `$vcf`. For instance, the following would be a valid TSV:

```
sample1	A
sample2	A
sample3	B
sample4	B
sample5	B
sample6	C
sample7	C
```

which will assign samples 1-2 to A, samples 3-5 to population B, and samples 6-7 to population C. To avoid problems with wildcards matching, and for better output file names, avoid putting punctuation in the population names. The populations TSV need not contain all samples used for genotype calling; to ignore some samples, simple leave them out of the population TSV.

After creating the population TSV and giving its path in the config file, an SFS may be created according to the pattern `results/animal/sfs/Animal_${pops}.sfs`, where `${pops}` is a `-`-delimited string of one or more population names **in alphabetically sorted order**. For instance, `results/animal/sfs/Animal_A-B-C.sfs` refers to the three-dimensional SFS of A, B, and C.

The output SFS format has two lines: the first line contains a short header giving the dimensions of the SFS, while the second line gives the values of the SFS in flat, row-major order. In other words, the second line of the output is the same format as is used in `realSFS`.

*Important*: Currently, the SFS pipeline will not consider any sites with missing genotypes. For example, when making the joint SFS of A and B, any site that has missing data in either population A or B will simply be ignored. This will typically leave enough data to make meaningful spectra, but this should be manually verified. That is, check the count of sites in the raw output SFS itself (e.g. `tail -n 1 $sfs | awk '{ sum=0; for (i=1; i<=NF; i++) { sum+= $i } print sum }'`), or inspect the [plots](#Plotting).

### Folding

By replacing the suffix `.sfs` with `.folded.sfs`, the folded SFS will be created instead. For instance, `results/animal/sfs/Animal_A-B.folded.sfs` will create the folded SFS of A and B.

### Plotting

Folded spectra in one and two dimensions can be visualised by asking for a `.pdf` file instead of `.sfs`. That is, `results/animal/sfs/Animal_A-B.pdf` will plot the folded two-dimensional SFS of A and B.

### Statistics

Various summary statistics can be calculated from one and two dimensional spectra by using the `.stats.txt` prefix.

#### Fst

In particular, by asking for the file `results/animal/sfs/Animal.fst.csv`, all two-dimensional spectra will be created and the CSV file will contain Fst values for all pairs in long format. The estimator is the (corrected) Hudson's Fst. Plotting is also provided by using the pattern `results/animal/sfs/Animal.fst.pdf`.
