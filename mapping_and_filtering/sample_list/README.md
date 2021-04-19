# Creating sample list for use with the PALEOMIX BAM pipeline

The script `sample_tsv_to_yaml.py` takes a TSV file containing a table of samples and associated FASTQ files and outputs these in YAML format. This output can append to a PALEOMIX BAM pipeline YAML file to create a runnable project:

	paleomix bam makefile > project.yaml
	./scripts/sample_tsv_to_yaml.py input/samples.tsv >> project.yaml
