.PHONEY = all

THREADS := 2 # Speed is largely limited by decompression rate
INSTANCES := 5 # Total threads = THREADS * INSTANCES


all:
	# Comparing to BGISeq adapters
	python3 scripts/identify_adapters.py input/project.yaml \
		--threads $(THREADS) \
		--adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
		--adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
		| parallel -j $(INSTANCES) -u
