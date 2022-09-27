.PHONEY = all

THREADS := 2 # Speed is largely limited by decompression rate
INSTANCES := 5 # Total threads = THREADS * INSTANCES


batch_1:
	# Comparing to BGISeq adapters
	python3 scripts/identify_adapters.py input/project.yaml \
                --output output.batch_1 \
		--threads $(THREADS) \
		--adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
		--adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
		| parallel -j $(INSTANCES) -u

batch_3:
	# Comparing to BGISeq adapters
	python3 scripts/identify_adapters.py input/project.batch2022.yaml \
                --output output.batch_3 \
		--threads $(THREADS) \
		--adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
		--adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
		| parallel -j $(INSTANCES) -u
