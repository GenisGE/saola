.PHONEY = all

THREADS := 3 # Speed is largely limited by decompression rate
INSTANCES := 4 # Total threads = THREADS * INSTANCES


all:
	python3 scripts/identify_adapters.py input/project.yaml --threads $(THREADS) | parallel -j $(INSTANCES) -u
