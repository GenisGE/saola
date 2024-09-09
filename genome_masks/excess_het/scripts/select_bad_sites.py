import numpy as np
import pandas as pd
import sys

inprefix = sys.argv[1]
outbed = sys.argv[2]
min_f = float(sys.argv[3])
min_lrt = float(sys.argv[4])

lrt = np.load("{}.lrt.sites.npy".format(inprefix))
f = np.load("{}.inbreed.sites.npy".format(inprefix))

badmask = np.logical_and(lrt>min_lrt, f<min_f)

with open("{}.sites".format(inprefix), "r") as fh:
    sites = [x.strip() for x in fh.readlines()]


sites=np.asarray(sites)
badsites = sites[badmask]

chrom = ["_".join(x.split("_")[:-1]) for x in badsites]
start = [int(x.split("_")[-1])-1 for x in badsites]
end = [int(x.split("_")[-1]) for x in badsites]

bed = pd.DataFrame({"chr":chrom, "start":start, "end":end})
bed.to_csv(outbed, header=False, index=False, sep="\t")
