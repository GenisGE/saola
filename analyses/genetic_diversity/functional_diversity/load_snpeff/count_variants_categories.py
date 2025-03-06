import pysam
import numpy as np
import pandas as pd
import sys
import warnings

def update_category(entry, gt):
    """ update resutls dic for a given site and sample"""
    entry["total_alleles"] += 2
    entry["total_derived"] += gt
    entry["total_sites"] += 1
    if gt == 0:
        entry["total_hom_ref_sites"] += 1
    elif gt == 1:
        entry["total_het_sites"] += 1
    elif gt == 2:
        entry["total_hom_der_sites"] += 1
    return entry


    
def count_categories(vcf):

    res = {}
    categories = ["MODIFIER", "LOW", "MODERATE", "HIGH"]
        
    with pysam.VariantFile(vcf) as fh:
        samples = list(fh.header.samples)

        # initialize results dic with entries for each samples and for each sample for each categories
        for x in samples:
            res[x] = {}
            for y in categories:
                res[x][y] = {"total_alleles": 0,
                             "total_derived": 0,
                             "total_sites": 0,
                             "total_hom_ref_sites": 0,
                             "total_hom_der_sites": 0,
                             "total_het_sites": 0}
        
        for rec in fh.fetch():
            
            c = rec.info["ANN"][0].split("|")[2]
            
            if c not in categories:
                # for when something weird happens
                warnings.warn("Found non-standard category " + c + ", will ignore and continue.")
                continue
            
            for s in samples:

                if rec.samples[s]['GT'][0] == None or rec.samples[s]['GT'][1] == None:
                # skip this one if genotype is missing
                    continue
                
                gt = sum(rec.samples[s]["GT"])
                
                if gt > 2:
                    # in case there are multiallelics
                    warnings.warn("Found what looks like multiallelic site, will skip and continue. But there might be others undetected, don't trust the results.")
                    continue
                
                res[s][c] = update_category(res[s][c], gt)
    return res




if __name__ == "__main__":

    vcf = sys.argv[1]
    outfile = sys.argv[2]

    outdic = count_categories(vcf)

    # https://stackoverflow.com/a/54300940
    outdf = pd.concat({k: pd.DataFrame(v).T for k,v in outdic.items()}, axis=0)
    outdf.to_csv(outfile, sep="\t", index_label = ["sample", "category"])
