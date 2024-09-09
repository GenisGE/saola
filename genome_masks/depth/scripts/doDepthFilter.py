# adpated from frederik ffs script. this one creates a bed list with sites to remove!!
# python3 doDepthFilter.py indir medianfile chromfile outprefix
#      indir: folder with input .pos.gz files
#      medianfile: file with a single number which is median depth
#      chromfile: file with list of chromosomes depth has been stimated for
#      outprefix: prefix of outpuf file. will create outprefix.bed and outprefix.log

import gzip
import sys
import glob
import os

def is_within_thresholds(x, min_thres, max_thres):
    global nout
    global nin
    
    if x == 'NA\n':
        nout += 1
        return False

    x = int(x)
    if x < min_thres or x > max_thres:
        nout += 1
        return False
    
    nin += 1
    return True


def get_bad_ranges(infile, out, min_thres, max_thres):

    #out = open(outfile, 'w+')
    is_out = False
    with gzip.open(infile, 'rt') as f:
        first_line = True
        for line in f:
            if first_line:
                first_line = False
                continue
            chrr, pos, d = line.split('\t')
            if not is_out:
                if is_within_thresholds(d, min_thres, max_thres):
                    continue
                elif not is_within_thresholds(d, min_thres, max_thres):
                    start = int(pos) - 1
                    is_out = True
            elif is_out:
                if is_within_thresholds(d, min_thres, max_thres):
                    end = int(pos) - 1
                    out.write(f'{chrr}\t{start}\t{end}\n')
                    is_out = False

                elif not is_within_thresholds(d, min_thres, max_thres):
                    continue

        if is_out:
             end = int(pos)
             out.write(f'{chrr}\t{start}\t{end}\n')
             
    #out.close()



indir = sys.argv[1] # sys.argv[0] is the name of the script!   
medianfile = sys.argv[2]
chromfile = sys.argv[3]
outpre = sys.argv[4]

outfile = outpre + ".bed"
logfile = outpre + ".log"

out = open(outfile, 'w+')
log = open(logfile, 'w+')

with open(medianfile, 'r') as fh:
    median = int(fh.read().rstrip())

# hardcoded thersholds. might want to change that so it is an argument to the script
minl = median * 0.5
maxl = median * 1.5

log.write(f"# Will create bed list with regions to remove {outfile}.\n")
log.write(f"# Median total depth per site is {median}, will only keep sites with depth between {minl} and {maxl}.\n")
log.write("#Summary\n\nChromosome\tKept\tRemoved\tTotal\n")

with open(chromfile, 'rt') as f:
    for line in f:
        
        nout = 0
        nin = 0
         
        chrr = line.strip('\n')
        infile = glob.glob(os.path.join(indir, chrr + "_*.pos.gz"))[0]

        get_bad_ranges(infile, out, minl, maxl)
 

        log.write(f"{chrr}\t{nin}\t{nout}\t{nin+nout}\n")


out.close()
log.close()

