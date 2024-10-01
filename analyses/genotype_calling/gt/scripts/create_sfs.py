#!/usr/bin/env python3

import gzip
import numpy as np
import sys


def max_an_from_samplesfile(path):
    with open(path, "r") as f:
        return 2 * len(f.read().splitlines())


def create_1d_sfs(path, max_an):
    sfs = np.zeros(max_an + 1)
    missing = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            ac, an = [int(x) for x in line.split(",")]
            if an == max_an:
                sfs[ac] += 1
            else:
                missing += 1
    print(f"{missing} sites had missing data", file=sys.stderr)
    return sfs


def create_2d_sfs(paths, max_ans):
    sfs = np.zeros([x + 1 for x in max_ans] )
    missing = 0
    with gzip.open(paths[0], "rt") as f1, gzip.open(paths[1], "rt") as f2:
        for line1, line2 in zip(f1, f2):
            ac1, an1 = [int(x) for x in line1.split(",")]
            ac2, an2 = [int(x) for x in line2.split(",")]
            if an1 == max_ans[0] and an2 == max_ans[1]:
                sfs[ac1, ac2] += 1
            else:
                missing += 1
    print(f"{missing} sites had missing data", file=sys.stderr)
    return sfs


def create_3d_sfs(paths, max_ans):
    sfs = np.zeros([x + 1 for x in max_ans] )
    missing = 0
    with gzip.open(paths[0], "rt") as f1, \
         gzip.open(paths[1], "rt") as f2, \
         gzip.open(paths[2], "rt") as f3:
        for line1, line2, line3 in zip(f1, f2, f3):
            ac1, an1 = [int(x) for x in line1.split(",")]
            ac2, an2 = [int(x) for x in line2.split(",")]
            ac3, an3 = [int(x) for x in line3.split(",")]
            if an1 == max_ans[0] and an2 == max_ans[1] and an3 == max_ans[2]:
                sfs[ac1, ac2, ac3] += 1
            else:
                missing += 1
    print(f"{missing} sites had missing data", file=sys.stderr)
    return sfs


def create_4d_sfs(paths, max_ans):
    sfs = np.zeros([x + 1 for x in max_ans] )
    missing = 0
    with gzip.open(paths[0], "rt") as f1, \
         gzip.open(paths[1], "rt") as f2, \
         gzip.open(paths[2], "rt") as f3, \
         gzip.open(paths[3], "rt") as f4:
        for line1, line2, line3, line4 in zip(f1, f2, f3, f4):
            ac1, an1 = [int(x) for x in line1.split(",")]
            ac2, an2 = [int(x) for x in line2.split(",")]
            ac3, an3 = [int(x) for x in line3.split(",")]
            ac4, an4 = [int(x) for x in line4.split(",")]
            if an1 == max_ans[0] and an2 == max_ans[1] and an3 == max_ans[2] and an4 == max_ans[3]:
                sfs[ac1, ac2, ac3, ac4] += 1
            else:
                missing += 1
    print(f"{missing} sites had missing data", file=sys.stderr)
    return sfs


def write_sfs(sfs, file):
    fmt_dim = "/".join([str(x) for x in sfs.shape])
    header = f"# Alleles: {fmt_dim}"
    fmt_sfs = [f"{x:.0f}" for x in sfs.flatten(order="C")]
    flat_sfs = ' '.join(fmt_sfs)
    print(f"{header}\n{flat_sfs}", file=file)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        # 1D SFS; first argument is samples file, second is counts file
        max_an = max_an_from_samplesfile(sys.argv[1])
        path = sys.argv[2]
        sfs = create_1d_sfs(path, max_an)
    elif len(sys.argv) == 5:
        # 2D SFS; first two arguments are samples files, last two are counts files in same order
        max_ans = [max_an_from_samplesfile(arg) for arg in sys.argv[1:3]]
        paths = sys.argv[3:5]
        sfs = create_2d_sfs(paths, max_ans)
    elif len(sys.argv) == 7:
        # 2D SFS; first three arguments are samples files, last three are counts files in same order
        max_ans = [max_an_from_samplesfile(arg) for arg in sys.argv[1:4]]
        paths = sys.argv[4:7]
        sfs = create_3d_sfs(paths, max_ans)
    elif len(sys.argv) == 9:
        # 2D SFS; first four arguments are samples files, last four are counts files in same order
        max_ans = [max_an_from_samplesfile(arg) for arg in sys.argv[1:5]]
        paths = sys.argv[5:9]
        sfs = create_4d_sfs(paths, max_ans)
    else:
        print(f"Wrong arguments", file=sys.stderr)
        sys.exit(1)

    write_sfs(sfs, sys.stdout)
