import numpy as np
import sys


def read_sfs(path):
    with open(path) as f:
        header = next(f).strip()
        shape = [int(x) for x in header.removeprefix("# Alleles: ").split("/")]
        raw = next(f).strip()
        sfs = np.array([int(x) for x in raw.split(" ")], dtype="float64").reshape(shape)
    return sfs


def fold_sfs(sfs):
    total_samples = np.sum(sfs.shape) - len(sfs.shape)
    total_per_entry = np.sum(np.indices(sfs.shape), axis=0)

    where_folded_out = total_per_entry > int(total_samples / 2)

    reversed = reverse_array(np.where(where_folded_out, sfs, 0))
    folded = np.ma.masked_array(sfs + reversed)
    folded.data[where_folded_out] = 0

    where_ambiguous = (total_per_entry == total_samples / 2.0)
    ambiguous = np.where(where_ambiguous, sfs, 0)
    folded += -0.5 * ambiguous + 0.5 * reverse_array(ambiguous)

    return folded


def reverse_array(arr):
    reverse_slice = tuple(slice(None, None, -1) for _ in arr.shape)
    return arr[reverse_slice]


def write_folded_sfs(sfs, file):
    fmt_dim = "/".join([str(x) for x in sfs.shape])
    header = f"# Alleles: {fmt_dim} (folded)"
    fmt_sfs = [f"{x:.1f}" for x in sfs.flatten(order="C")]
    flat_sfs = ' '.join(fmt_sfs)
    print(f"{header}\n{flat_sfs}", file=file)


if __name__ == "__main__":
    inpath = sys.argv[1]
    unfolded = read_sfs(inpath)
    folded = fold_sfs(unfolded)
    write_folded_sfs(folded, sys.stdout)

