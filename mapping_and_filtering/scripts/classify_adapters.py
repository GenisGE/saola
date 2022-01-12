#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import re
import sys

from pathlib import Path

KMER_LENGTH = 9

# Illumina: https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
# BGI:      https://en.mgi-tech.com/Download/download_file/id/71

ADAPTER_1 = {
    "BGI": "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
    "ILLUMINA": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
}

ADAPTER_2 = {
    "BGI": "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG",
    "ILLUMINA": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
}


def collect_sequences(filepath):
    sequences = []
    with filepath.open("rt") as handle:
        for line in handle:
            if " Consensus: " in line:
                _, consensus = line.split()

                sequences.append({"consensus": consensus, "kmers": []})
                continue

            match = re.search(r" [0-9]: ([ACGT]+) ", line)
            if match is not None:
                (kmer,) = match.groups()

                sequences[-1]["kmers"].append(kmer)

    return sequences


def match_seq(read, adapter, max_errors=1):
    score = 0
    for idx, (nt_1, nt_2) in enumerate(zip(read, adapter), start=1):
        if nt_1 == nt_2:
            score += 1
        elif max_errors <= 0:
            break
        else:
            max_errors -= 1

    return score


def best_match(read, adapters, min_diff=2, max_errors=1):
    best_score = 0
    best_name = "-"

    for name, sequence in adapters.items():
        score = match_seq(read, sequence, max_errors=max_errors)
        if score - best_score > min_diff:
            best_score = score
            best_name = name
        elif abs(score - best_score) <= min_diff:
            best_score = max(score, best_score)
            best_name = "-"

    return best_name, best_score


def best_candidate(data, adapters, min_diff=2, max_errors=1, min_score=9):
    best_name, best_score = best_match(
        data["consensus"],
        adapters,
        min_diff=min_diff,
        max_errors=max_errors,
    )

    # No need to check kmers if they cannot match better
    if best_score >= KMER_LENGTH + min_diff and best_score >= min_score:
        return best_name, best_score

    assert all(len(kmer) == KMER_LENGTH for kmer in data["kmers"]), data
    for name, sequence in adapters.items():
        kmer = sequence[:KMER_LENGTH]
        score = KMER_LENGTH

        if kmer in data["kmers"]:
            if score - best_score > min_diff:
                best_score = score
                # Lowercase indicates that the kmer was the best match
                best_name = name.lower()
            elif abs(score - best_score) <= min_diff:
                # Collaborating matches are expected
                if name == best_name:
                    # Lowercase indicates that the kmer was the best match
                    if score > best_score:
                        best_name = best_name.lower()
                # Different kmer than the consensus matched
                elif name != best_name:
                    best_name = "-"

                best_score = max(score, best_score)

    if best_score < min_score:
        best_name = "-"

    return best_name, best_score


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--min-overlap", default=15, type=int)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    for filepath in sorted(args.files):
        name = filepath.name.split(".", 1)[0]
        sequences = collect_sequences(filepath)

        name_1, score_1 = best_candidate(sequences[0], ADAPTER_1)
        name_2, score_2 = best_candidate(sequences[1], ADAPTER_2)

        print(f"{name}\t{name_1}\t{score_1}\t{name_2}\t{score_2}")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
