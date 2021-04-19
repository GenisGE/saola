#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import csv
import sys

from pathlib import Path


def read_table(filepath):
    with filepath.open(newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')



class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("files", nargs="+", type=Path)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    samples = {}
    for filepath in args.files:
        for row in read_table(filepath):
            table = samples
            for key in ('sample', 'library', 'run'):
                table = table.setdefault(row[key], {})

            table[row["pair"]] = row

    for sample, libraries in sorted(samples.items()):
        print(f"{sample!r}:")
        print(f"  {sample!r}:")

        for library, runs in sorted(libraries.items()):
            print(f"    {library!r}:")
            
            for run, data in sorted(runs.items()):
                keys = sorted(data)

                if keys == ["0"]:
                    fqpath = data['0']['fqpath']
                elif keys == ["1", "2"]:
                    fpath_1 = data["1"]["fqpath"]
                    fpath_2 = data["2"]["fqpath"]
                    if fpath_1 == fpath_2:
                        raise RuntimeError(data)

                    consensus = []
                    for ca, cb in zip(fpath_1, fpath_2):
                        if ca == cb:
                            consensus.append(ca)
                        else:
                            consensus.append("{Pair}")
                    
                    fqpath = "".join(consensus)
                else:
                    raise RuntimeError(data)
                
                print(f"      {run!r}: {fqpath!r}")
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
