#!/usr/bin/env python3
import argparse
import os
import sys
import warnings

from shlex import quote
from pathlib import Path

import ruamel.yaml


def safe_load(stream):
    yaml = ruamel.yaml.YAML(typ="safe", pure=True)
    yaml.version = (1, 1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ruamel.yaml.error.MantissaNoDotYAML1_1Warning)

        return yaml.load(stream)


def filter_options(dd):
    for key, value in dd.items():
        if key not in ("Options", "Prefixes", "Genomes"):
            yield key, value


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("yaml", type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output", default="output")
    parser.add_argument("--adapter1")
    parser.add_argument("--adapter2")

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    os.makedirs(args.output, exist_ok=True)

    data = safe_load(args.yaml)
    for group, samples in filter_options(data):
        for sample, libraries in filter_options(samples):
            for library, lanes in filter_options(libraries):
                for lane, tmpl in lanes.items():
                    if "{Pair}" in tmpl:
                        output = os.path.join(
                            args.output, "{}_{}_{}.txt".format(sample, library, lane)
                        )

                        if os.path.exists(output):
                            print(
                                f"INFO: Skipping {output!r}; already exists ..",
                                file=sys.stderr,
                            )
                            continue

                        command = [
                            "AdapterRemoval",
                            "--identify-adapters",
                            "--threads",
                            str(args.threads),
                            "--file1",
                            quote(tmpl.replace("{Pair}", "1")),
                            "--file2",
                            quote(tmpl.replace("{Pair}", "2")),
                        ]

                        if args.adapter1:
                            command += ["--adapter1", quote(args.adapter1)]

                        if args.adapter2:
                            command += ["--adapter2", quote(args.adapter2)]

                        command += [
                            ">",
                            quote(output + ".tmp"),
                            "&&",
                            "mv",
                            "-v",
                            quote(output + ".tmp"),
                            quote(output),
                        ]

                        print(" ".join(command))
                    else:
                        print(
                            "WARNING: Skipping SE lane {} > {} > {} > {}".format(
                                group, sample, library, lane
                            ),
                            file=sys.stderr,
                        )

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
