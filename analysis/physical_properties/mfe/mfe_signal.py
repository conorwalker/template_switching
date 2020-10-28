#!/usr/bin/env python3

"""
Copyright (C) 2020 EMBL - European Bioinformatics Institute
Contact: goldman@ebi.ac.uk, cwalker@ebi.ac.uk

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import argparse

from sys import argv
from sh import RNAfold  # , Fold, mfold


def parse_args(args):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", help="Output directory.")
    parser.add_argument("--seq", help="Input sequence.")
    parser.add_argument("--id", help="ID for output filename and tsv tagging.")
    parser.add_argument("--window", help="""Number of bases in sliding mfe
                        window (default=20).""", default=50)
    return parser.parse_args()


def calculate_gc(seq):
    """
    Calculate GC content of a nucleotide sequence.
    """
    a, c, g, t = 0, 0, 0, 0
    for i in seq.upper():
        if i == "A":
            a += 1
        elif i == "C":
            c += 1
        elif i == "G":
            g += 1
        elif i == "T":
            t += 1
    try:
        return float((g+c)) / (a+c+g+t)
    except ZeroDivisionError:
        return 0


def mfe_sliding_window(pos_i, line_seq, win_size, outfile):
    """
    Slide over sequence fragments and scan windows for minimum free energy.
    """
    # Nearest neighbour parameters for ViennaRNA - ships with ViennaRNA
    dna_parameters = "dna_mathews2004.par"

    # Retrieve sliding window size
    window_size = int(win_size)
    step_size = int(win_size) - 1

    # Generate windows of sequences upstream/downstream
    seq_chunks = [line_seq[i:i+window_size] for i in range(len(line_seq)-step_size)]

    list_mfes = []

    # Calculate mfes and gc content for overlapping sequence windows
    for seq in seq_chunks:
        # ViennaRNA RNAfold mfes
        output = RNAfold(noLP=True, gquad=True, noPS=True, noconv=True,
                         P=dna_parameters, noGU=True, _in=seq)
        mfe = output.split("\n")[1].partition(" ")[2].partition("(")[2].partition(")")[0]

        # Calculate chunk GC content
        seq_gc = calculate_gc(seq)
        if "N" in seq.upper():
            contains_n = 1
        else:
            contains_n = 0
        list_mfes.append((mfe.strip(), seq_gc, contains_n))

    # Write to downstream/upstream files
    with open(outfile, "w+") as upstream_file:
        for i in list_mfes:
            to_write = str(pos_i) + "\t" + str(i[0]) + "\t" + str(i[1]) + \
                       "\t" + str(i[2]) + "\n"
            upstream_file.write(to_write)
    return None


def process_mfe_seq(mfe_args):
    """
    Processes each line in the known event file, calls other functions.
    """
    input_sequence = mfe_args.seq
    event_id = mfe_args.id
    window_size = mfe_args.window

    if mfe_args.outdir[-1] == "/":
        output_filename = str(mfe_args.outdir) + str(event_id) + ".tsv"
    else:
        output_filename = str(mfe_args.outdir) + "/" + str(event_id) + ".tsv"

    mfe_sliding_window(event_id, input_sequence, window_size, output_filename)
    return None


def main():
    # Process command line arguments for mfe scanning
    mfe_args = parse_args(argv[1:])
    process_mfe_seq(mfe_args)


if __name__ == "__main__":
    main()
