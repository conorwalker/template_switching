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

import random
import os
import numpy as np
from sh import RNAfold
from sys import argv


def generate_random_gc_seq(seq_len, gc_perc):
    gc = ["G", "C"]
    at = ["A", "T"]
    random_seq = []
    required_gc = int(gc_perc * seq_len)
    for _ in range(required_gc):
        random_seq.append(random.choice(gc))
    for _ in range(seq_len - required_gc):
        random_seq.append(random.choice(at))
    random_seq = list(random_seq)
    random.shuffle(random_seq)
    return "".join(random_seq)


def main():
    random.seed(42)

    # make output directories
    if not os.path.exists("all_output"):
        os.makedirs("all_output")
    if not os.path.exists("mean_output"):
        os.makedirs("mean_output")

    target_gc = float(argv[1])
    if target_gc > 1:
        raise ValueError
    target_seq_len = 50
    n_iters = 10000

    file_lines = []
    mfes = []

    for n in range(n_iters):
        random_seq = generate_random_gc_seq(target_seq_len, target_gc)
        output = RNAfold(noLP=True, gquad=True, noPS=True, noconv=True,
                         P="dna_mathews2004.par", noGU=True, _in=random_seq)
        mfe = output.split("\n")[1].partition(" ")[2].partition("(")[2].partition(")")[0]
        file_lines.append("{0}\t{1}\t{2}\n".format(target_gc, mfe, random_seq))
        mfes.append(float(mfe))

    with open("all_output/random_seq_MFEs_gc_{0}.tsv".format(target_gc), "w") as outfi:
        outfi.write("GC_content\tMFE\trandom_sequence\n")
        for fl in file_lines:
            outfi.write(fl)

    mean_mfe = np.mean(mfes)
    with open("mean_output/average_mfe_gc_{0}.tsv".format(target_gc), "w") as outfi:
        outfi.write("{0}\t{1}\n".format(target_gc, mean_mfe))


if __name__ == "__main__":
    main()
