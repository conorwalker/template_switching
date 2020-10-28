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

import os
import numpy as np
from sys import argv
from sh import RNAfold


def main():
    if not os.path.exists("all_output_ws{0}".format(argv[1])):
        os.makedirs("all_output_ws{0}".format(argv[1]))
    if not os.path.exists("mean_output_ws{0}".format(argv[1])):
        os.makedirs("mean_output_ws{0}".format(argv[1]))

    if argv[2][-1] == "0":
        infi = argv[2][:-1]
    else:
        infi = argv[2]

    with open("ws" + argv[1] + "/" + infi + ".tsv", "r") as seq_fi:
        seqs = [i.strip().split("\t")[1] for i in seq_fi.readlines()]

    file_lines = []
    mfes = []
    for seq in seqs:
        output = RNAfold(noLP=True, gquad=True, noPS=True, noconv=True,
                         P="dna_mathews2004.par", noGU=True, _in=seq)
        mfe = output.split("\n")[1].partition(" ")[2].partition("(")[2].partition(")")[0]
        file_lines.append("{0}\t{1}\t{2}\n".format(argv[2], mfe, seq))
        mfes.append(float(mfe))

    with open("all_output_ws{0}/random_seq_MFEs_gc_{1}.tsv".format(argv[1], argv[2]), "w") as outfi:
        outfi.write("GC_content\tMFE\trandom_sequence\n")
        for fl in file_lines:
            outfi.write(fl)

    mean_mfe = np.mean(mfes)
    with open("mean_output_ws{0}/average_mfe_gc_{1}.tsv".format(argv[1], argv[2]), "w") as outfi:
        outfi.write("{0}\t{1}\n".format(argv[2], mean_mfe))


if __name__ == "__main__":
    main()
