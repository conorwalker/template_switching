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

import tempfile
import subprocess

from sys import argv
from random import randint
from Bio import AlignIO


def main():
    """
    Generate a distribution of EPO alignment-wide unidirectional probabilities
    to produce a per-base alignment quality threshold.
    """

    # path to unidirectional pairHMM code
    unidirectional_pairhmm = ""

    fi_path = argv[1]
    current_align = AlignIO.read(fi_path, "fasta")
    homo_seq, pan_seq, gor_seq = "", "", ""
    for record in current_align:
        if "homo" in record.id:
            homo_seq = str(record.seq)
        elif "troglodytes" in record.id:
            pan_seq = str(record.seq)
        elif "gorilla" in record.id:
            gor_seq = str(record.seq)

    homo_len = len(homo_seq)
    pan_len = len(pan_seq)
    gor_len = len(gor_seq)

    homopan_len = min(homo_len, pan_len)
    homogor_len = min(homo_len, gor_len)
    pangor_len = min(pan_len, gor_len)

    for i in range(1, 21):
        homopan_coord = randint(500, homopan_len-500)
        homogor_coord = randint(500, homogor_len-500)
        pangor_coord = randint(500, pangor_len-500)

        homopan_1 = homo_seq[homopan_coord:homopan_coord+300]
        homopan_2 = pan_seq[homopan_coord:homopan_coord+300]

        homogor_1 = homo_seq[homogor_coord:homogor_coord+300]
        homogor_2 = gor_seq[homogor_coord:homogor_coord+300]

        pangor_1 = pan_seq[pangor_coord:pangor_coord+300]
        pangor_2 = gor_seq[pangor_coord:pangor_coord+300]

        with tempfile.NamedTemporaryFile(prefix=fi_path[26:56], dir="/scratch/") as tf:
            homopan_cont = ">seq1\n{}\n>seq2\n{}".format(homopan_1, homopan_2)
            tf.write(homopan_cont.encode())
            tf.flush()
            with open("homopan_alignment_logprobs.txt", "a+") as outfi:
                subprocess.run([unidirectional_pairhmm, "--pair", tf.name],
                               stdout=outfi)

        with tempfile.NamedTemporaryFile(prefix=fi_path[26:56], dir="/scratch/") as tf:
            homopan_cont = ">seq1\n{}\n>seq2\n{}".format(homogor_1, homogor_2)
            tf.write(homopan_cont.encode())
            tf.flush()
            with open("homogor_alignment_logprobs.txt", "a+") as outfi:
                subprocess.run([unidirectional_pairhmm,
                                "--pair", tf.name,
                                "--divergence", "0.016"],
                               stdout=outfi)

        with tempfile.NamedTemporaryFile(prefix=fi_path[26:56], dir="/scratch/") as tf:
            homopan_cont = ">seq1\n{0}\n>seq2\n{1}".format(pangor_1, pangor_2)
            tf.write(homopan_cont.encode())
            tf.flush()
            with open("pangor_alignment_logprobs.txt", "a+") as outfi:
                subprocess.run([unidirectional_pairhmm,
                                "--pair", tf.name,
                                "--divergence", "0.016"],
                               stdout=outfi)


if __name__ == "__main__":
    main()
