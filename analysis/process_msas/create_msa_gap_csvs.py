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
import csv

from sys import argv, exit
from Bio import SeqIO, AlignIO


def get_finame(fir3, fi_split):
    for f in fi_split:
        if f[:3] == fir3:
            return f


def main():
    # Read in species names, can be truncated, 2 must be query sequence, 3 ref
    spec1 = str(argv[3])
    spec2 = str(argv[4])

    # parse output directory name
    outdir = str(argv[2])

    # create directory
    try:
        os.mkdir(outdir)
    except OSError:
        pass

    # Create list to store records to be written to file (2/3 matches)
    recs_to_write = []

    # Split file name at - to determine new filename
    finame_split = argv[1].split("-")
    finame_split[-1] = finame_split[-1][:-4]
    if spec1[:3] == "hom":
        finame1 = get_finame("hom", finame_split)
    elif spec1[:3] == "pan":
        finame1 = get_finame("pan", finame_split)
    elif spec1[:3] == "gor":
        finame1 = get_finame("gor", finame_split)

    if spec2[:3] == "hom":
        finame2 = get_finame("hom", finame_split)
    elif spec2[:3] == "pan":
        finame2 = get_finame("pan", finame_split)
    elif spec2[:3] == "gor":
        finame2 = get_finame("gor", finame_split)

    finame = "-".join((finame1, finame2))

    temp_rec_dic = {}

    # Retrieve records of interest to write
    with open(argv[1], "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id[:3] == spec1[:3]:
                temp_rec_dic["species1"] = record
            elif record.id[:3] == spec2[:3]:
                temp_rec_dic["species2"] = record

    recs_to_write.append(temp_rec_dic["species1"])
    recs_to_write.append(temp_rec_dic["species2"])

    # exit if alignment block in species 2 missing
    if len(recs_to_write) != 2:
        exit()

    # Create temp file which contains both sequences with gap-only columns
    # still present
    with open("/tmp/" + finame + ".tmp", "w") as temp_file:
        SeqIO.write(recs_to_write, temp_file, "fasta")

    # keep count of gaps in seq 1 (query)
    spec1_vect = []
    spec1_gapcount = 0
    spec2_vect = []
    spec2_gapcount = 0

    # read temp alignment
    aln = AlignIO.read("/tmp/" + finame + ".tmp", "fasta")
    aln_len = aln.get_alignment_length()

    # two separate loops required
    for i in range(aln_len):
        col_nucs = [sr.seq[i] for sr in aln]
        if col_nucs[0] == "-":
            spec1_gapcount += 1
            continue
        spec1_vect.append(spec1_gapcount)
    for i in range(aln_len):
        col_nucs = [sr.seq[i] for sr in aln]
        if col_nucs[1] == "-":
            spec2_gapcount += 1
            continue
        spec2_vect.append(spec2_gapcount)

    both_specs_vect = [spec1_vect, spec2_vect]

    # save gaps in seq 1 vector as csv file
    with open(os.path.join(outdir, finame+".csv"), "w") as csv_fi:
        writer = csv.writer(csv_fi)
        writer.writerows(both_specs_vect)

    # remove temp alignment with gap only cols in
    os.remove("/tmp/" + finame+".tmp")


if __name__ == "__main__":
    main()
