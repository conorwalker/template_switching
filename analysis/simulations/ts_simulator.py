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
import subprocess

from Bio import SeqIO
from shutil import copyfile
from sys import argv


# global variables
DISTRIBUTION_PATH = "coord_distances.csv"
np.random.seed(42)


def load_distribution_data(dist_file, dist_path=DISTRIBUTION_PATH):
    """
    Load switch point distributions into a dataframe.
    """
    return pd.read_csv(os.path.join(dist_path, dist_file),
                       sep="\t",
                       header=0,
                       index_col=False)


def load_alignment_block(sim_fi):
    """
    Load an alignment FASTA.
    """
    copyfile(sim_fi, "ancestral.fasta")
    with open(sim_fi, "r") as sim_file:
        read_record = SeqIO.read(sim_file, "fasta")
    return read_record.id, str(read_record.seq)


def chunk_base_frequencies(sequence):
    """
    Get the base frequencies of a sequence.
    """
    sequence = sequence.upper()
    A = sequence.count("A") / float(len(sequence))
    C = sequence.count("C") / float(len(sequence))
    G = sequence.count("G") / float(len(sequence))
    T = sequence.count("T") / float(len(sequence))
    return [A, C, G, T]


def comp(seq):
    """
    Return the complement of an input nucleotide sequence.
    """
    return str.translate(seq, str.maketrans("ACGT", "TGCA"))


def replace_ancestral_id(ancestral_seq="ancestral.fasta"):
    """
    Change FASTA header.
    """
    contents = SeqIO.read(ancestral_seq, "fasta")
    with open(ancestral_seq, "w") as f:
        f.write(">tempseq\n")
        f.write(str(contents.seq))


def make_sim_directories():
    """
    Make output directory for simulated sequences.
    """
    cwd = os.getcwd()
    if not os.path.exists(cwd + "/simulated_sequences/"):
        os.makedirs("simulated_sequences")


def create_sim_fasta(chr_seq, chrom):
    """
    Sample a large chromosome fasta file and create a small fasta file for
    template switch simulation.
    """
    while True:
        start_coord = np.random.choice(len(chr_seq))
        subseq = chr_seq[start_coord:start_coord+1000]
        if "N" not in subseq and subseq.isupper():
            break
    seq_head = "GRCh38.{0}.{1}-{2}".format(str(chrom), str(start_coord),
                                           str(start_coord + 1000))
    with open("ancestral.fasta", "w") as fasta_file:
        fasta_file.write(">tempseq\n")
        fasta_file.write(subseq[:1000])
    return seq_head, subseq


def create_control_file(time, sim_seq_len):
    """
    Creates a new control.txt file for indelible.
    """
    try:
        os.remove("control.txt")
    except FileNotFoundError:
        pass
    l1 = "[TYPE] NUCLEOTIDE 1\n"
    l2 = "[MODEL]    hky_indels\n"
    l3 = "[submodel] HKY 2.1\n"
    l4 = "[indelmodel] POW 1.5 500\n"
    l5 = "[indelrate] 0.14\n"
    l6 = "[statefreq] 0.3 0.2 0.3 0.2\n\n"
    l7 = "[TREE] single_spec_tree   (tempseq:{0});\n".format(time)
    l8 = "[PARTITIONS] tempseq\n"
    l9 = "[single_spec_tree hky_indels {0}]\n\n".format(sim_seq_len)
    l10 = "[EVOLVE] tempseq 1 simulated_seq\n"
    with open("control.txt", "w") as f:
        f.writelines([l1, l2, l3, l4, l5, l6, l7, l8, l9, l10])


def main():
    divergence = float(argv[1]) / 100
    make_sim_directories()

    chr_dirs = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"
    chr_dirs = chr_dirs.split(",")

    # set working directory
    base_chrom_dir = ""
    # set indelible path
    indelible_path = ""
    # set nw path
    needleman_path = ""

    for chrom in chr_dirs:
        # read chromosome sequence
        with open(base_chrom_dir + chrom + ".fa", "rU") as fasta_fi:
            chr_seq = str(SeqIO.read(fasta_fi, "fasta").seq)
        # keep track of number of sequences simulated for this chromosome
        sim_file_count = 0
        # create a small fasta file by sampling the chromosome, also contains
        # the sequence and final id
        while sim_file_count != 1000:
            divide_divergence_by = np.random.choice(range(0,1000),1)[0] / 1000
            # load data
            coord_df = load_distribution_data(DISTRIBUTION_PATH)
            homo_seq = create_sim_fasta(chr_seq, chrom)
            coords_fragment_trimmed = coord_df[coord_df["two_three"] >= -100]

            # Determine location of template-switch and
            # time intervals of mutations
            ts_index = np.random.choice(range(200, len(homo_seq[1]) - 200, 1))
            time1 = divergence * divide_divergence_by
            time2 = divergence - time1

            time1 = "{:f}".format(time1)
            time2 = "{:f}".format(time2)

            # Edit control file for time1
            create_control_file(str(time1), str(len("".join(homo_seq[1]))))

            # Call indelible
            try:
                subprocess.call(indelible_path,
                                stderr=subprocess.STDOUT,
                                timeout=5)
            except:
                continue
            # Determine co-ordinates for switch-points
            switch_coords = []

            # Randomly select a row of the df to get coords from
            ts_coord_row = np.random.choice(len(coords_fragment_trimmed.index))

            point_1 = ts_index
            point_2 = ts_index + coords_fragment_trimmed.iloc[ts_coord_row]["one_two"]
            point_3 = point_2 + coords_fragment_trimmed.iloc[ts_coord_row]["two_three"]
            point_4 = ts_index + coords_fragment_trimmed.iloc[ts_coord_row]["one_four"]
            switch_coords = ([point_1, point_2+1, point_3, point_4])

            # Retrieve new seq from time1 mutations, add template switch
            try:
                interim_seq = SeqIO.read("simulated_seq.fas", "fasta")
            except ValueError:
                continue

            seq_to_edit = list(str(interim_seq.seq))
            ts_fragment = seq_to_edit[switch_coords[2]:switch_coords[1]]

            if len(ts_fragment) == 0:
                continue

            ts_fragment = [comp(base) for base in ts_fragment][::-1]
            del seq_to_edit[switch_coords[0]+1:switch_coords[3]]
            seq_to_edit.insert(switch_coords[0]+1, ts_fragment)
            seq_to_edit = [j for i in seq_to_edit for j in i]

            # Update sequence for second round of indelible
            with open("ancestral.fasta", "w") as f:
                f.write(">tempseq" + "\n")
                f.write("".join(seq_to_edit))

            create_control_file(str(time2), str(len("".join(seq_to_edit))))
            # second rounds of snps & indels
            try:
                subprocess.call(indelible_path,
                                stderr=subprocess.STDOUT,
                                timeout=5)
            except:
                continue

            try:
                simulated_sequence = str(SeqIO.read("simulated_seq.fas",
                                                    "fasta").seq)
            except:
                continue

            # create a single string of the two sequences for alignment
            # to pass as STDIN to subprocess module
            try:
                p = subprocess.Popen([needleman_path,
                                      simulated_sequence,
                                      str(homo_seq[1])],
                                     stdout=subprocess.PIPE)
            except OSError:
                continue

            align_out = p.communicate()[0]
            aligned_seqs = str(align_out).split("\\n")
            simulated_seq = aligned_seqs[0][2:]
            original_seq = aligned_seqs[1]

            new_seq_file = "ts_simulated_chrom_" + str(chrom) + \
                           "_fraglen_" + str(len(ts_fragment)) + \
                           "_time_" + str(time1) + "_len_" + \
                           str(len(seq_to_edit)) + "_seqid_" + \
                           str(sim_file_count) + ".fa"

            file_to_scan = "simulated_sequences/" + new_seq_file

            sw1, sw2 = str(switch_coords[0]), str(switch_coords[1])
            sw3, sw4 = str(switch_coords[2]), str(switch_coords[3])

            sim_head = ">simulated" + str(sim_file_count) + "." + \
                       str(chrom) + "_time_" + str(time1) + \
                       "_fraglen_" + str(len(ts_fragment)) + \
                       "_len_" + str(len(seq_to_edit)) + \
                       "_" + sw1 + "-" + sw2 + "-" + sw3 + "-" + sw4

            with open(file_to_scan, "w") as new_file:
                new_file.write(sim_head + "\n")
                new_file.write(simulated_seq + "\n")
                new_file.write(">" + str(homo_seq[0]) + "\n")
                new_file.write(original_seq + "\n")
            sim_file_count += 1


if __name__ == "__main__":
    main()
