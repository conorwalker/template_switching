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


import subprocess
import argparse
import operator

from sys import argv
from os import mkdir, path, remove
from Bio import AlignIO


def get_args(argv):
    """
    Parse command line arguments. Requires an input multiple alignment in FASTA
    format, optionally specify if old filtering scheme should be used.
    """
    parser = argparse.ArgumentParser(prog="event_compendium")
    parser.add_argument("--filename",
                        type=str,
                        default="a",
                        help="Input fasta file.")
    parser.add_argument("--logprob",
                        type=float,
                        default=9.0,
                        help="Use only log-probability and large deletions to filter events. (default)")
    parser.add_argument("--deletions",
                        type=int,
                        default=50,
                        help="Maximum number of deletions corrected by a template switch event.")
    parser.add_argument("--masking",
                        type=int,
                        default=0,
                        help="Boolean, allow masked bases in the 2->3 fragment? 0=False, 1=True")
    parser.add_argument("--nuc_sum",
                        type=int,
                        default=4,
                        help="Required sum of unique nucleotides in the 2->3 fragment (1-4)")
    return parser.parse_args(argv)


def parse_name(finame):
    """
    Determine which species are present in the MSA.
    """
    name_split = finame.split("-")
    if name_split[1][:5] == "pongo":
        return 0
    else:
        name_split[0] = name_split[0].split("/")[1]
        finames = [i.split(".")[0] for i in name_split]
        return finames


def generate_finame_dic(finame):
    """
    Determine filenames for split alignments.
    """
    finame = finame.split("/")[1][:-4]
    finame_dic = {}
    for spec in finame.split("-"):
        finame_dic[spec.split(".")[0]] = spec
    return finame_dic


def generate_fasta(head1, head2, seq1, seq2, finame):
    """
    Generate alignment fasta files.
    """
    with open(finame, "w+") as f:
        f.write(">" + head1 + "\n")
        f.write(seq1 + "\n")
        f.write(">" + head2 + "\n")
        f.write(seq2)
    return


def remove_gap_cols(seq1, seq2):
    """
    Remove gap only columns from the pairwise alignment.
    """
    keep_cols = []
    assert len(seq1) == len(seq2)
    for i in range(len(seq1)):
        if seq1[i] == "-" and seq2[i] == "-":
            continue
        else:
            keep_cols.append(i)
    seq1 = "".join([seq1[i] for i in keep_cols])
    seq2 = "".join([seq2[i] for i in keep_cols])
    return seq1, seq2


def convert_to_global_coords(fields, gapfile):
    """
    Converts the coordinates of the pairwise FPA output into coordinates
    of the original MSA before gap only columns were removed.
    """
    with open(gapfile, "r") as fi:
        csvlines = fi.readlines()
    qry_gaps = [int(i) for i in csvlines[0].strip().split(",")]
    ref_gaps = [int(i) for i in csvlines[1].strip().split(",")]
    try:
        for i in range(3, 10):
            if i < 6:
                fields[i] = str(int(fields[i]) + qry_gaps[int(fields[i])])
            else:
                fields[i] = str(int(fields[i]) + ref_gaps[int(fields[i])])
    except IndexError:
        print("gapfi:", gapfile)
        print("i:", i)
        print("fields:", fields)
        print("len fields:", len(fields))
        print("len qry gaps:", len(qry_gaps))
        print("len ref gaps:", len(ref_gaps))
    return fields


def generate_global_file(scanfile, gapfi):
    """
    Convert coordinates of scanned lines to global coordinates.
    """
    if not path.exists(scanfile):
        return
    with open(scanfile, "r") as scnfi:
        for line in scnfi:
            fields = line.strip().split(",")
            with open(".".join(scanfile.split(".")[:-1]) + ".global.csv", "a+") as outfi:
                if fields[0] == "chrom" or fields[0] == "scan finished":
                    outfi.write(line)
                else:
                    outfi.write(",".join(convert_to_global_coords(fields, gapfi)) + "\n")
    return


def realign_short_fragments(scanfi, line, diverg, fasfi):
    """
    Realign (using --print-file) template switch events with short
    L->1, 2->3, 4->R fragments.
    Return the line and a boolean indicator of no longer short or still short.
    """
    # create temp fi containing the single csv line for re-alignment
    # using --print-file
    print("realigning...")
    tmp_align_fi = "/tmp/temp_realign_{0}".format(scanfi.split("/")[-1])
    with open(tmp_align_fi, "w") as realign_f:
        realign_f.write(line)

    # re-align, output to temp file
    tmp_out_fi = "/tmp/temp_phmm_out_{0}".format(scanfi.split("/")[-1])
    with open(tmp_out_fi, "w") as outfi:
        subprocess.run(["",  # TSA pairHMM path
                        "--pair", fasfi,
                        "--divergence", diverg,
                        "--print-file", tmp_align_fi,
                        "--scan-flank", "400",
                        "--long-output"], stdout=outfi)
    try:
        # get new line from pairHMM output
        with open(tmp_out_fi, "r") as expanded_pairHMM_out:
            out_lines = [i.strip() for i in expanded_pairHMM_out.readlines()]
            pairHMM_out = out_lines[[i+1 for i, j in enumerate(out_lines) if j.startswith('chrom,clus_start')][0]]

        # file cleanup
        remove(tmp_align_fi)
        remove(tmp_out_fi)

        # check fragment lengths
        pairHMM_fields = pairHMM_out.split(",")
        f1 = int(pairHMM_fields[30])
        f2 = int(pairHMM_fields[31])
        f3 = int(pairHMM_fields[32])
        if all(i >= 2 for i in [f1, f2, f3]):
            return [pairHMM_out + "\n", 1]
        else:
            return [0, 0]
    # hit the end of the input FASTA with the expanded boundary, causing
    # blank pairHMM output. skip resolving this event.
    except IndexError:
        return [0, 0]


def filter_scanned(scanfi, filterfi, args):
    """
    Filter FPA output files.
    """
    filter_fi_list = []
    with open(scanfi, "r") as scnfi:
        for line in scnfi:
            fields = line.strip().split(",")
            if fields[0] == "chrom" or fields[0] == "scan finished":
                continue
            point2, point3 = int(fields[7]), int(fields[8])
            down_id, up_id = float(fields[10]), float(fields[12])
            masked = int(fields[16])
            sum_nuc = int(fields[20])
            fwd_score, local_prob = float(fields[25]), float(fields[26])
            clus_del = int(fields[23])
            if point2 - point3 > args.length and \
               down_id >= args.downstream and \
               up_id >= args.upstream and \
               masked == args.masking and \
               sum_nuc == args.nuc_sum and \
               local_prob - fwd_score >= args.logprob and \
               clus_del <= args.deletions:
                filter_finame = filterfi + "{0}-{1}.csv".format(str(int(fields[3])),
                                                                str(int(fields[4])))
                filter_fi_list.append(filter_finame)
                with open(filter_finame, "w") as outfi:
                    outfi.write(line)
    return filter_fi_list


def logprob_scanned(scanfi, filterfi, args, per_base_logprob_thresh, diverg, fasfi):
    """
    Use minimal filtering for log-probability ratio threshold, max deletion cutoff,
    no masking and minimum number of bases for filtering events.
    """
    filter_fi_list = []
    filtered_events = []
    with open(scanfi, "r") as scnfi:
        for line in scnfi:
            fields = line.strip().split(",")
            if fields[0] == "chrom" or fields[0] == "scan finished":
                continue
            fwd_score, local_prob = float(fields[25]), float(fields[26])
            qry_len = float(fields[29])
            clus_del = int(fields[23])
            masked = int(fields[16])
            sum_nuc = int(fields[20])
            ts_logprob = float(fields[26])
            ts_per_base_logprob = (ts_logprob + 14.78259) / qry_len
            f1 = int(fields[30])
            f2 = int(fields[31])
            f3 = int(fields[32])

            if local_prob - fwd_score >= args.logprob and \
               masked == args.masking and \
               sum_nuc == args.nuc_sum and \
               ts_per_base_logprob >= per_base_logprob_thresh and \
               clus_del <= args.deletions:
                # check if all alignment fragments >= 2 in length
                if all(i >= 2 for i in [f1, f2, f3]):
                    filter_finame = filterfi + "{0}-{1}.csv".format(str(int(fields[3])),
                                                                    str(int(fields[4])))
                    filter_fi_list.append(filter_finame)
                    filtered_events.append(line)
                # if one fragment is v. short, try to re-align with expanded search space
                # if still short, continue, otherwise include as a filtered event
                else:
                    # output is line, boolean
                    realigned = realign_short_fragments(scanfi, line, diverg, fasfi)
                    fields = realigned[0].split(",")
                    if realigned[1] == 1:
                        if fields[0] == "chrom" or fields[0] == "scan finished":
                            continue
                        fwd_score, local_prob = float(fields[25]), float(fields[26])
                        qry_len = float(fields[29])
                        clus_del = int(fields[23])
                        masked = int(fields[16])
                        sum_nuc = int(fields[20])
                        ts_logprob = float(fields[26])
                        ts_per_base_logprob = ts_logprob / qry_len
                        if local_prob - fwd_score >= args.logprob and \
                           masked == args.masking and \
                           sum_nuc == args.nuc_sum and \
                           ts_per_base_logprob >= per_base_logprob_thresh and \
                           clus_del <= args.deletions:
                            # use a different filename to indicate that a larger scan-flank
                            # setting was used, will be needed later to print the correct
                            # output
                            filter_finame = filterfi + "{0}-{1}.large_flank.csv".format(str(int(fields[3])),
                                                                                        str(int(fields[4])))
                            filter_fi_list.append(filter_finame)
                            # realigned[0] replaces "line" above
                            filtered_events.append(realigned[0])

    # remove overlapping events, keep only the event with the largest cluster
    if len(filtered_events) > 1:
        cluster_lens = {}
        event_count = 0
        coord_list = []
        for event in filtered_events:
            e_vals = event.strip().split(",")
            e_vals = [float(i) for i in e_vals[1:]]
            cluster_len = e_vals[3] - e_vals[2]
            cluster_lens[event_count] = cluster_len
            points = sorted(e_vals[5:9])
            coord_list.append((int(points[0]), int(points[-1])))
            event_count += 1
        if len(set.intersection(*(set(range(start, finish+1)) for start, finish in coord_list))) != 0:
            largest_clust_key = max(cluster_lens.items(), key=operator.itemgetter(1))[0]
            with open(filter_finame, "w") as outfi:
                outfi.write(filtered_events[largest_clust_key])
        else:
            for line in filtered_events:
                with open(filter_finame, "w") as outfi:
                    outfi.write(line)
    else:
        for line in filtered_events:
            with open(filter_finame, "w") as outfi:
                outfi.write(line)
    return filter_fi_list


def main():
    args = get_args(argv[1:])

    # change to the path of the TSA pairHMM
    tsa_pairhmm_path = ""

    parsed_name = parse_name(args.filename)

    with open(args.filename, "r") as alignfi:
        msa = AlignIO.read(alignfi, "fasta")

    finame_dic = generate_finame_dic(args.filename)

    record_dic = {}

    for record in msa:
        record_dic[record.id.split(".")[0]] = [record.id, str(record.seq)]

    blockdir = args.filename.split("/")[1][:-4]

    if not path.exists("alignment_blocks/"):
        mkdir("alignment_blocks")
    if not path.exists("alignment_blocks/" + blockdir):
        mkdir("alignment_blocks/" + blockdir)

    pairs = [("homo_sapiens", "pan_troglodytes"),
             ("homo_sapiens", "gorilla_gorilla"),
             ("pan_troglodytes", "gorilla_gorilla")]

    # 20th percentile as determined using random alignment permutations
    hp_pb, hg_pb, pg_pb = -0.148574, -0.17974167, -0.18545034

    per_bases = {
        "homopan": hp_pb,
        "homogor": hg_pb,
        "pangor": pg_pb
    }

    for comp in pairs:
        trunc_name = comp[0].split("_")[0]
        trunc_name2 = comp[1].split("_")[0]
        if trunc_name == "gorilla":
            trunc_name = "gor"
        if trunc_name2 == "gorilla":
            trunc_name2 = "gor"
        species_dir = trunc_name + trunc_name2
        species_dir_swapped = trunc_name2 + trunc_name
        if comp[0] in parsed_name and comp[1] in parsed_name:
            if trunc_name == "homo" and trunc_name2 == "pan":
                diverg = "0.01"
            else:
                diverg = "0.016"
            mkdir("alignment_blocks_np/" + blockdir + "/" + species_dir)
            mkdir("alignment_blocks_np/" + blockdir + "/" + species_dir_swapped)
            first_head = record_dic[comp[0]][0]
            second_head = record_dic[comp[1]][0]
            nuc_seq, nuc_seq_swapped = remove_gap_cols(record_dic[comp[0]][1],
                                                       record_dic[comp[1]][1])

            finame1 = "alignment_blocks_np/{0}/{1}/{2}.fas".format(blockdir,
                                                                    species_dir,
                                                                    finame_dic[comp[0]])
            finame2 = "alignment_blocks_np/{0}/{1}/{2}.fas".format(blockdir,
                                                                    species_dir_swapped,
                                                                    finame_dic[comp[1]])

            generate_fasta(first_head, second_head, nuc_seq, nuc_seq_swapped, finame1)
            generate_fasta(second_head, first_head, nuc_seq_swapped, nuc_seq, finame2)

            scanfi1 = finame1[:-4] + ".scanned.csv"
            with open(scanfi1, "w") as outfi:
                subprocess.run([tsa_pairhmm_path,
                                "--scan", "--pair", finame1,
                                "--divergence", diverg], stdout=outfi)

            scanfi2 = finame2[:-4] + ".scanned.csv"
            with open(scanfi2, "w") as outfi:
                subprocess.run([tsa_pairhmm_path,
                                "--scan", "--pair", finame2,
                                "--divergence", diverg], stdout=outfi)

            filterfi1 = finame1[:-4] + ".filtered."
            filterfi2 = finame2[:-4] + ".filtered."

            global_gap_file = ".".join(args.filename.split("/")[-1].split(".")[:-1]) + ".csv"

            gapfi1 = species_dir + "_gaps/" + global_gap_file
            gapfi2 = species_dir_swapped + "_gaps/" + global_gap_file

            filtered1 = logprob_scanned(scanfi1, filterfi1, args, per_bases[species_dir],
                                        diverg, finame1)  # params passed for realignment
            filtered2 = logprob_scanned(scanfi2, filterfi2, args, per_bases[species_dir],
                                        diverg, finame2)

            for fi in filtered1:
                generate_global_file(fi, gapfi1)
            for fi in filtered2:
                generate_global_file(fi, gapfi2)

            generate_global_file(scanfi1, gapfi1)
            generate_global_file(scanfi2, gapfi2)


if __name__ == "__main__":
    main()
