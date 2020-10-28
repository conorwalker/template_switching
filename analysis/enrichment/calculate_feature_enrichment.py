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
import os
import numpy as np

from random import shuffle
from sys import argv


def generate_random_coords(len_to_sample):
    """
    Generate random GRCh38 coordinates, checking they do not intersect
    with gapped regions of the genome.
    """
    chrom_lens = "GRCh38.p13.chrom_lengths.noY.tsv"
    gap_file = ""  # path to file containing gaps in GRCh38
    bed_random = subprocess.Popen(["bedtools",
                                   "random",
                                   "-l",
                                   "0",
                                   "-n",
                                   str(int(len_to_sample)*2),
                                   "-g",
                                   chrom_lens],
                                  stdout=subprocess.PIPE)
    bed_gap_inter = subprocess.Popen(["bedtools",
                                      "intersect",
                                      "-v",
                                      "-a",
                                      "stdin",
                                      "-b",
                                      gap_file],
                                     stdin=bed_random.stdout,
                                     stdout=subprocess.PIPE)

    bed_sort_out = [i.decode("utf-8").strip() for i in bed_gap_inter.stdout]
    bed_sort_write = bed_sort_out[:int(len_to_sample)]
    shuffle(bed_sort_write)

    with open("/tmp/tstout_unsorted", "w") as outfi:
        outfi.write("\n".join(bed_sort_write))

    with open("/tmp/tstout", "w") as outfi:
        subprocess.Popen(["bedtools", "sort", "-i", "/tmp/tstout_unsorted"],
                         stdout=outfi)
    return


def feature_permutations(feature_bed, len_to_sample, out_path, n_perms):
    """
    Randomly generate n_perms coordinates and check for intersection size with
    feature bed file.
    """
    n_intersect_all = []

    with open(out_path, "w+") as out_f:
        for i in range(n_perms):
            generate_random_coords(len_to_sample)
            bed_intersect = subprocess.Popen(["bedtools",
                                              "intersect",
                                              "-a",
                                              "/tmp/tstout",
                                              "-b",
                                              feature_bed,
                                              "-wa"],
                                             stdout=subprocess.PIPE)
            bed_intersect_out = [i.decode("utf-8").strip() for i in bed_intersect.stdout]
            n_intersect = len(bed_intersect_out)
            n_intersect_all.append(n_intersect)
            out_f.write(str(n_intersect) + "\n")
    return n_intersect_all


def get_event_intersect(feature_bed):
    """
    Get intersection size of feature file with gold-standard events.
    """
    bed_intersect = subprocess.Popen(["bedtools",
                                      "intersect",
                                      "-a",
                                      argv[1],
                                      "-b",
                                      feature_bed,
                                      "-wa"],
                                     stdout=subprocess.PIPE)
    bed_intersect_out = [i.decode("utf-8").strip() for i in bed_intersect.stdout]
    n_intersect = len(bed_intersect_out)
    return n_intersect


def empirical_p_value(event_n, permutations, n_perms):
    """
    Calculate empirical p-value enriched intersections.
    """
    r = 0
    for i in permutations:
        if i >= event_n:
            r += 1
    return (r+1)/(n_perms+1)


def empirical_p_value_neg(event_n, permutations, n_perms):
    """
    Calculate empirical p-value for depleted intersections.
    """
    r = 0
    for i in permutations:
        if i <= event_n:
            r += 1
    return (r+1)/(n_perms+1)


def test_feature(len_to_sample, bedfile, outfile):
    """
    Generate summary statistics and associated p-values for
    enrichment/depletion of events with a specific genomic feature.
    """
    constant = 0.001
    out_path = "random_permutations/" + outfile
    event_intersection_n = get_event_intersect(bedfile)
    n_permutations = int(argv[2])
    permutations = feature_permutations(bedfile,
                                        len_to_sample,
                                        out_path,
                                        n_permutations)
    perm_fold_changes = []

    for perm in permutations:
        pfc = np.log2((event_intersection_n+constant) / (perm+constant))
        perm_fold_changes.append(pfc)
    mean_log2_fold_change = np.log2(event_intersection_n+constant) - np.log2(np.mean(permutations)+constant)
    std_log2_fold_change = np.std(perm_fold_changes, dtype=np.float64)

    if mean_log2_fold_change < 0:
        p_val = empirical_p_value_neg(event_intersection_n,
                                      permutations,
                                      n_permutations)
    else:
        p_val = empirical_p_value(event_intersection_n,
                                  permutations,
                                  n_permutations)
    return [outfile.split(".")[0],
            mean_log2_fold_change,
            std_log2_fold_change,
            p_val,
            event_intersection_n,
            np.mean(permutations)]


def main():
    # make output directory
    if not os.path.isdir("random_permutations"):
        os.makedirs("random_permutations")

    # determine number of samples per run
    human_desc_events = argv[1]
    with open(human_desc_events, "r") as f:
        ts_events = [i.strip() for i in f.readlines()]
    len_to_sample = str(len(ts_events))

    # perform tests
    enhancers = test_feature(len_to_sample,
                             "computationally_predicted_enhancers/enhancer_atlas_enhancers.40celltissues.merged.lifted.bed",
                             "enhancer_atlas_enhancers_{0}.txt".format(argv[1]))
    
    super_enhancers = test_feature(len_to_sample,
                                   "super_enhancers/all_hg19_bed.lifted.merged.gt5.bed",
                                   "super_enhancers_{0}.txt".format(argv[1]))

    HARs = test_feature(len_to_sample,
                        "HARs/HARs_hg38lift_merged.bed",
                        "HARs_{0}.txt".format(argv[1]))

    PARs = test_feature(len_to_sample,
                        "PARs/PARs_hg38lift_merged.bed",
                        "PARs_{0}.txt".format(argv[1]))

    promoters = test_feature(len_to_sample,
                             "gencode_annotations/processed/promoters.requiring.protein_coding.bed",
                             "promoters_{0}.txt".format(argv[1]))

    TSS = test_feature(len_to_sample,
                       "gencode_annotations/processed/transcription_start_sites.requiring.protein_coding.bed",
                       "TSS_{0}.txt".format(argv[1]))

    intergenic = test_feature(len_to_sample,
                              "gencode_annotations/processed/intergenic.bed",
                              "intergenic_{0}.txt".format(argv[1]))

    introns = test_feature(len_to_sample,
                           "gencode_annotations/processed/introns.requiring.protein_coding_transcript.bed",
                           "introns_{0}.txt".format(argv[1]))

    lnc_rnas = test_feature(len_to_sample,
                            "gencode_annotations/processed/lncrna.requiring.transcript.merged.bed",
                            "lncrnas_{0}.txt".format(argv[1]))

    pseudogenes = test_feature(len_to_sample,
                               "gencode_annotations/processed/pseudogenes.requiring.transcript.merged.bed",
                               "pseudogenes_{0}.txt".format(argv[1]))

    UTRs = test_feature(len_to_sample,
                        "gencode_annotations/processed/UTR.requiring.protein_coding.merged.bed",
                        "UTRs_{0}.txt".format(argv[1]))

    exons = test_feature(len_to_sample,
                         "gencode_annotations/processed/exons.requiring.protein_coding.merged.bed",
                         "exons_{0}.txt".format(argv[1]))

    coding_regions = test_feature(len_to_sample,
                                  "gencode_annotations/processed/cds.requiring.protein_coding.merged.bed",
                                  "cds_{0}.txt".format(argv[1]))

    TFBS = test_feature(len_to_sample,
                        "transcription_factor_binding_sites/encRegTfbsClusteredWithCells.hg38.gt200.sorted.merged.gt4.cut.bed",
                        "transcription_factor_binding_sites_{0}.txt".format(argv[1]))

    full_list = [
        enhancers, super_enhancers, HARs,
        promoters, TSS, intergenic, introns,
        lnc_rnas, pseudogenes, UTRs, exons,
        coding_regions, TFBS, PARs
    ]

    outfi = "{0}_genomic_feature_df.tsv".format(argv[1].split(".")[0])

    with open(outfi, "w+") as f:
        f.write("feature\tmean_log2_foldchange\tstd_log2_foldchange\tp_val\tevent_intersect_n\tpermutation_mean_n\n")
        for feature in full_list:
            f.write(str(feature[0]) + "\t" + str(feature[1]) + "\t" +
                    str(feature[2]) + "\t" + str(feature[3]) + "\t" +
                    str(feature[4]) + "\t" + str(feature[5]) + "\n")


if __name__ == "__main__":
    main()
