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

import pandas as pd
import subprocess

from os import scandir, makedirs, path
from Bio import AlignIO
from glob import glob

import direction_dicts as dd


def find_min_max_coord(list_of_lines):
    """
    Get minimum/maximum coordinate from switch points or mutation cluster.
    Used for checking overlap across pairwise comparisons.
    """
    min_coord, max_coord = 100000000, 0
    list_of_lines = [i[3:10] for i in list_of_lines]
    for sublist in list_of_lines:
        for i in sublist:
            if int(i) < min_coord:
                min_coord = int(i)
            if int(i) > max_coord:
                max_coord = int(i)
    return min_coord, max_coord


def switch_point_annotations(sps, align_len, min_coord, max_coord):
    """
    Generates switch point annotations (1,2,3,4) for printed output.
    """
    if len(sps) == 0:
        return ""
    sps = [int(i) for i in sps]
    spaces = list(" " * align_len)
    p1, p2, p3, p4 = sps[0], sps[1], sps[2], sps[3]
    spaces[p1] = "1"
    spaces[p2] = "2"
    spaces[p3] = "3"
    spaces[p4] = "4"
    return "".join(spaces)[min_coord-20:max_coord+20]


def cluster_annotation(clus_start, clus_end, align_len, min_coord, max_coord):
    """
    Generates cluster boundary annotations ([,]) for printed output.
    """
    spaces = list(" " * align_len)
    spaces[clus_start+2] = "["
    spaces[clus_end] = "]"
    return "".join(spaces)[min_coord-20:max_coord+20]


def print_align_seq(species_seq, start_pos, end_pos):
    """
    Generates subsequence for species in printed output.
    """
    if len(species_seq) == 0:
        return ""
    return species_seq[start_pos-20:end_pos+20]


def print_arrow(present, updown):
    """
    Adds an up/down arrow in the printed output to indicate direction.
    """
    if len(present) == 0:
        return "          "
    if updown == "up":
        return u"        \u2191 "
    else:
        return u"        \u2193 "


def resolve_start_and_end(finame_dic):
    """
    Resolves the start/end coordinates when checking for overlaps across
    multiple comparisons for events without a congruent direction.
    """
    resolve_start_coords, resolve_end_coords = [], []
    coord_list = []
    for fis in finame_dic:
        with open(fis, "r") as f:
            ficont = [i.strip() for i in f.readlines()]
        ficont = ficont[0].split(",")
        resolve_start_coords.append(int(ficont[3]))
        resolve_end_coords.append(int(ficont[4]))
        append_indices = [3, 4, 6, 7, 8, 9]
        for ai in append_indices:
            coord_list.append(int(ficont[ai]))
    # include switch points into overlap range
    return min(coord_list), max(coord_list)


def ils_interest_dirs(direction):
    """
    If directionality not established for some of the single set events, check
    the comparisons that would indicate incomplete lineage sorting before
    assigning the event as direction uninferrable.
    """
    dir_dic = {"1": "gorpan",
               "2": "gorhomo",
               "4": "panhomo",
               "6": "homopan"}
    return dir_dic[direction]


def determine_failure_reason(line, direction, comp_dir):
    """
    Checks the pairHMM output and looks for values that failed the various
    thresholds, records these as reasons for failure.
    """
    comp_dic = {"homopan": 1,
                "panhomo": 2,
                "homogor": 3,
                "gorhomo": 4,
                "pangor": 5,
                "gorpan": 6}
    inferred_direction = [int(i) for i in list(direction)]
    inferred_direction.append(comp_dic[comp_dir])
    inferred_direction = "".join([str(j) for j in sorted(inferred_direction)])
    line = line.strip().split(",")
    point2, point3 = int(line[7]), int(line[8])
    up_id, down_id = float(line[10]), float(line[12])
    masked = int(line[16])
    sum_nuc = int(line[20])
    fwd_logprob = float(line[25])
    ts_logprob = float(line[26])
    clus_del = int(line[23])
    qry_len = float(line[29])
    # factor out o1 -> M2 and M2 -> o3 transitions to get per-base
    # log probability, described in Methods
    ts_per_base_logprob = (ts_logprob + 14.78259) / qry_len

    rows = [direction, inferred_direction, "na", "na", "na", "na", "na",
            "na", "na", "na", "na"]
    if point2 - point3 <= 8:
        rows[2] = point2 - point3
    if up_id < 0.95:
        rows[3] = up_id
    if down_id < 0.95:
        rows[4] = down_id
    if masked != 0:
        rows[5] = "yes"
    if sum_nuc != 4:
        rows[6] = sum_nuc
    if ts_logprob - fwd_logprob <= 9:
        rows[7] = ts_logprob - fwd_logprob
    if clus_del >= 50:
        rows[8] = clus_del
    rows[10] = ts_per_base_logprob
    return rows


def scan_global_file(global_scanned, target_set, direction, comp_dir):
    """
    In the case that an event is significant in one comparison, but not
    congruent in direction across comparisons, the global coordinates of the
    significant event are used to scan for event output across all pairwise
    comparisons that would rectify the congruence.
    """
    with open(global_scanned, "r") as scan_fi:
        scan_lines = scan_fi.readlines()
    for line in scan_lines:
        try:
            clust_coords = [int(i) for i in line.strip().split(",")[3:5]]
            scan_coords = clust_coords
        except ValueError:
            continue
        try:
            assert len(clust_coords) > 1
        except AssertionError:
            continue
        qry_range = range(min(scan_coords), max(scan_coords))
        intersect_result = len(target_set.intersection(qry_range))
        if intersect_result != 0:
            fail_row = determine_failure_reason(line, direction, comp_dir)
            return fail_row
    return None


def inferred_direction(direction, inferred_dirs):
    """
    From the list of pairwise comparisons the event is found in, produce the
    inferred direction string (e.g. homopan, homogor becomes 13).
    """
    comp_dic = {"homopan": 1,
                "panhomo": 2,
                "homogor": 3,
                "gorhomo": 4,
                "pangor": 5,
                "gorpan": 6}
    inferred_direction = [int(i) for i in list(direction)]
    for i in inferred_dirs:
        inferred_direction.append(comp_dic[i])
    inferred_direction = "".join([str(j) for j in sorted(inferred_direction)])
    return inferred_direction


def combine_fail_reason(row1, row2):
    """
    Combine various reasons for failing to pass thresholds if failure was
    observed in 2 comparisons.
    """
    combined_row = row1
    for i, j in enumerate(row2):
        if row1[i] == "na" and row2[i] != "na":
            combined_row[i] = row2[i]
    return combined_row


def combine_fail_reason_three(row1, row2, row3):
    """
    Combine various reasons for failing to pass thresholds if failure was
    observed in 3 comparisons.
    """
    combined_row = row1
    for i, j in enumerate(row2):
        if row1[i] == "na" and row2[i] != "na":
            combined_row[i] = row2[i]
    for i, j in enumerate(row3):
        if row1[i] == "na" and row3[i] != "na":
            combined_row[i] = row2[i]
    return combined_row


def resolve_missing(direction, target_set, block_path):
    """
    Check for event across rectifying comparisons from the global coord pairHMM
    scan output. Try to identify inferred direction and reason for failure.
    """
    global_scanned, global_scanned1, global_scanned2 = None, None, None
    int_dirs = dd.determine_interest_dirs(direction)
    if direction == "3" or direction == "5":
        for comp_pair in int_dirs:
            for name in glob("{}/{}/*scanned.global.csv".format(block_path,
                                                                int_dirs[0])):
                global_scanned1 = name
            for name in glob("{}/{}/*scanned.global.csv".format(block_path,
                                                                int_dirs[1])):
                global_scanned2 = name
            if not global_scanned1 and not global_scanned2:
                return direction, [direction, "missing", "na", "na", "na",
                                   "na", "na", "na", "na", "missing", "na"]
            else:
                scan_result1 = scan_global_file(global_scanned1, target_set,
                                                direction, int_dirs[0])
                scan_result2 = scan_global_file(global_scanned2, target_set,
                                                direction, int_dirs[1])
                if scan_result1:
                    fail_row = scan_result1
                    inf_dir = inferred_direction(direction, [int_dirs[0]])
                    fail_row[1] = inf_dir
                    return inf_dir, fail_row
                elif scan_result2:
                    fail_row = scan_result2
                    inf_dir = inferred_direction(direction, [int_dirs[1]])
                    fail_row[1] = inf_dir
                    return inf_dir, fail_row
            return None

    if len(direction) == 1 or len(direction) == 3:
        for name in glob("{0}/{1}/*scanned.global.csv".format(block_path,
                                                              int_dirs)):
            global_scanned = name
        if not global_scanned:
            return direction, [direction, "missing", "na", "na", "na", "na",
                               "na", "na", "na", "missing", "na"]
        else:
            scan_result = scan_global_file(global_scanned, target_set,
                                           direction, int_dirs)

            if scan_result:
                fail_row = scan_result
                inf_dir = inferred_direction(direction, [int_dirs])
                fail_row[1] = inf_dir
                return inf_dir, fail_row

    elif direction in ["12", "34", "56"]:
        for comp_pair in int_dirs:
            for name in glob("{0}/{1}/*scanned.global.csv".format(block_path,
                                                                  comp_pair[0])):
                global_scanned1 = name
            for name in glob("{0}/{1}/*scanned.global.csv".format(block_path,
                                                                  comp_pair[1])):
                global_scanned2 = name
            if not global_scanned1 or not global_scanned2:
                return direction, [direction, "missing", "na", "na", "na",
                                   "na", "na", "na", "na", "missing", "na"]
            else:
                scan_result1 = scan_global_file(global_scanned1, target_set,
                                                direction, comp_pair[0])
                scan_result2 = scan_global_file(global_scanned2, target_set,
                                                direction, comp_pair[1])
                if scan_result1 and scan_result2:
                    fail_row = combine_fail_reason(scan_result1, scan_result2)
                    inf_dir = inferred_direction(direction,
                                                 [comp_pair[0], comp_pair[1]])
                    fail_row[1] = inf_dir
                    return inf_dir, fail_row

    # formerly contained "16" [C(H,G)] and "24" [H(C,G)], removed these as
    # ILS events can be considered congruent
    elif direction in ["15", "23", "36", "45", "26", "14"]:
        for name in glob("{0}/{1}/*scanned.global.csv".format(block_path,
                                                              int_dirs[0])):
            global_scanned1 = name
        for name in glob("{0}/{1}/*scanned.global.csv".format(block_path,
                                                              int_dirs[1])):
            global_scanned2 = name
        if not global_scanned1 or not global_scanned2:
            return direction
        else:
            scan_result1 = scan_global_file(global_scanned1, target_set,
                                            direction, int_dirs[0])
            scan_result2 = scan_global_file(global_scanned2, target_set,
                                            direction, int_dirs[1])
            if scan_result1 and scan_result2:
                fail_row = combine_fail_reason(scan_result1, scan_result2)
                inf_dir = inferred_direction(direction,
                                             [int_dirs[0], int_dirs[1]])
                fail_row[1] = inf_dir
                return inf_dir, fail_row
    # Before assigning unresolved, check relevant direction events for ILS
    if direction in ["1", "2", "4", "6"]:
        ils_int_dir = ils_interest_dirs(direction)
        for name in glob("{0}/{1}/*scanned.global.csv".format(block_path,
                                                              ils_int_dir)):
            global_scanned = name
        if not global_scanned:
            return direction, [direction, "missing", "na", "na", "na", "na",
                               "na", "na", "na", "missing", "na"]
        else:
            scan_result = scan_global_file(global_scanned, target_set,
                                           direction, ils_int_dir)
            if scan_result:
                fail_row = scan_result
                inf_dir = inferred_direction(direction, [ils_int_dir])
                fail_row[1] = inf_dir
                return inf_dir, fail_row
    return None


def generate_comparison_code(homo_seq, pan_seq, gor_seq):
    """
    Compares the three muscle aligned sequenced in a 3-way zip tuple,
    assigns a value V={0,1,3,5} along the multiple alignment length N.
    """
    code_list = []
    for a, b, c in zip(homo_seq, pan_seq, gor_seq):
        if a == b == c:
            code_list.append(0)
        elif a == b and a != c:
            code_list.append(1)
        elif a == c and b != c:
            code_list.append(3)
        elif b == c and a != b:
            code_list.append(5)
    return code_list


def muscle_align_segs(h, p, g):
    """
    Use muscle to align clusters in the case that initial resolution attempt
    was unsuccessful.
    """
    seqs = {"h": "", "p": "", "g": ""}
    muscle = "/nfs/research1/goldman/conor/tools/muscle3.8.31/muscle"
    with open("tmpfas.fa", "w") as f:
        f.write(">h\n{0}\n>p\n{1}\n>g\n{2}".format(h, p, g))
    try:
        musc_out = subprocess.check_output([muscle, "-in",
                                            "tmpfas.fa", "-quiet"])
        musc_out = musc_out.decode("utf-8").split("\n")[:-1]
    except subprocess.CalledProcessError:
        return h, p, g
    for i, line in enumerate(musc_out):
        if line.startswith(">h"):
            seqs["h"] = musc_out[i+1]
        if line.startswith(">p"):
            seqs["p"] = musc_out[i+1]
        if line.startswith(">g"):
            seqs["g"] = musc_out[i+1]
    longest_seq = max([len(v) for v in seqs.values()])
    for k, v in seqs.items():
        if len(seqs[k]) == 0:
            seqs[k] = "-" * longest_seq
    return seqs["h"], seqs["p"], seqs["g"]


def align_cluster_resolve(sequence_dic, direction, min_coord, max_coord):
    """
    Try to resolve event direction by aligning clusters using muscle.
    """
    homo_seq = sequence_dic["homo_sapiens"][min_coord+2:max_coord]
    pan_seq = sequence_dic["pan_troglodytes"][min_coord+2:max_coord]
    gor_seq = sequence_dic["gorilla_gorilla"][min_coord+2:max_coord]

    homo_seq, pan_seq, gor_seq = muscle_align_segs(homo_seq,
                                                   pan_seq,
                                                   gor_seq)

    muscle_align_segs(homo_seq, pan_seq, gor_seq)

    sim_code = generate_comparison_code(homo_seq, pan_seq, gor_seq)
    similarity_code = []
    for i in sim_code:
        if i != 0:
            similarity_code.append(i)
    if direction in ["3", "5", "12", "34", "56"]:
        needed_codes = dd.determine_needed_code(direction)
        for need in needed_codes:
            if not any(i != need for i in similarity_code):
                if direction == "3" and need == 1:
                    return "35"
                if direction == "3" and need == 5:
                    return "13"
                if direction == "5" and need == 1:
                    return "35"
                if direction == "5" and need == 3:
                    return "25"
                if direction == "12" and need == 3:
                    return "1256"
                if direction == "12" and need == 5:
                    return "1234"
                if direction == "34" and need == 5:
                    return "1234"
                if direction == "34" and need == 1:
                    return "3456"
                if direction == "56" and need == 3:
                    return "1256"
                if direction == "56" and need == 1:
                    return "3456"
    else:
        needed_code = dd.determine_needed_code(direction)
        if not any(i != needed_code for i in similarity_code):
            return dd.resolved_direction(str(direction))
    return [None, "\n".join([homo_seq, pan_seq, gor_seq]), similarity_code]


def annotate_fail(fail_row):
    """
    Annotate reason for failure with verbose string.
    """
    fail_dic = {}
    fail_dic["Masked:"] = fail_row[5]
    fail_dic["Sum of nucleotides:"] = fail_row[6]
    fail_dic["Log-probability ratio:"] = fail_row[7]
    fail_dic["Deletions in cluster:"] = fail_row[8]
    fail_dic["Per-base log-probability:"] = fail_row[10]
    fail_dic = {k: v for k, v in fail_dic.items() if v != "na"}
    return fail_dic


def determine_anc_desc_seqs(direction, anc_desc_dic, list_of_fis):
    """
    Take in a congruent direction, the dictionary to map this
    direction to ancestral and descendant sequences, and a list of
    filtered event files. From these files, construc the pairwise
    alignment FASTA file name, and extract the ancestral and
    descendant sequences from this. This will be used to generate
    the upstream/downstream sequences for secondary structure analysis.
    """
    # determine ancestor/descendant
    anc_desc = anc_desc_dic[direction]
    # generate fasta file names
    fas_fis = []
    for csv_fi in list_of_fis:
        if "large_flank" in csv_fi.split("."):
            fas_fis.append(".".join(csv_fi.split(".")[:-5]) + ".fas")
        else:
            fas_fis.append(".".join(csv_fi.split(".")[:-4]) + ".fas")
    # determine prefix for descendant sequences
    desc_prefix = anc_desc[0][0]
    # go through alignments, extract sequences
    for fas_fi in fas_fis:
        if fas_fi.split("/")[2][:3] == desc_prefix[:3]:
            pair_al = AlignIO.read(fas_fi, "fasta")
            for rec in pair_al:
                if rec.name.split(".")[0][:3] == desc_prefix[:3]:
                    desc_seq = str(rec.seq).replace("-", "")
                else:
                    anc_seq = str(rec.seq).replace("-", "")
            break
    return desc_seq, anc_seq


def determine_switch_points(anc_desc, list_of_fis, anc_or_desc):
    """
    Take in anc_desc, list of finames, a well as indicator for
    if the desired switch points are ancestral of descendant. For ancestor,
    return points 1-4 coordinates, for ref, just return point 1 coordinates.
    """
    desc_prefix = anc_desc[0][0]
    for csv_fi in list_of_fis:
        if csv_fi.split("/")[2][:3] == desc_prefix[:3]:
            csv_fi_local = csv_fi.split(".")
            csv_fi_local = ".".join(csv_fi_local[:-2]) + ".csv"
            with open(csv_fi_local, "r") as event_f:
                csv_read = event_f.readlines()[-1]
                csv_read_split = csv_read.strip().split(",")
            if anc_or_desc == "anc":
                p1 = int(csv_read_split[6])
                p2 = int(csv_read_split[7])
                p3 = int(csv_read_split[8])
                p4 = int(csv_read_split[9])
                return [p1, p2, p3, p4]
            elif anc_or_desc == "desc":
                return int(csv_read_split[5])


def determine_inf_switch_points(anc_desc, list_of_fis, anc_or_desc):
    """
    Take in anc_desc, list of finames, a well as indicator for
    if the desired switch points are ancestral of descendant. For ancestor,
    return points 1-4 coordinates, for ref, just return point 1 coordinates.
    """
    for desc_prefix in anc_desc[0]:
        for csv_fi in list_of_fis:
            if csv_fi.split("/")[2][:3] == desc_prefix[:3]:
                csv_fi_local = csv_fi.split(".")
                csv_fi_local = ".".join(csv_fi_local[:-2]) + ".csv"
                with open(csv_fi_local, "r") as event_f:
                    csv_read = event_f.readlines()[-1]
                    csv_read_split = csv_read.strip().split(",")
                if anc_or_desc == "anc":
                    p1 = int(csv_read_split[6])
                    p2 = int(csv_read_split[7])
                    p3 = int(csv_read_split[8])
                    p4 = int(csv_read_split[9])
                    return [p1, p2, p3, p4]
                elif anc_or_desc == "desc":
                    return int(csv_read_split[5])


def determine_point_ordering(anc_points):
    """
    Return indices of ancestral switch points.
    """
    points = [1, 2, 3, 4]
    return "".join([str(x) for _, x in sorted(zip([int(i)
                   for i in anc_points], points))])


def main():
    root_dir = "hominid_template_switch_output"
    # make output directories
    if not path.exists(root_dir):
        makedirs(root_dir)
    if not path.exists("{}/pairHMM_output".format(root_dir)):
        makedirs("{}/pairHMM_output".format(root_dir))
    if not path.exists("{}/directions".format(root_dir)):
        makedirs("{}/directions".format(root_dir))
    if not path.exists("{}/fail_dataframe".format(root_dir)):
        makedirs("{}/fail_dataframe".format(root_dir))
    if not path.exists("{}/bed_output".format(root_dir)):
        makedirs("{}/bed_output".format(root_dir))
    if not path.exists("{}/sequences".format(root_dir)):
        makedirs("{}/sequences".format(root_dir))
    if not path.exists("{}/troubleshooting".format(root_dir)):
        makedirs("{}/troubleshooting".format(root_dir))
    if not path.exists("{}/point_orders".format(root_dir)):
        makedirs("{}/point_orders".format(root_dir))

    # add header to bed files
    out_bed = "{}/bed_output/template_switch_events.bed".format(root_dir)
    with open(out_bed, "w") as bed_f:
        bed_f.write("chrom\tchromStart\tchromEnd\tname\torientation\tsetAnnotation\thumanAncOrDesc\teventType\n")

    # used to check fragment content later
    nucs = ["A", "T", "C", "G"]

    # various dictionaries for checking directions
    chrom_len_dic = dd.chrom_len_dic
    comp_dic = dd.comp_dic
    comp_translate = dd.comp_translate
    acceptable = dd.acceptable

    # keep track of events processed, used to assign event ID
    total_event_count = 0

    # lists to store directions
    congruent_directions = []
    filter_fail_directions = []
    aligned_cluster_directions = []
    unresolved_directions = []
    missing_directions = []

    # all directions (congruent, inferred and unresolved)
    final_direction_list = []

    # df to hold fail reasons per event
    fail_columns = ["direction", "inferred_direction", "ts_len", "up_id",
                    "down_id", "masked", "sum_nuc", "logprob_ratio",
                    "clus_del", "missing", "logprob_perbase"]
    fail_df = pd.DataFrame(columns=fail_columns)

    # path to pairHMM
    pairHMM = "tools/tsa_pairhmm/tsa_pairhmm"

    with scandir("alignment_blocks_np") as block_dir:
        # process each alignment block
        for block in block_dir:
            sequence_dic = {"homo_sapiens": "",
                            "pan_troglodytes": "",
                            "gorilla_gorilla": "",
                            "pongo_abelii": ""}
            # get name of 4-way species concat from directory name
            current_block = block.name
            # Read alignment
            alignment_fasta_path = "fas4/{0}.fas".format(current_block)
            # read the 4-way alignment
            with open(alignment_fasta_path, "r") as fasfi:
                alignment_read = AlignIO.read(fasfi, "fasta")
            # Get alignment length for cluster output later
            alignment_len = len(alignment_read[0].seq)
            # Get seq for each species if present
            for record in alignment_read:
                if "homo_sapiens" in record.id:
                    sequence_dic["homo_sapiens"] = str(record.seq)
                if "pan_troglodytes" in record.id:
                    sequence_dic["pan_troglodytes"] = str(record.seq)
                if "gorilla_gorilla" in record.id:
                    sequence_dic["gorilla_gorilla"] = str(record.seq)
                if "pongo_abelii" in record.id:
                    sequence_dic["pongo_abelii"] = str(record.seq)

            # keep track of number of events that passed filtering from
            # find_events.py
            filter_count = 0
            with scandir("alignment_blocks_np/" + current_block) as curr_dir:
                df_rows = []
                # process each of the 6 pairwise blocks separately, ie.
                # chimp <-> human, human <-> gorilla, chimp <-> gorilla
                for comparison in curr_dir:
                    with scandir("alignment_blocks_np/{0}/{1}".format(current_block, comparison.name)) as comp:
                        # get all csvs that passed filters using the global
                        # coordinates ie. global coords = those that correspond
                        # to multiple alignment before removing pairwise
                        # gap-only columns
                        for filterfi in comp:
                            # for each event csv that passes filters, get
                            # event CSV and event FASTA
                            if "filtered" in filterfi.name and "global" in filterfi.name:
                                spl_name = filterfi.name.split(".")
                                # now indexes based on position of "filtered" in file name,
                                # ape regions tagged with CAB* had extra fields, messed up printing
                                if "large_flank" in filterfi.name:
                                    local_filterfi = ".".join(spl_name[:spl_name.index("filtered")+3]) + ".csv"
                                else:
                                    local_filterfi = ".".join(spl_name[:spl_name.index("filtered")+2]) + ".csv"

                                local_filterfi = "alignment_blocks_np/{0}/{1}/{2}".format(current_block,
                                                                                          comparison.name,
                                                                                          local_filterfi)
                                # account for "chrom" names with "." in name
                                pair_fasfi = spl_name[:5]
                                pair_fasfi = ".".join(list(filter(lambda b: b != "filtered", pair_fasfi))) + ".fas"
                                filter_count += 1
                                pair_fasfiname_format = "alignment_blocks_np/{0}/{1}/{2}".format(current_block,
                                                                                                 comparison.name,
                                                                                                 pair_fasfi)
                                finame_format = "alignment_blocks_np/{0}/{1}/{2}".format(current_block,
                                                                                         comparison.name,
                                                                                         filterfi.name)
                                with open(finame_format, "r") as fi:
                                    # retrieve some of the values from the
                                    # event CSV corresponding to coordinates
                                    # that can be used to check for overlaps
                                    # across species comparisons in global
                                    # coordinates
                                    csv_vals = [int(i) for i in fi.readlines()[-1].strip().split(",")[1:5]]
                                    # get some meta-information about the
                                    # event and store associated files
                                    df_rows.append([comparison.name,
                                                    finame_format,
                                                    csv_vals[1],
                                                    csv_vals[2],
                                                    csv_vals[3],
                                                    pair_fasfiname_format,
                                                    local_filterfi])

                # initialise DFs and dicts for processing all events in
                # each multiple alignment block
                df_cols = ["pair_comp", "path", "chrom_start", "start_coord",
                           "end_coord", "pair_fasfi", "local_filterfi"]
                df_ref = pd.DataFrame(df_rows, columns=df_cols)
                df_pop = pd.DataFrame(df_rows, columns=df_cols)
                unique_events = {}
                unique_events_finames = {}
                unique_events_local_finames = {}
                unique_event_cluster_start_end = {}
                unique_events_pair_fastas = {}

                # keep track of event number (not final assigned event
                # number, used for tracking temporarily)
                event_n = 0
                dir_matched = []
                for index, row in df_ref.iterrows():
                    event_n += 1
                    # rename for brevity below
                    comp, event_path, chrom_start, start, end = row.pair_comp, row.path, row.chrom_start, row.start_coord, row.end_coord
                    # don't process dir if already seen
                    if str(row.path) in dir_matched:
                        continue

                    # go through rows for this alignment block
                    for ind_pop, row_pop in df_pop.iterrows():
                        # if the start coords across comparisons match, add
                        # to dictionaries
                        if row_pop.pair_comp == comp and (row_pop.start_coord == start or row_pop.end_coord == end or row_pop.chrom_start == chrom_start):
                            unique_events["event_{0}".format(event_n)] = [comp_dic[comp]]
                            unique_events_finames["event_{0}".format(event_n)] = [event_path]
                            unique_events_local_finames["event_{0}".format(event_n)] = {}
                            unique_events_local_finames["event_{0}".format(event_n)][comp] = row_pop.local_filterfi
                            unique_events_pair_fastas["event_{0}".format(event_n)] = {}
                            unique_events_pair_fastas["event_{0}".format(event_n)][comp] = row_pop.pair_fasfi
                            df_pop.drop(ind_pop)
                        # if start coords don't match, check for overlaps
                        else:
                            # get range of coords that events are allowed to overlap
                            # to be called as same event
                            qry_range = range(start, end)
                            target_range = range(row_pop.start_coord, row_pop.end_coord)
                            set_range = set(qry_range)
                            # if some overlap observed, process adding events
                            if len(set_range.intersection(target_range)) != 0:
                                # either create dict value or append to it if exists for this event
                                # unique_events stores pairwise comparisons observed in
                                # unique_events_finames stores the filenames associated with these
                                try:
                                    unique_events["event_{0}".format(event_n)].append(comp_dic[row_pop.pair_comp])
                                    unique_events_finames["event_{0}".format(event_n)].append(row_pop.path)
                                except KeyError:
                                    unique_events["event_{0}".format(event_n)] = [comp_dic[row_pop.pair_comp]]
                                    unique_events_finames["event_{0}".format(event_n)] = [row_pop.path]
                                comp_to_add = row_pop.pair_fasfi.split("/")[2]
                                try:
                                    unique_events_local_finames["event_{0}".format(event_n)][comp_to_add] = row_pop.local_filterfi
                                except KeyError:
                                    unique_events_local_finames["event_{0}".format(event_n)] = {}
                                    unique_events_local_finames["event_{0}".format(event_n)][comp_to_add] = row_pop.local_filterfi
                                try:
                                    unique_events_pair_fastas["event_{0}".format(event_n)][comp_to_add] = row_pop.pair_fasfi
                                except KeyError:
                                    unique_events_pair_fastas["event_{0}".format(event_n)] = {}
                                    unique_events_pair_fastas["event_{0}".format(event_n)][comp_to_add] = row_pop.pair_fasfi

                                # record cluster range of matched events
                                start = min(start, row_pop.start_coord)
                                end = max(end, row_pop.end_coord)
                                unique_event_cluster_start_end["event_{0}".format(event_n)] = [start, end]

                                # store path for blocks with matched events
                                dir_matched.append(str(row_pop.path))

                # process all unique events found
                for k, v in unique_events.items():
                    # store the ascertained direction of the event for printing
                    # and checking congruence
                    final_direction = "".join([str(i) for i in sorted(v)])
                    # keep track of unique events processed
                    total_event_count += 1

                    # list to store all cluster boundaries for printing
                    clus_start, clus_end = [], []
                    comp_switch_point_dic = {"homopan": [],
                                             "panhomo": [],
                                             "homogor": [],
                                             "gorhomo": [],
                                             "pangor": [],
                                             "gorpan": []}

                    # dictionary to identify descendant and ancestral species
                    # from the final, congruent direction
                    # first tuple = descendants, second = ancestors
                    anc_desc_dic = {
                        "13": [["homo"], ["pan", "gor"]],
                        "25": [["pan"], ["homo", "gor"]],
                        "46": [["gor"], ["pan", "homo"]],
                        "35": [["homo", "pan"], ["gor"]],
                        "1234": [["homo"], ["pan", "gor"]],
                        "1256": [["pan"], ["homo", "gor"]],
                        "3456": [["homo", "pan", "gor"], ["homo", "pan", "gor"]],
                        # ILS C(H,G)
                        "16": [["homo", "gor"], ["pan"]],
                        # ILS H(C,G)
                        "24": [["pan", "gor"], ["homo"]]
                    }

                    set_annotation = ""
                    # determine if larger scan-flank parameter needed for
                    # the TSA pairHMM
                    large_flank = 0
                    acceptable_line_dic = {}

                    ############################################
                    # PROCESS EVOLUTIONARILY CONSISTENT EVENTS #
                    ############################################
                    if final_direction in acceptable:
                        if final_direction not in ["1234", "1256", "3456"]:
                            set_annotation = "gold_standard"
                        else:
                            set_annotation = "ambiguous_placement"
                        # per congruent event, record one set of switch
                        # point ordering as well as reversibility
                        # for reporting later
                        if final_direction not in ["1234", "1256", "3456"]:
                            ancestral_switch_point_output = determine_switch_points(anc_desc_dic[final_direction], unique_events_finames[k], "anc")
                            ancestral_switch_point_order = determine_point_ordering(ancestral_switch_point_output)
                            # if congruent direction in 4 comparisons, it is
                            # reversible retrieve the descendant switch point
                            # ordering
                            if len(unique_events_finames[k]) == 4:
                                reversible_switch_point_order = "1"
                            else:
                                reversible_switch_point_order = "0"
                            with open("hominid_template_switch_output/point_orders/congruent_points.tsv", "a+") as point_fi:
                                point_fi.write(ancestral_switch_point_order + "\t" + reversible_switch_point_order + "\n")

                        # get path of alignment block
                        block_path = "/".join(unique_events_finames[k][0].split("/")[:2])

                        # store this event as congruent
                        congruent_directions.append(final_direction)

                        # don't try to resolve, print blocks, clusters and
                        # pairHMM output for each comparison with event
                        final_direction_list.append(final_direction)

                        # check here for ancestra/descendant state of human
                        anc_desc = anc_desc_dic[final_direction]
                        if "homo" in anc_desc[0]:
                            # human_desc_bool = 1
                            anc_or_desc = "descendant"
                        elif "homo" in anc_desc[1]:
                            # human_anc_bool = 1
                            anc_or_desc = "ancestral"
                        # if a 3456 event is found, resolve manually
                        # register human as descendant for downstream
                        # processing, although impossible to tell with
                        # 3 pairwise comparisons
                        if final_direction == "3456":
                            # human_desc_bool = 1
                            anc_or_desc = "ambiguous"

                        for fis in unique_events_finames[k]:
                            # gets the non global file for neg coord generation
                            fi_for_bed = ".".join(fis.split(".")[:-2]) + ".csv"
                            # determine if found on negative strand, changes
                            # how coordinates are parsed, given that all
                            # coordinates should be presented on + strand for
                            # associating with genomic features and retrieving
                            # sequences from databases
                            if "large_flank" in fis.split("."):
                                large_flank = 1
                                pair_fasfi = ".".join(fis.split(".")[:-5]) + ".fas"
                            else:
                                pair_fasfi = ".".join(fis.split(".")[:-4]) + ".fas"

                            # get switch points for this comparison for
                            # printing in output
                            with open(fis, "r") as event_f:
                                # get last line as above
                                csv_read = event_f.readlines()[-1]
                                csv_read_split = csv_read.strip().split(",")
                                # add csv output to dic or append if exists
                                try:
                                    acceptable_line_dic[k].append(csv_read_split)
                                except KeyError:
                                    acceptable_line_dic[k] = [csv_read_split]
                                # get switch point order
                                comp_switch_point_dic[fis.split("/")[2]] = csv_read_split[6:10]
                                if len(csv_read_split) == 2:
                                    csv_read_split = csv_read_split[-1]
                                clus_start.append(int(csv_read_split[3]))
                                clus_end.append(int(csv_read_split[4]))

                            with open(fi_for_bed, "r") as local_event_fi:
                                local_csv_read = local_event_fi.readlines()[-1]

                            # store output of pairHMM for congruent events
                            with open("hominid_template_switch_output/pairHMM_output/congruent_pairHMM_output.csv", "a+") as f:
                                f.write(str(total_event_count) + "," + final_direction + "," + local_csv_read)

                            # check if human bed file needs writing to, or
                            # if the human ancestral or descendant coords
                            # have already been recorded
                            anc_points = determine_switch_points(anc_desc_dic[final_direction], unique_events_finames[k], "anc")
                            event_type = determine_point_ordering(anc_points)
                        # GET START/END CLUSTER POSITIONS IN HUMAN GENOME
                        fas4_fi = "fas4/" + fis.split("/")[1] + ".fas"
                        with open(fas4_fi, "r") as f:
                            fas4_lines = f.readlines()
                            for i_line, line in enumerate(fas4_lines):
                                if line.startswith(">homo"):
                                    fas4_homo_header = line.strip()[1:]
                                    fas4_homo_seq = fas4_lines[i_line+1].strip()
                                    break
                        human_chrom = fas4_homo_header.split(".")[1]
                        human_block_start = int(fas4_homo_header.split(".")[2])
                        if fas4_homo_header.split("_")[-1] == "neg":
                            strand = "reverse"
                            start_gap_count = fas4_homo_seq[:max(clus_end)].count("-")
                            end_gap_count = fas4_homo_seq[:(min(clus_start)+2)].count("-")
                            chrom_len = chrom_len_dic[human_chrom]
                            bed_start = chrom_len - human_block_start - max(clus_end)+1 + start_gap_count
                            bed_end = chrom_len - human_block_start - min(clus_start)-1 + end_gap_count
                        else:
                            strand = "forward"
                            start_gap_count = fas4_homo_seq[:(min(clus_start)+2)].count("-")
                            end_gap_count = fas4_homo_seq[:max(clus_end)].count("-")
                            bed_start = human_block_start + min(clus_start)+2 - start_gap_count
                            bed_end = human_block_start + max(clus_end) - end_gap_count

                        bed_line = "\t".join(
                            [str(i) for i in [human_chrom,
                             bed_start,
                             bed_end,
                             total_event_count,
                             strand,
                             set_annotation,
                             anc_or_desc,
                             event_type + "\n"]]
                        )
                        with open("hominid_template_switch_output/bed_output/template_switch_events.bed", "a") as bed_f:
                            bed_f.write(bed_line)
                        # function which takes the congruent direction,
                        # anc_desc_dic to determine anc and descendant seq,
                        # and the list of file names for this event to output
                        # the descendant and ancestral sequences for use in
                        # local sequence signal analysis
                        if final_direction not in ["1234", "1256", "3456"]:
                            desc_seq, anc_seq = determine_anc_desc_seqs(final_direction, anc_desc_dic, unique_events_finames[k])
                            anc_points = determine_switch_points(anc_desc_dic[final_direction], unique_events_finames[k], "anc")
                            desc_point = determine_switch_points(anc_desc_dic[final_direction], unique_events_finames[k], "desc")
                            with open(unique_events_finames[k][0], "r") as chrom_f:
                                chrom = chrom_f.readlines()[-1].strip().split(",")[0]

                            # get +/- 2k sequences around switch points
                            anc_p1 = anc_seq[anc_points[0]-500:anc_points[0]+551]
                            anc_p2 = anc_seq[anc_points[1]-500:anc_points[1]+551]
                            anc_p3 = anc_seq[anc_points[2]-500:anc_points[2]+551]
                            anc_p4 = anc_seq[anc_points[3]-500:anc_points[3]+551]
                            desc_p1 = desc_seq[desc_point-500:desc_point+551]

                            # get +/- 10k sequences around point 1
                            anc_p1_10k = anc_seq[anc_points[0]-10000:anc_points[0]+10001]
                            desc_p1_10k = desc_seq[desc_point-10000:desc_point+10001]

                            # get +/- 100 sequences around point 1
                            anc_p1_150 = anc_seq[anc_points[0]-150:anc_points[0]+151]
                            desc_p1_150 = desc_seq[desc_point-150:desc_point+151]

                            # meta information for writing to sequence file
                            event_twothree_len = (anc_points[1] - anc_points[2]) + 1
                            event_type = determine_point_ordering(anc_points)
                            event_direction = final_direction

                            out_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(total_event_count, chrom, event_twothree_len, event_direction, event_type, desc_p1, anc_p1, anc_p2, anc_p3, anc_p4)
                            out_line_10k = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(total_event_count, chrom, event_twothree_len, event_direction, event_type, desc_p1_10k, anc_p1_10k)
                            out_line_150 = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(total_event_count, chrom, event_twothree_len, event_direction, event_type, desc_p1_150, anc_p1_150)

                            with open("hominid_template_switch_output/sequences/congruent_event_seqs.tsv", "a+") as fasta_f:
                                fasta_f.write(out_line)
                            with open("hominid_template_switch_output/sequences/congruent_event_seqs_10k.tsv", "a+") as fasta_f:
                                fasta_f.write(out_line_10k)
                            with open("hominid_template_switch_output/sequences/congruent_event_seqs_150.tsv", "a+") as fasta_f:
                                fasta_f.write(out_line_150)

                        # PRINT OUTPUT OF CONGRUENT EVENTS
                        min_coord, max_coord = find_min_max_coord(acceptable_line_dic[k])
                        event_title = "Event {0}".format(total_event_count)
                        if total_event_count != 1:
                            print("\n")
                        print("=" * len(event_title))
                        print(event_title)
                        print("=" * len(event_title))
                        if final_direction in ["1234", "1256"]:
                            print("Phylogenetic placement: ambiguous, species tree or ILS consistent")
                        elif final_direction in ["3456"]:
                            print("Phylogenetic placement: ambiguous, species tree consistent")
                        elif final_direction in ["16", "24"]:
                            print("Phylogenetic placement: ILS consistent")
                        else:
                            print("Phylogenetic placement: species tree consistent")
                        print("Significant comparisons: {0}".format(final_direction))
                        print("Detected in comparisons: {0}".format(final_direction))
                        print("\nMultiple sequence alignment:\n")
                        print("    Human:", print_align_seq(sequence_dic["homo_sapiens"], min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["panhomo"], "down") + switch_point_annotations(comp_switch_point_dic["panhomo"], alignment_len, min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["homopan"], "up") + switch_point_annotations(comp_switch_point_dic["homopan"], alignment_len, min_coord, max_coord))
                        print("    Chimp:", print_align_seq(sequence_dic["pan_troglodytes"], min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["gorpan"], "down") + switch_point_annotations(comp_switch_point_dic["gorpan"], alignment_len, min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["pangor"], "up") + switch_point_annotations(comp_switch_point_dic["pangor"], alignment_len, min_coord, max_coord))
                        print("  Gorilla:", print_align_seq(sequence_dic["gorilla_gorilla"], min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["homogor"], "down") + switch_point_annotations(comp_switch_point_dic["homogor"], alignment_len, min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["gorhomo"], "up") + switch_point_annotations(comp_switch_point_dic["gorhomo"], alignment_len, min_coord, max_coord))
                        print("    Human:", print_align_seq(sequence_dic["homo_sapiens"], min_coord, max_coord))
                        print("\n  Cluster:" + cluster_annotation(min(clus_start), max(clus_end), alignment_len, min_coord, max_coord))
                        print("\n")
                        print("Orangutan:", print_align_seq(sequence_dic["pongo_abelii"], min_coord, max_coord))
                        for comp_k in comp_switch_point_dic.keys():
                            if len(comp_switch_point_dic[comp_k]) != 0:
                                # try:
                                diverg = dd.pairwise_divergence(unique_events_pair_fastas[k][comp_k])
                                if large_flank == 1:
                                    scan_flank_par = "100"
                                else:
                                    scan_flank_par = "40"
                                pairHMM_out = subprocess.check_output([pairHMM, "--pair", unique_events_pair_fastas[k][comp_k],
                                                                       "--print-file", unique_events_local_finames[k][comp_k],
                                                                       "--divergence", diverg,
                                                                       "--scan-flank", scan_flank_par])
                                title = comp_translate[comp_k]
                                print("\n\n{0}".format(title))
                                print("-" * len(title))
                                print(pairHMM_out.decode().lstrip())

                    
                    #################################
                    # PROCESS PARTIALLY PASS EVENTS #
                    #################################
                    else:
                        # track if event has been resolved
                        filter_resolved = False
                        cluster_align_resolved = False
                        block_missing_resolved = False
                        resolved_bool = "unresolved"
                        phylo_place = ""
                        unacceptable_line_dic = {}

                        # try to find event or align clusters and determine direction,
                        # print alignment block output with cluster annotation and sp annotation,
                        # as well as inferred direction and reason(s) for failure
                        # first check if all gapped or block missing

                        # record pairHMM output for events without a congruent direction
                        # used to determine most stringent reason for failure
                        # most of the time this is the LPR threshold
                        unres_fi = open("hominid_template_switch_output/troubleshooting/noncongruent_event_output.csv", "a+")
                        unres_fi.write(str(final_direction))
                        for fis in unique_events_finames[k]:
                            for scanned_dirs in unique_events_finames[k]:
                                split_scanned_dir = scanned_dirs.split(".")
                                if "large_flank" in split_scanned_dir:
                                    large_flank = 1
                            with open(fis, "r") as event_f:
                                csv_read = event_f.readlines()[0]
                                csv_read_split = csv_read.strip().split(",")
                                try:
                                    unacceptable_line_dic[k].append(csv_read_split)
                                except KeyError:
                                    unacceptable_line_dic[k] = [csv_read_split]
                                comp_switch_point_dic[fis.split("/")[2]] = csv_read_split[6:10]
                                clus_start.append(int(csv_read_split[3]))
                                clus_end.append(int(csv_read_split[4]))
                            unres_fi.write("\t" + str(fis))
                        unres_fi.write("\n")
                        unres_fi.close()

                        # get global coord range for event
                        min_coord, max_coord = find_min_max_coord(unacceptable_line_dic[k])
                        resolve_start, resolve_end = resolve_start_and_end(unique_events_finames[k])
                        target_set = set(range(resolve_start, resolve_end))

                        # get path to alignment block for current event
                        block_path = "/".join(unique_events_finames[k][0].split("/")[:2])
                        block_missing = None

                        # prelim checking for failure reason of sequence in the multiple alignment
                        # don't process orangutan, only included for manual outgroup inspection
                        for blank_k in sequence_dic.keys():
                            if "pongo" in blank_k:
                                continue
                            # store sequence for species
                            seq_of_interest = print_align_seq(sequence_dic[blank_k], min_coord, max_coord)
                            # assess if sequence is all gaps, N, 0 len, or contains no nucs
                            if seq_of_interest == "N" * len(seq_of_interest):
                                block_missing = True
                            if seq_of_interest == "n" * len(seq_of_interest):
                                block_missing = True
                            if seq_of_interest == "-" * len(seq_of_interest):
                                block_missing = True
                            if len(seq_of_interest) == 0:
                                block_missing = True
                            if not any(x in seq_of_interest for x in nucs):
                                block_missing = True
                        inferred_direction = ""

                        # resolved if confirmed alignment block missing
                        if block_missing:
                            inferred_direction = final_direction
                            resolved_bool = "inferred (missing alignment block)"
                            phylo_place = "ambiguous, missing/gapped/masked alignment block"
                            missing_directions.append(final_direction)
                            block_missing_resolved = True

                        # if alignment block, try to resolve
                        else:
                            try:
                                resolve_status = resolve_missing(final_direction, target_set, block_path)
                                inferred_direction, fail_row = resolve_status[0], resolve_status[1]
                                #fail_annotation = annotate_fail(fail_row)
                                filter_fail_directions.append(inferred_direction)
                                resolved_bool = "inferred (failed filter)"
                                if inferred_direction in ["1234","1256"]:
                                    phylo_place = "ambiguous, species tree or ILS consistent"
                                elif inferred_direction in ["3456"]:
                                    phylo_place = "ambiguous, species tree consistent"
                                elif inferred_direction in ["16","24"]:
                                    phylo_place = "ILS consistent"
                                else:
                                    phylo_place = "species tree consistent"
                                fail_df.loc[total_event_count] = fail_row
                                filter_resolved = True
                            # if can't access resolve_status by index, only None has been output
                            # e.g. fail, no matching line found in global scan file
                            # then try to align the mutation cluster using clustal.
                            # this step is rarely used and provides a minor improvement in resolution
                            # for the associated upset plot vs. not performing this step
                            except TypeError:
                                resolve_status = align_cluster_resolve(sequence_dic, final_direction, min(clus_start), max(clus_end))
                                if resolve_status[0]:
                                    aligned_cluster_directions.append(resolve_status)
                                    inferred_direction = resolve_status
                                    resolved_bool = "inferred (cluster aligned)"
                                    if inferred_direction in ["1234","1256"]:
                                        phylo_place = "ambiguous, species tree or ILS consistent"
                                    elif inferred_direction in ["3456"]:
                                        phylo_place = "ambiguous, species tree consistent"
                                    elif inferred_direction in ["16","24"]:
                                        phylo_place = "ILS consistent"
                                    else:
                                        phylo_place = "species tree consistent"
                                    cluster_align_resolved = True
                            # if everything fails, leave as unresolved
                            except KeyError as keyerr:
                                with open("key_errors.txt", "a") as f:
                                    f.write(str(keyerr))
                                inferred_direction = "unresolved"
                                phylo_place = "unresolved"
                        
                        # GET START/END CLUSTER POSITIONS IN HUMAN GENOME
                        if resolved_bool == "inferred (missing alignment block)":
                            set_annotation = "incomplete_detection"
                        elif resolved_bool == "inferred (failed filter)":
                            set_annotation = "partially_significant"
                        elif resolved_bool == "inferred (cluster aligned)":
                            set_annotation = "partially_significant"
                        else:
                            set_annotation = "unresolved"
                        try:
                            anc_desc = anc_desc_dic[inferred_direction]
                            if "homo" in anc_desc[0]:
                                # human_desc_bool = 1
                                anc_or_desc = "descendant"
                            elif "homo" in anc_desc[1]:
                                # human_anc_bool = 1
                                anc_or_desc = "ancestral"
                        except KeyError:
                            anc_or_desc = "uncertain"

                        fas4_fi = "fas4/" + fis.split("/")[1] + ".fas"
                        event_type = "."
                        fas4_homo_header = ""
                        with open(fas4_fi, "r") as f:
                            fas4_lines = f.readlines()
                            for i_line, line in enumerate(fas4_lines):
                                if line.startswith(">homo"):
                                    fas4_homo_header = line.strip()[1:]
                                    fas4_homo_seq = fas4_lines[i_line+1].strip()
                                    break
                        if fas4_homo_header != "":
                            human_chrom = fas4_homo_header.split(".")[1]
                            human_block_start = int(fas4_homo_header.split(".")[2])
                            if fas4_homo_header.split("_")[-1] == "neg":
                                strand = "reverse"
                                start_gap_count = fas4_homo_seq[:max(clus_end)].count("-")
                                end_gap_count = fas4_homo_seq[:(min(clus_start)+2)].count("-")
                                chrom_len = chrom_len_dic[human_chrom]
                                bed_start = chrom_len - human_block_start - max(clus_end)+1 + start_gap_count
                                bed_end = chrom_len - human_block_start - min(clus_start)-1 + end_gap_count
                            else:
                                strand = "forward"
                                start_gap_count = fas4_homo_seq[:(min(clus_start)+2)].count("-")
                                end_gap_count = fas4_homo_seq[:max(clus_end)].count("-")
                                bed_start = human_block_start + min(clus_start)+2 - start_gap_count
                                bed_end = human_block_start + max(clus_end) - end_gap_count
                                
                            bed_line = "\t".join(
                                [str(i) for i in [human_chrom,
                                 bed_start,
                                 bed_end,
                                 total_event_count,
                                 strand,
                                 set_annotation,
                                 anc_or_desc,
                                 event_type + "\n"]]
                            )
                            with open("hominid_template_switch_output/bed_output/template_switch_events.bed", "a") as bed_f:
                                bed_f.write(bed_line)
                        else:
                            bed_line = "\t".join(
                                [str(i) for i in [".",
                                 ".",
                                 ".",
                                 total_event_count,
                                 ".",
                                 set_annotation,
                                 ".",
                                 "." + "\n"]]
                            )
                            with open("hominid_template_switch_output/bed_output/template_switch_events.bed", "a") as bed_f:
                                bed_f.write(bed_line)

                        # PRINT INFERRED/UNRESOLVED EVENTS
                        event_title = "Event {0}".format(total_event_count)
                        if total_event_count != 1:
                            print("\n\n")
                        print("=" * len(event_title))
                        print(event_title)
                        print("=" * len(event_title))
                        print("Phylogenetic placement: {0}".format(phylo_place))
                        print("Significant comparisons: {0}".format(final_direction))
                        if inferred_direction:
                            print("Detected in comparisons: {0}".format(inferred_direction))
                        else:
                            print("Detected in comparisons: {0}".format(final_direction))
                        print("\nMultiple sequence alignment:\n")
                        print("    Human:", print_align_seq(sequence_dic["homo_sapiens"], min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["panhomo"], "down") + switch_point_annotations(comp_switch_point_dic["panhomo"], alignment_len, min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["homopan"], "up") + switch_point_annotations(comp_switch_point_dic["homopan"], alignment_len, min_coord, max_coord))
                        print("    Chimp:", print_align_seq(sequence_dic["pan_troglodytes"], min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["gorpan"], "down") + switch_point_annotations(comp_switch_point_dic["gorpan"], alignment_len, min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["pangor"], "up") + switch_point_annotations(comp_switch_point_dic["pangor"], alignment_len, min_coord, max_coord))
                        print("  Gorilla:", print_align_seq(sequence_dic["gorilla_gorilla"], min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["homogor"], "down") + switch_point_annotations(comp_switch_point_dic["homogor"], alignment_len, min_coord, max_coord))
                        print(print_arrow(comp_switch_point_dic["gorhomo"], "up") + switch_point_annotations(comp_switch_point_dic["gorhomo"], alignment_len, min_coord, max_coord))
                        print("    Human:", print_align_seq(sequence_dic["homo_sapiens"], min_coord, max_coord))
                        print("\n  Cluster:" + cluster_annotation(min(clus_start), max(clus_end), alignment_len, min_coord, max_coord))
                        print("\n")
                        print("Orangutan:", print_align_seq(sequence_dic["pongo_abelii"], min_coord, max_coord))
                        for comp_k in comp_switch_point_dic.keys():
                            if len(comp_switch_point_dic[comp_k]) != 0:
                                try:
                                    diverg = dd.pairwise_divergence(unique_events_pair_fastas[k][comp_k])

                                    if large_flank == 1:
                                        scan_flank_par = "100"
                                    else:
                                        scan_flank_par = "40"
                                    pairHMM_out = subprocess.check_output([pairHMM, "--pair", unique_events_pair_fastas[k][comp_k],
                                                                           "--print-file", unique_events_local_finames[k][comp_k],
                                                                           "--divergence", diverg,
                                                                           "--scan-flank", scan_flank_par])

                                    title = comp_translate[comp_k]
                                    print("\n\n{0}".format(title))
                                    print("-" * len(title))
                                    print(pairHMM_out.decode().lstrip())
                                except subprocess.CalledProcessError:
                                    pass
                        if resolved_bool == "unresolved":
                            unresolved_directions.append(final_direction)
                        else:
                            final_direction_list.append(inferred_direction)
                        if filter_resolved:
                            with open("hominid_template_switch_output/pairHMM_output/filter_resolved_pairHMM_output.csv", "a+") as f:
                                f.write(str(total_event_count) + "," + final_direction + "," + csv_read)
                        elif cluster_align_resolved:
                            with open("hominid_template_switch_output/pairHMM_output/clusteralign_resolved_pairHMM_output.csv", "a+") as f:
                                f.write(str(total_event_count) + "," + final_direction + "," + csv_read)
                        elif block_missing_resolved:
                            with open("hominid_template_switch_output/pairHMM_output/missing_resolved_pairHMM_output.csv", "a+") as f:
                                f.write(str(total_event_count) + "," + final_direction + "," + csv_read)
                        else:
                            with open("hominid_template_switch_output/pairHMM_output/unresolved_pairHMM_output.csv", "a+") as f:
                                f.write(str(total_event_count) + "," + final_direction + "," + csv_read)

    # write final direction output as csvs, used for plotting later
    # this information is also included in various other outputs in a more
    # human readable format
    with open("hominid_template_switch_output/directions/final_directions.csv", "w") as f:
        f.write(",".join(final_direction_list))
    with open("hominid_template_switch_output/directions/congruent_directions.csv", "w") as f:
        f.write(",".join(congruent_directions))
    with open("hominid_template_switch_output/directions/filter_fail_directions.csv", "w") as f:
        f.write(",".join(filter_fail_directions))
    with open("hominid_template_switch_output/directions/cluster_aligned_directions.csv", "w") as f:
        f.write(",".join(aligned_cluster_directions))
    with open("hominid_template_switch_output/directions/unresolved_directions.csv", "w") as f:
        f.write(",".join(unresolved_directions))
    with open("hominid_template_switch_output/directions/missing_directions.csv", "w") as f:
        f.write(",".join(missing_directions))

    # csv of fail reasons
    fail_df.to_csv("hominid_template_switch_output/fail_dataframe/event_fail_dataframe.csv", sep=",")


if __name__ == "__main__":
    main()
