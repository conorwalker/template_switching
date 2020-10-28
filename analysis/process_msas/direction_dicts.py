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


comp_dic = {"homopan":  1,
            "panhomo":  2,
            "homogor":  3,
            "gorhomo":  4,
            "pangor":  5,
            "gorpan":  6}


comp_translate = {"homopan":  "Chimp > Human",
                  "panhomo":  "Human > Chimp",
                  "homogor":  "Gorilla > Human",
                  "gorhomo":  "Human > Gorilla",
                  "pangor":  "Gorilla > Chimp",
                  "gorpan":  "Chimp > Gorilla"}


acceptable = ["13",
              "46",
              "35",
              "25",
              "1256",
              "3456",
              "1234",
              "16",
              "24"]


chrom_len_dic = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415
}


def resolved_direction(direct):
    dir_dic = {"1":  "13",
               "2":  "25",
               "4":  "46",
               "6":  "46",
               "15":  "1256",
               "26":  "1256",
               "14":  "1234",
               "23":  "1234",
               "36":  "3456",
               "45":  "3456",
               "156":  "1256",
               "256":  "1256",
               "125":  "1256",
               "126":  "1256",
               "134":  "1234",
               "234":  "1234",
               "123":  "1234",
               "124":  "1234",
               "456":  "3456",
               "356":  "3456",
               "345":  "3456",
               "346":  "3456"}
    try:
        return dir_dic[direct]
    except KeyError:
        return "invalid_key"


def determine_interest_dirs(direct):
    dir_dic = {"1":  "homogor",
               "2":  "pangor",
               "3":  ["homopan", "pangor"],
               "4":  "gorpan",
               "5":  ["panhomo", "homogor"],
               "6":  "gorhomo",
               "16":  ["panhomo", "pangor"],
               "15":  ["panhomo", "gorpan"],
               "26":  ["homopan", "pangor"],
               "24":  ["homopan", "homogor"],
               "14":  ["panhomo", "homogor"],
               "23":  ["homopan", "gorhomo"],
               "36":  ["gorhomo", "pangor"],
               "45":  ["homogor", "gorpan"],
               "12":  [("homogor", "gorhomo"), ("pangor", "gorpan")],
               "34":  [("homopan", "panhomo"), ("pangor", "gorpan")],
               "56":  [("homopan", "panhomo"), ("homogor", "gorhomo")],
               "156":  "panhomo",
               "256":  "homopan",
               "125":  "gorpan",
               "126":  "pangor",
               "134":  "panhomo",
               "234":  "homopan",
               "123":  "gorhomo",
               "124":  "homogor",
               "456":  "homogor",
               "356":  "gorhomo",
               "345":  "gorpan",
               "346":  "pangor"}
    return dir_dic[direct]


def determine_needed_code(direct):
    """
    Interprets the sequence code along multiple alignment length N, characters
    are V={1,3,5}. For the query direction, it specifies the required
    single characters required along the mutation cluster to resolve
    the direction.
    """
    dir_dic = {"1":  3,
               "2":  5,
               "3":  [1, 5],
               "4":  1,
               "5":  [1, 3],
               "6":  1,
               "15":  3,
               "23":  5,
               "36":  1,
               "45":  1,
               "26":  3,
               "14":  5,
               "12":  [3, 5],
               "34":  [5, 1],
               "56":  [3, 1],
               "156":  3,
               "256":  3,
               "125":  3,
               "126":  3,
               "134":  5,
               "234":  5,
               "123":  5,
               "124":  5,
               "456":  1,
               "356":  1,
               "345":  1,
               "346":  1}
    return dir_dic[direct]
