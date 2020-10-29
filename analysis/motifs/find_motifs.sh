# Copyright (C) 2020 EMBL - European Bioinformatics Institute
# Contact: goldman@ebi.ac.uk, cwalker@ebi.ac.uk

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.

# You should have received a copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.

meme gold_standard_ancestral_sequences_point1_pm150.fa -dna -mod zoops -minw 6 -maxw 10 -p 4 -objfun de -neg ape_sequences_301bp_samples.fa -revcomp -markov_order 0 -seed 42 -nmotifs 50 -o zoops_differential_enrichment_ancestral_pm150_min6_max10
meme gold_standard_ancestral_sequences_point1_pm150.fa -dna -mod zoops -minw 10 -maxw 20 -p 4 -objfun de -neg ape_sequences_301bp_samples.fa -revcomp -markov_order 0 -seed 42 -nmotifs 50 -o zoops_differential_enrichment_ancestral_pm150_min10_max20
meme gold_standard_ancestral_sequences_point1_pm150.fa -dna -mod zoops -minw 20 -maxw 50 -p 4 -objfun de -neg ape_sequences_301bp_samples.fa -revcomp -markov_order 0 -seed 42 -nmotifs 50 -o zoops_differential_enrichment_ancestral_pm150_min20_max50