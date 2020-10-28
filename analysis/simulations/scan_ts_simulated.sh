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


tsa_path="";

for f in simulated_sequences/*; do
    switch=$(head -n 1 $f | cut -d"_" -f8);
    fpa_out=$(tsa_path --scan --pair $f --divergence $1 | grep -v "scan" | grep -v "chrom");
    for fpa_line in $fpa_out; do
        echo "$f","$switch","$fpa_line";
    done;
done > ts_simulated_scanned.csv
