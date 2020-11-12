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

require(DNAshapeR)
require(gridExtra)
require(ggplot2)
require(ggpubr)


# read sequence files, predict DNA MGW, helical twist, roll, propeller twist
desc_filename="sp1_sequences_pm500nt/gold_standard_event_seqs_desc.tsv"
desc_preds = getShape(desc_filename, shapeType = 'Default', parse = TRUE,
                      methylate = FALSE, methylatedPosFile = NULL)

anc_filename="sp1_sequences_pm500nt/gold_standard_event_seqs_anc.tsv"
anc_preds = getShape(anc_filename, shapeType = 'Default', parse = TRUE,
                     methylate = FALSE, methylatedPosFile = NULL)

random_filename="random_seqs/grch38_random_sequences_concat.fa"
random_preds = getShape(random_filename, shapeType = 'Default', parse = TRUE,
                        methylate = FALSE, methylatedPosFile = NULL)


# save predictions to csvs
write.table(desc_preds$MGW, file="desc_MGW.csv", sep=",",col.names=F,row.names=F)
write.table(desc_preds$HelT, file="desc_HelT.csv", sep=",",col.names=F,row.names=F)
write.table(desc_preds$Roll, file="desc_Roll.csv", sep=",",col.names=F,row.names=F)
write.table(desc_preds$ProT, file="desc_ProT.csv", sep=",",col.names=F,row.names=F)

write.table(anc_preds$MGW, file="anc_MGW.csv", sep=",",col.names=F,row.names=F)
write.table(anc_preds$HelT, file="anc_HelT.csv", sep=",",col.names=F,row.names=F)
write.table(anc_preds$Roll, file="anc_Roll.csv", sep=",",col.names=F,row.names=F)
write.table(anc_preds$ProT, file="anc_ProT.csv", sep=",",col.names=F,row.names=F)

write.table(random_preds$MGW, file="rand_MGW.csv", sep=",",col.names=F,row.names=F)
write.table(random_preds$HelT, file="rand_HelT.csv", sep=",",col.names=F,row.names=F)
write.table(random_preds$Roll, file="rand_Roll.csv", sep=",",col.names=F,row.names=F)
write.table(random_preds$ProT, file="rand_ProT.csv", sep=",",col.names=F,row.names=F)
