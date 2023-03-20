cd output/plots

convert edit_distances_AdapterRemoval.png edit_distances_bbmerge.png +append 1.png
convert edit_distances_ClipAndMerge.png edit_distances_fastp.png +append 2.png
convert edit_distances_leeHom.png edit_distances_SeqPrep.png +append 3.png
convert 1.png 2.png 3.png edit_distances_seqtk_adna_trim.png -append final/merging_insert_lengths.png

rm 1.png 2.png 3.png  
