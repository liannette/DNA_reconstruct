cd output/plots

convert fraglen_distributions_AdapterRemoval.png fraglen_distributions_bbmerge.png +append 1.png
convert fraglen_distributions_ClipAndMerge.png fraglen_distributions_fastp.png +append 2.png
convert fraglen_distributions_leeHom.png fraglen_distributions_SeqPrep.png +append 3.png
convert 1.png 2.png 3.png fraglen_distributions_seqtk_adna_trim.png -append final/merging_insert_length_distributions.png

rm 1.png 2.png 3.png  
