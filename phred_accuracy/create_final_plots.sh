cd output/plots

convert AdapterRemoval_combined.png bbmerge_combined.png +append 1.png
convert ClipAndMerge_combined.png fastp_combined.png +append 2.png
convert leeHom_combined.png SeqPrep_combined.png +append 3.png
convert 1.png 2.png 3.png seqtk_adna_trim_combined.png -append final/phred_accuracy_combined.png
rm 1.png 2.png 3.png  

convert AdapterRemoval_0.png AdapterRemoval_-5.png AdapterRemoval_-10.png AdapterRemoval_-15.png AdapterRemoval_-20.png -append final/phred_accuracy_AdapterRemoval.png
convert bbmerge_0.png bbmerge_-5.png bbmerge_-10.png bbmerge_-15.png bbmerge_-20.png -append final/phred_accuracy_bbmerge.png
convert ClipAndMerge_0.png ClipAndMerge_-5.png ClipAndMerge_-10.png ClipAndMerge_-15.png ClipAndMerge_-20.png -append final/phred_accuracy_ClipAndMerge.png
convert fastp_0.png fastp_-5.png fastp_-10.png fastp_-15.png fastp_-20.png -append final/phred_accuracy_fastp.png
convert leeHom_0.png leeHom_-5.png leeHom_-10.png leeHom_-15.png leeHom_-20.png -append final/phred_accuracy_leeHom.png
convert SeqPrep_0.png SeqPrep_-5.png SeqPrep_-10.png SeqPrep_-15.png SeqPrep_-20.png -append final/phred_accuracy_SeqPrep.png
convert seqtk_adna_trim_0.png seqtk_adna_trim_-5.png seqtk_adna_trim_-10.png seqtk_adna_trim_-15.png seqtk_adna_trim_-20.png -append final/phred_accuracy_seqtk_adna_trim.png
