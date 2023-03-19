cd output/plots

convert AdapterRemoval_match.png bbmerge_match.png +append match1.png
convert ClipAndMerge_match.png fastp_match.png +append match2.png
convert leeHom_match.png SeqPrep_match.png +append match3.png
convert match1.png match2.png  -append final/match_1.png
convert match3.png seqtk_adna_trim_match.png -append final/match_2.png
      
convert AdapterRemoval_mismatch.png bbmerge_mismatch.png +append mismatch1.png
convert ClipAndMerge_mismatch.png fastp_mismatch.png +append mismatch2.png
convert leeHom_mismatch.png SeqPrep_mismatch.png +append mismatch3.png
convert mismatch1.png mismatch2.png  -append final/mismatch_1.png
convert mismatch3.png seqtk_adna_trim_mismatch.png -append final/mismatch_2.png

rm match1.png match2.png match3.png  
rm mismatch1.png mismatch2.png mismatch3.png  