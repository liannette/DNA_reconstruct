#!/usr/bin/perl


use strict;
use warnings;

#>4:+:155073817:155073849:32
#ACCTCCCAAGCTCAAGTGATCCTCCCACCTCA
my $seqf="ACCTCCCAAGCTCAAGTGATCCTCCCACCTC";
my $seqr="GAGGTGGGAGGATCACTTGAGCTTGGGAGGT";

my $adpt1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
my $adpt2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";

my $filef="qscores_s1.fq";
my $filer="qscores_s2.fq";

open(FILEF,">".$filef) or die "cannot open $filef";
open(FILER,">".$filer) or die "cannot open $filer";


#identical
for(my $i=0;$i<=41;$i+=1){
    for(my $j=0;$j<=41;$j+=1){
	#def line
	print FILEF "\@i_G".$i."_G".$j."/1\n";
	print FILER "\@i_G".$i."_G".$j."/2\n";
	#seq
	print FILEF substr($seqf,0,15)."G".substr($seqf,15+1,15)."".substr($adpt1,0,44)."\n";
	print FILER substr($seqr,0,15)."C".substr($seqr,15+1,15)."".substr($adpt2,0,44)."\n";
	#+
	print FILEF "+\n";
	print FILER "+\n";
	#QS
	print FILEF (chr(41+33)x15).chr($i+33).(chr(41+33)x15).(chr(41+33)x44)."\n";
	print FILER (chr(41+33)x15).chr($j+33).(chr(41+33)x15).(chr(41+33)x44)."\n";

	warn $i."\t".$j."\n";
	
    }
}

#identical
for(my $i=0;$i<=41;$i+=1){
    for(my $j=0;$j<=41;$j+=1){
	#def line
	print FILEF "\@d_G".$i."_T".$j."/1\n";
	print FILER "\@d_G".$i."_T".$j."/2\n";
	#seq
	print FILEF substr($seqf,0,15)."G".substr($seqf,15+1,15).substr($adpt1,0,44)."\n";
	print FILER substr($seqr,0,15)."A".substr($seqr,15+1,15).substr($adpt2,0,44)."\n";
	#+
	print FILEF "+\n";
	print FILER "+\n";
	#QS
	print FILEF (chr(41+33)x15).chr($i+33).(chr(41+33)x15).(chr(41+33)x44)."\n";
	print FILER (chr(41+33)x15).chr($j+33).(chr(41+33)x15).(chr(41+33)x44)."\n";

	warn $i."\t".$j."\n";
	
    }
}


close(FILEF);
close(FILER);


#open(FILE,$ARGV[0]) or die "cannot open ".$ARGV[0];
#while(my $line = <FILE>){
#  chomp($line);
#  print $line;
#}
#close(FILE);

