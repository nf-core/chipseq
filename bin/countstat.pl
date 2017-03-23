#!/usr/bin/env perl
use strict;
use warnings;

open(OUTPUT, ">read_count_statistics.txt");
my @fileList = @ARGV;

print OUTPUT "File\tTotalCounts\tUniqueCounts\tUniqueStartCounts\tUniqueRatio\tUniqueStartRatio\n";

foreach my $f(@fileList){

 open(IN,"<$f")||die $!;
 my $Tcnt=0;
 my $prev="NA";
 my $lcnt=0;
 my $Tcnt_2=0;
 my $prev_2="NA";

 while(<IN>){
   chomp;
   my @line=split("\t",$_);
   $lcnt++;
   my $t = join("_",@line[0..2]);
   $Tcnt++ unless($t eq $prev);
   $prev=$t;

   my $t_2 = join("_",@line[0..1]);
   $Tcnt_2++ unless($t_2 eq $prev_2);
   $prev_2=$t_2;
 }

 print OUTPUT "$f\t$lcnt\t$Tcnt\t$Tcnt_2\t".($Tcnt/$lcnt)."\t".($Tcnt_2/$lcnt)."\n";
 close(IN);
}
close(OUTPUT);