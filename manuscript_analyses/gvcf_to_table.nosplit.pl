#!/usr/bin/env perl
use warnings;
#  The Broad Institute
#  SOFTWARE COPYRIGHT NOTICE AGREEMENT
#  This software and its documentation are copyright (YEAR) by the
#  Broad Institute/Massachusetts Institute of Technology. All rights are
#  reserved.
#
#  This software is supplied without any warranty or guaranteed support
#  whatsoever. Neither the Broad Institute nor MIT can be responsible for its
#  use, misuse, or functionality.

use strict;
use Data::Dumper;
use Getopt::Long;

my $DEBUG=1;

########################################################################
# Usage
########################################################################
sub Usage {

print <<USAGE;

Usage: $0 [-h][-v][-m min_vaf][-pass i][-vafcol i] < infile

Input GVCF file eg

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  GTEX-1117F-0003-SM-6WBT7        GTEX-111CU-0003-SM-6WBUD        GTEX-111FC-0001-SM-6WBTJGTEX-111VG-0004-SM-6WBTS GTEX-111YS-0004-SM-6WBTN        GTEX-1122O-0004-SM-6WBTE
MT      10      .       T       C       .       .       AF=.    GT:DP:HL:MQ:TLOD:FT     ./.:.:.:.       ./.:.:.:.       ./.:.:.:.       ./.:.:.:.       ./.:.:.:.       ./.:.:.:.
MT      16      .       A       T       .       .       AF=.    GT:DP:VL:VQ     ./.:.:.:.       ./.:.:.:.       ./.:.:.:.       ./.:.:.:.       ./.:.:.:.       ./.:.:.:.

Output file (one line per sample-allele pair)

SAMPLE_ID POS.REF.ALT FILTER HL DP

Flags:
h:     print this message
v:     verbose
m:     minimum VAF to output
pass:  if set, this is the 0-based column that should be PASS, otherwise no variants are output
USAGE
}

########################################################################
# Command Line Args
########################################################################
my $help=0;
my $verbose;
my $min_vaf=0;
my $passcol=undef;
my $coldp=1;
my $colvaf=2;
my $colfilter=5;
GetOptions(     'h' => \$help,
		'v' => \$verbose,
		'm=f' => \$min_vaf,
		'pass=i' => \$passcol,
		'vafcol=i' => \$colvaf,
		'filtercol=i' => \$colfilter,
		'dpcol=i' => \$coldp,
	  );
if ($help or scalar(@ARGV)< 0) {Usage;exit;}

########################################################################
# GLOBALS
########################################################################

my %name2col;
my %col2name;

# 0-based columns within genotype field
#my $coldp=1;
#my $colvaf=2;
#my $colfilter=5;

# 0-based column per line
my $colpos=1;
my $colref=3;
my $colalt=4;

########################################################################
# Main
########################################################################

my $line;
my @fields;
my @headers;
my $num_fields;
my $num_samples;
my $genotype_start_col;
my ($str,$allele);
my @genos;
#print "SAMPLE_ID\tSAMPLE_ID:ALLELE\tALLELE\tCHROM\tPOS\tREF\tALT\tGT\tDP\tVL\tVQ\tVAF_ALT1\tVAF_ALT2\tVAF_ALT3\n";
print "SAMPLE_ID\tPOS.REF.ALT\tFILTER\tHL\tDP\n";
while ($line=<>) {
  chomp($line);
  next if ($line =~ /^\#\#/);
  @fields=split("\t",$line);
  if ($line =~ /^#CHROM/) {
    @headers=@fields;
    $num_fields=scalar(@headers);
    for (my $i=0; $i<$num_fields; $i++) {
      $name2col{$headers[$i]}=$i;
      $col2name{$i}=$headers[$i];
      if ($headers[$i] eq "FORMAT") {
	$genotype_start_col=$i+1;
	$num_samples=$num_fields-$i-1;
      }
    }
    $name2col{"CHROM"}=$name2col{"\#CHROM"};
  } else {
    next if (defined($passcol) && ($fields[$passcol] ne "PASS"));
      
    $allele=sprintf("%s.%s.%s",$fields[$colpos],$fields[$colref],$fields[$colalt]);
    for (my $i=$genotype_start_col; $i<$num_fields; $i++) {
      $str=$fields[$i];
      #print STDERR "str=$str\n";
      if (($str ne "./.") && ($str ne ".")){
	if (substr($str,0,3) ne "./.") {
	  @genos=split(":",$str,-1);
	  if ($genos[$colvaf] >= $min_vaf) {
	    # output: SAMPLE_ID ALLELE FILTER VAF DP
	    printf("%s\t%s\t%s\t%s\t%s\n",
		   $col2name{$i}, $allele,$genos[$colfilter],$genos[$colvaf],$genos[$coldp]);
	  }
	}
      }
    } 
  }
}

#print STDERR Dumper(\%name2col);
#print STDERR "num_samples=$num_samples\n"

########################################################################
# Subroutines
########################################################################
