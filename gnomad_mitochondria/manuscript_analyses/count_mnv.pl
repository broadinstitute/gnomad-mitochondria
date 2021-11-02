#!/usr/bin/env perl
use warnings;
#  The Broad Institute
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

Usage: $0 [-h][-v] < infile

Script to identify multi-nucleotide variants (MNVs);  Input is tab-delimited file (one line per variant per person) assumed to be sorted by SAMPLE_ID then by POS;  Output is set of unique variants (POS.REF.ALT) with the allele count of that variant in the dataset (AC), the allele count of that variant as part of a multi-nucleotide variant (AC_MNV) and the percentage (AC_MNV/AC).  

infile is tab-del file:
SAMPLE_ID       POS.REF.ALT     FILTER  HL      DP      POS

Output rows per unique POS.REF.ALT
POS.REF.ALT AC AC_MNV PERC_AC_MNV

Example: adjacent variants 5185.G.A (present 1 sample) and 5186.A.T (present in 80 samples) were observed together in 1 sample and thus two lines would be output:
POS.REF.ALT AC AC_MNV PERC_AC_MNV
5185.G.A 1 1 1
5186.A.T 80 1 1/80

Flags:
h:     print this message
v:     verbose
USAGE
}

########################################################################
# Command Line Args
########################################################################
my $help=0;
my $verbose;

GetOptions(     'h' => \$help,
		'v' => \$verbose,
	  );
if ($help or scalar(@ARGV)< 0) {Usage;exit;}

########################################################################
# GLOBALS
########################################################################

# how many times have we seen this variant ID
my %id2count;
# how many times have we seen this variant ID in an MNV with an adjacent variant in this sample
my %id2mnv;
# make sure we don't double count; mark samples and IDs where we have seen an MNV 
my %sample2id2count;

########################################################################
# Main
########################################################################

# STEP 1: build hash tables with AC and AC_MNV
my $line=<>;
# parse line: s=SAMPLE_ID, id=POS.REF.ALT, filter=FILTER, HL=heteroplasmy level aka VAF,dep=depth, pos=position in chrM
my ($s,$id,$filter,$hl,$dp,$pos);
my $lastline="";
my $lastpos="";
my $lasts="";
my $lastid="";
while($line=<>) {
  chomp($line);
  ($s,$id,$filter,$hl,$dp,$pos)=split("\t",$line,-1);
  $id2count{$id}++; # increment # times have we seen this variant ID
  if (($s eq $lasts) && (defined($pos)) && (defined($lastpos)) && ($lastpos +1 == $pos)) {
    $id2mnv{$lastid}++ unless ($sample2id2count{$lasts}->{$lastid});
    $id2mnv{$id}++;
    $sample2id2count{$lasts}->{$lastid}=1;
    $sample2id2count{$s}->{$id}=1;
  }
  $lastline=$line;
  $lasts=$s;
  $lastpos=$pos;
  $lastid=$id;
}

# STEP 2: output this table of AC and AC_MNV and PERC_AC_MNV
print "POS.REF.ALT\tAC\tAC_MNV\tFRACTION_AC_MNV\n";
foreach $id (keys %id2count) {
  $id2mnv{$id}=0 unless ($id2mnv{$id});
  printf("%s\t%d\t%d\t%0.2f\n", $id, $id2count{$id}, $id2mnv{$id},$id2mnv{$id}/$id2count{$id})
}
