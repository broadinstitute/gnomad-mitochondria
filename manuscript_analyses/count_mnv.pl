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

Usage: $0 [-h][-v] < infile

infile is tab-del file:
SAMPLE_ID       POS.REF.ALT     FILTER  HL      DP      POS

with rows sorted by sample and then POS

Output rows per unique POS.REF.ALT
POS.REF.ALT AC AC_MNV PERC_AC_MNV

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
#make sure we don't double count; mark samples and IDs where we have seen an MNV 
my %sample2id2count;

########################################################################
# Main
########################################################################

# STEP 1: build hash tables with AC and AC_MNV
my $line=<>;
my ($s,$id,$filter,$hl,$dp,$pos);
my $lastline="";
my $lastpos="";
my $lasts="";
my $lastid="";
while($line=<>) {
  chomp($line);
  ($s,$id,$filter,$hl,$dp,$pos)=split("\t",$line,-1);
  $id2count{$id}++;
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


########################################################################
# Subroutines
########################################################################
