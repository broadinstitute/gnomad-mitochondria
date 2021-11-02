#!/usr/bin/env perl
use warnings;

use strict;
use Getopt::Long;

########################################################################
# Usage
########################################################################
sub Usage {

print <<USAGE;

Usage: $0 [-h][-v][-header][-nowarn][-byname][-x lookup_file=ref_col,id_col,merge_cols] ref.txt

Given a tab-del text file ref.txt, merge in data from lookup_file. This is similar in spirit to Excel vlookup function, or a database join query.

For each -x lookup_file=ref_col,id_col,merge_cols
 load in lookup_file and hash by id_col
 for each row in ref.txt, if ref_col = some val in id_col
 add cols to ref.txt for each col in merge_cols

NOTE: ref_col is column number of the key column in file ref.txt
      id_col is the column number of the key column in lookup_file
      merge_col is the column number of all other columns in lookup_file that we want to merge


EG. ref.txt has cols (probe,ensp) and human.ensp.enst.desc has cols (ensp,enst,sym,desc)
 merge_sets.pl -x human.ensp.enst.desc,1,0,1,2,3 ref.txt 
will produce file
 probe,ensp,enst,sym,desc
for all probes in ref.txt

Note, if the xref file has multiple rows with the same id, the first
row will be taken (and a warning will be given) by default.  There are special
characters to indicate other options when there are multiple rows with the 
same ID, eg use the row with the largest value (>), use the row with the 
greatest text value (g), use the row with the least text value (l).

If -header, then assume all files have headers, and merge these as well

Flags:
h:     print this message
v:     verbose
x:     file=ref_col,id_col,merge_cols
       NOTE: if the first of merge_cols start with > then use this key to
       sort (if duplicate matches, choose the one with the highest value)
       NOTE2: if the first of the merge_cols starts with {g | l} then
       choose duplicate with (g=greaterthan, l=lessthan) text value.
       NOTE3: last is replaced with the colnum of the last col in header row
       NOTE4: nlast is replaced with the colnum of the next-to-last col
       NOTE5: -x file=ref_col,id_col,merge_cols=default, then use default
       as the default value to use if no cross-reference is found.
       Usually default is null.
nowarn: do not give warnings of keys with multiple values
byname: if specified, then -x argument takes column NAMES not IDs. Note
        -header flag must be specified, and both key file and each -x file
        must have header line.  Note also that if byname, then cannot use
        the {g | l} version since this would get in the way of the headernames
USAGE
}

########################################################################
# Command Line Args
########################################################################
my $help=0;
my $verbose;
my @filecols;
my $nowarn="";
my $header="";
my $byname="";

GetOptions(     'h' => \$help,
		'v' => \$verbose,
		'x=s' => \@filecols,
		'nowarn' => \$nowarn,
		'header' => \$header,
		'byname' => \$byname,
	  );
if ($help or scalar(@ARGV)< 1) {Usage;exit;}
if ($byname && !$header) {
  warn "-header flag must be set when byname flag specified!\n"; exit; }
my $inf=shift;

########################################################################
# GLOBALS
########################################################################

my %key_colname2id;
my $next_key_colid=0;

########################################################################
# Main
########################################################################
my ($err,$ra_refrows,$ra_headers)=parse_rows($inf,$header,'same_row_len');
die "$err\n" if $err;

if ($header and $byname) {
  # map column names to ids
  my $num_headers=scalar(@$ra_headers);
  for (my $i=0; $i < $num_headers; $i++) {
    $key_colname2id{$ra_headers->[$i]}=$i;
  }
  $next_key_colid=$num_headers; # need this because we want to add -x names to this global hash
}

my ($xfile,$colstr,$ref_col,$id_col,$ra,$ra_find_row,$ra_xheader,$default_val);
my @merge_cols;
my @nulls_to_merge;
my %id2row;
my $sort_by_col;
my $sort_by_type;
foreach my $xfilecols (@filecols) {
  ($xfile,$colstr,$default_val)=split("=",$xfilecols);
  $default_val="" unless (defined($default_val));

  # if colstr contains "last" for the reference col replace with the correct number
  if ($colstr =~ /^(n*)last/) { # last or nlast
    my $last_col=scalar(@$ra_headers)-1;
    my $nlast_col=$last_col-1;
    $colstr =~ s/^nlast/$nlast_col/g;
    $colstr =~ s/^last/$last_col/g;
  }
  # parse xfile
  ($err,$ra,$ra_xheader)=parse_rows($xfile,$header);
  die "$err\n" if $err;
  my %x_colname2id;
  if ($header and $byname) {
    # map column names to ids
    my $num_headers=scalar(@$ra_xheader);
    for (my $i=0; $i < $num_headers; $i++) {
      $x_colname2id{$ra_xheader->[$i]}=$i;
    }
  }

  # colstr contains last for the xfile
  if ($colstr =~ /(n*)last/) {
    my $last_col=scalar(@$ra_xheader)-1;
    my $nlast_col=$last_col-1;
    $colstr =~ s/nlast/$nlast_col/g;
    $colstr =~ s/last/$last_col/g;
  }
  ($ref_col,$id_col,@merge_cols)=split(",",$colstr,-1);

  # handle special case of greater than or less than in merge_cols
  if ((scalar(@merge_cols)) && 
      ((($merge_cols[0] =~ /^([><])(.*)$/) && ($byname)) ||  # if byname, cannot use gl
       (($merge_cols[0] =~ /^([><gl])(.*)$/) && !($byname))
      )) {
    $merge_cols[0]=$2;
    $sort_by_col=$merge_cols[0];
    $sort_by_type=$1;
  } else {
    $sort_by_col=undef;
  }

  if ($byname) {
    # convert col names to IDs
    die "Cannot find key column named $ref_col!\n" unless (defined($key_colname2id{$ref_col}));
    die "Cannot find xfile $xfile column named $id_col!\n" unless (defined($x_colname2id{$id_col}));
    $ref_col=$key_colname2id{$ref_col};
    $id_col=$x_colname2id{$id_col};
    # now convert each merge_col
    for (my $i=0; $i < scalar(@merge_cols); $i++) {
      my $x_col=$merge_cols[$i];
      die "Cannot find xfile $xfile column named $x_col!\n" unless (defined($x_colname2id{$x_col}));
      $merge_cols[$i]=$x_colname2id{$x_col};
      # now, after we do the merge, the key file should now have this x_col
      unless ($key_colname2id{$x_col}) {
	$key_colname2id{$x_col}=$next_key_colid;
	$next_key_colid++;
      }
    }
  }

  # initialize null data
  @nulls_to_merge=();
  foreach my $col (@merge_cols) {
    push @nulls_to_merge,$default_val;
  }

  # build hash keyed by id_col
  %id2row=();
  my ($orig_val,$new_val);
  foreach my $ra_row (@$ra) {
    if ($id2row{$ra_row->[$id_col]}) {
      # we have already seen this id before
      if (defined($sort_by_col)) {
	$orig_val=$id2row{$ra_row->[$id_col]}->[$sort_by_col];
	$new_val=$ra_row->[$sort_by_col];
	# if sort_by_col, then we choose to replace original saved row (orig) with new row
	# if new_val > orig_val

	next unless (defined($new_val) && ($new_val ne "")); # use orig unless new_val defined
	
	# use new if new_val > orig_val or orig_val not defined
	if (!defined($orig_val) || ($orig_val eq "")) {
	  $id2row{$ra_row->[$id_col]}=$ra_row;
	} elsif (($sort_by_type eq ">") &&  ($new_val > $orig_val)) {
	  $id2row{$ra_row->[$id_col]}=$ra_row;
	} elsif (($sort_by_type eq "<") &&  ($new_val < $orig_val)) {
	  $id2row{$ra_row->[$id_col]}=$ra_row;
	} elsif (($sort_by_type eq "g") &&  ($new_val gt $orig_val)) {
	  $id2row{$ra_row->[$id_col]}=$ra_row;
	} elsif (($sort_by_type eq "l") &&  ($new_val lt $orig_val)) {
	  $id2row{$ra_row->[$id_col]}=$ra_row;
	} 
	next;
      } 
      # otherwise we simply use the first row we came across
      # but give a warning and ignore this second entry
      unless ($nowarn) {
	my $str1=join(",",@$ra_row);
	my $str2=join(",",@{$id2row{$ra_row->[$id_col]}});
	next if ($str1 eq $str2); # same!!!
	printf(STDERR "%s: conflicting entry for id %s:\n new=%s\n old=%s\n",
	       $xfile, $ra_row->[$id_col], 
	       $str1,$str2);
	# otherwise, just use the first row that we came across
      }
    } else {
      $id2row{$ra_row->[$id_col]}=$ra_row
	unless (defined($ra_row->[$id_col]) && ($ra_row->[$id_col] eq ""));
    }
  }

  # merge with ref
  foreach my $ra_row (@$ra_refrows) {
    if($ra_row->[$ref_col]) {
      # find the row in id2row corresponding to this refrow
      $ra_find_row=$id2row{$ra_row->[$ref_col]};   
    } else {
      $ra_find_row=undef;
    }
    if ($ra_find_row) {
      # merge in items
      # do this in order "" || val so that if val is 0 then this is used
      # but otherwise no error message is written
      push @$ra_row,map {"" || $ra_find_row->[$_]} @merge_cols;
    } else {
      push @$ra_row,@nulls_to_merge
    }
  }
  # merge in headers
  if ($header) {
    push @$ra_headers,map {$ra_xheader->[$_]} @merge_cols;
  }
}

# print output
if ($header) {
  print join("\t",@$ra_headers)."\n";
}
my $tab=""; 
foreach my $ra_row (@$ra_refrows) {
  $tab=""; # first time through no tab
  foreach my $val (@$ra_row) {
    $val = "" unless defined($val);
    print "$tab$val";
    $tab="\t"; # after the first time, print the tab
  }
  print "\n";
}

########################################################################
# Subroutines
########################################################################

######################################################################

=head2 parse_rows

 ($err,$ra_rows,$header)=parse_rows($file,$header)

 given tab-del $file,

 return a reference to an array of rows, where each row
 is a ref to an array of columns

 if $header, return a ref to the header row separately

=cut

######################################################################
sub parse_rows {
  my ($file,$header,$same_row_len,$header_only)=@_;

  open INF, "<$file" or return ("Cannot open $file for reading\n");
  my @rows;
  my @header_row;
  my $num_fields;
  while(<INF>) {
    chomp;
    my @fields=split("\t",$_,-1); # keep trailing null columns
    $num_fields=scalar(@fields) unless (defined($num_fields));
    if ($header) {
      @header_row=@fields;
      $header="";
      next;
    }
    if (defined($same_row_len)) {
      while ($num_fields > scalar(@fields)) {
	push @fields, "";
      }
    }
    push @rows,\@fields;
  }
  close(INF);
  return (undef,\@rows,\@header_row);
}
