#!/usr/bin/perl

# blast_mask.pl takes
# -query or -i query file
# -db or -d database/target file
# and converts the -q or -d file to -ci (case_input) if specified
# uses blast tabular output -m8 or -outfmt 6 and gives back masked fasta
# with bases as lower/upper case (specified with -co or -case_output)
# and a table to stderr with counts of masked, unmasked bases, and (optionally) a redundancy score
# -evalue filters by evalue (picks only those lines with evalues <= -evalue)
# -perc_identity filters by percentage identity (only those with -p > perc_identity)
# -length - only those with length >= -length
# -score - only those with bitscore >= -score
# -max_target_seqs n (emulates blast+ behaviour)
# -hsp_per_target n (-h 1 emulates clobb behaviour)
# -redundancy - reports a redundancy score for each sequence in -i xor -d to STDERR (1 implies a single sequence hitting fully)

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($dbfile, $queryfile, $evalue, $perc_identity, $length, $bitscore, $opt_r, $max_target_seqs, $hsp_per_target, $Length_proportion);
my ($case_output, $case_input) = ("L", "N");
GetOptions (
  "db:s" => \$dbfile,
  "fasta|query|input:s" => \$queryfile,
  "case_output|co:s" => \$case_output,
  "case_input|ci:s" => \$case_input,
  "P|perc_identity:f" => \$perc_identity,
  "evalue:f" => \$evalue,
  "length:i" => \$length,
  "Length_proportion:f" => \$Length_proportion,
  "score:f" => \$bitscore,
  "redundancy" => \$opt_r,
  "max_target_seqs:i" => \$max_target_seqs,
  "hsp_per_target:i" => \$hsp_per_target,
);
$case_input  = uc(substr($case_input,0,1));
$case_output = uc(substr($case_output,0,1));
$case_output = "U" if $case_input eq "L";

die "Usage: blast_mask.pl [[-i|-query] <query_fasta_file> | [-d|-db] <target_fasta_file>] <blast_results_m8_format> [-co {U|L}] [-ci {U|L}]\nOnly one of query or target file should be specified, not both\n"
	unless ($dbfile xor $queryfile);

&multilinefastafile2singleline($queryfile) if $queryfile;
my $qsequences = &fastafile2hash($queryfile, $case_input) if $queryfile;
&multilinefastafile2singleline($dbfile) if $dbfile;
my $dsequences = &fastafile2hash($dbfile, $case_input) if $dbfile;

my (%qhitcount, %dhitcount); #hash (qid/did) of hash (position) with counts of hits at that position

my ($current_target_seqs, $current_hsps, $qid, $did) = (0,0,"","");

while (<>)
{
	next if /^#/;
	my @blastm8 = split /\t/;
	next if scalar @blastm8 < 12;
	next if $evalue and eval($blastm8[10]) > $evalue;
	next if $length and $blastm8[3] < $length;
	next if $perc_identity and $blastm8[2] < $perc_identity;
	next if $bitscore and $blastm8[11] < $bitscore;
	$current_target_seqs = 1 if $qid ne $blastm8[0];
	$current_target_seqs++ if $qid eq $blastm8[0] and $did ne $blastm8[1];
	$current_hsps = 1 if "$qid $did" ne "$blastm8[0] $blastm8[1]";
	$current_hsps++ if "$qid $did" eq "$blastm8[0] $blastm8[1]";
	$qid = $blastm8[0];
	$did = $blastm8[1];
	next if $max_target_seqs and $current_target_seqs > $max_target_seqs;
	next if $hsp_per_target  and $current_hsps > $hsp_per_target;

	my ($qoffset, $qlength, $doffset, $dlength, $qst, $qen, $dst, $den);
	if ($blastm8[6] < $blastm8[7]) { $qst = $blastm8[6]; $qen = $blastm8[7] } else { $qst = $blastm8[7] ; $qen = $blastm8[6] }
	if ($blastm8[8] < $blastm8[9]) { $dst = $blastm8[8]; $den = $blastm8[9] } else { $dst = $blastm8[9] ; $den = $blastm8[8] }
	if ($qid =~ /^(.+?)_(\d+)_(\d+)$/) { $qst = $qst + $2 - 1; $qen = $qen + $2 - 1; $qid = $1 }
	if ($did =~ /^(.+?)_(\d+)_(\d+)$/) { $dst = $dst + $2 - 1; $den = $den + $2 - 1; $did = $1 }
	$qoffset = ($qst - 1); $qlength = ($qen - $qst + 1);
	$doffset = ($dst - 1); $dlength = ($den - $dst + 1);

	if ($queryfile)
	{
		next if not defined $$qsequences{$qid};
		next if $Length_proportion and $qlength/length($$qsequences{$qid}{seq}) <= $Length_proportion;
		if ($opt_r) { for ($qst..$qen) { $qhitcount{$qid}{$_}++ } }
		if ($case_output eq "U")
		{ substr($$qsequences{$qid}{seq}, $qoffset, $qlength, uc(substr($$qsequences{$qid}{seq}, $qoffset, $qlength)) ) }
		else
		{ substr($$qsequences{$qid}{seq}, $qoffset, $qlength, lc(substr($$qsequences{$qid}{seq}, $qoffset, $qlength)) ) }
	}
	if ($dbfile)
	{
		next if not defined $$dsequences{$did};
		next if $Length_proportion and $dlength/length($$dsequences{$did}{seq}) <= $Length_proportion;
		if ($opt_r) { for ($dst..$den) { $dhitcount{$did}{$_}++ } }
		if ($case_output eq "U")
		{ substr($$dsequences{$did}{seq}, $doffset, $dlength, uc(substr($$dsequences{$did}{seq}, $doffset, $dlength)) ) }
		else
		{ substr($$dsequences{$did}{seq}, $doffset, $dlength, lc(substr($$dsequences{$did}{seq}, $doffset, $dlength)) ) }
	}
}

if ($queryfile) {foreach (sort keys %{$qsequences}) {print ">$_$$qsequences{$_}{desc}\n$$qsequences{$_}{seq}\n"}}
if ($dbfile) {foreach (sort keys %{$dsequences}) {print ">$_$$dsequences{$_}{desc}\n$$dsequences{$_}{seq}\n"}}
if ($queryfile) { &print_prop_covered($qsequences, \%qhitcount)}
if ($dbfile) { &print_prop_covered($dsequences, \%dhitcount)}

###################################################################################################
sub multilinefastafile2singleline
{
	my $fastafile = shift @_;
	unless (`wc -l $fastafile | cut -f1 -d' '` == (`grep -c '>' $fastafile` * 2))
	{
		my $tmpfile = rand();
		open FA,  "<$fastafile" or die "Can't open $fastafile\n";
		open TMP, ">$tmpfile" or die "Can't open $tmpfile for multiline fasta to singleline\n";
		my $i = 0;
		while (<FA>) {
			$i++;
			if (/^>/) { print TMP "\n" unless $i == 1; print TMP }
			elsif (/^\s*$/) { next }
			elsif (/\d/) { chomp; print TMP; print TMP " " }
			else { chomp; print TMP }
		}
		print TMP "\n";
		close FA;
		close TMP;
		system ("mv $tmpfile $fastafile");
	}
}
###################################################################################################
sub fastafile2hash
{
	my $fastafile = shift @_;
	my $changecase = "N"; $changecase = shift @_ if @_;
	my %sequences;
	open FA, "<$fastafile" or die $!;
	my $seqid = "";
	while (<FA>)
	{
		if (/^>(\S+)(.*)/) {
		    $seqid = $1;
    		$sequences{$seqid}{desc} = $2;
		}
		elsif ($changecase eq "L") { chomp($sequences{$seqid}{seq} .= lc(<FA>)) }
		elsif ($changecase eq "U") { chomp($sequences{$seqid}{seq} .= uc(<FA>)) }
		elsif ($changecase eq "N") { chomp($sequences{$seqid}{seq} .= <FA>) }
	}
	return \%sequences;
}
###################################################################################################
sub revcomp
{
	my $str = shift @_;
	$str =~ tr/atgcATGC/tacgTACG/;
	return scalar reverse($str);
}
###################################################################################################
sub print_prop_covered
{
	my ($sequences, $hitcount) = @_;
	for my $seqid (sort keys %{$sequences})
	{
		my ($masked_count, $unmasked_count);
		if ($case_output eq "L")
		{
			$masked_count = $$sequences{$seqid}{seq} =~ tr/[a-z]/[a-z]/;
			$unmasked_count = $$sequences{$seqid}{seq} =~ tr/[A-Z]/[A-Z]/;
		} else
		{
			$masked_count = $$sequences{$seqid}{seq} =~ tr/[A-Z]/[A-Z]/;
			$unmasked_count = $$sequences{$seqid}{seq} =~ tr/[a-z]/[a-z]/;
		}
		print STDERR "$seqid\t$$sequences{$seqid}{desc}\t$masked_count\t$unmasked_count";
		if ($opt_r) 
		{
			my $hitcount_total = 0;
			my $redundancy_inverse = 0;
			for (1..($masked_count + $unmasked_count))
			{
				$hitcount_total += $$hitcount{$seqid}{$_} if $$hitcount{$seqid}{$_};
				$redundancy_inverse += exists $$hitcount{$seqid}{$_}?1/$$hitcount{$seqid}{$_}:0;
			}
			print STDERR "\t" . sqrt($redundancy_inverse * $masked_count/($masked_count+$unmasked_count)**2); #. "\t" . ($hitcount_total?$masked_count/$hitcount_total : 0);
		}
		print STDERR "\n";
	}
}
