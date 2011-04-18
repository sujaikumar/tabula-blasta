#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($queryfile, $dbfile) = ("","");
GetOptions (
  "i|query:s" => \$queryfile,
  "db:s" => \$dbfile,
);

my $query_seqs = &fastafile2hash($queryfile);
my $db_seqs = &fastafile2hash($dbfile);

while (<>) {
    chomp;
    next if /^#/;
    my @fields = split /\t/;
    my ($qid, $did, $perc_identity, $length, $gaps, $mismatches, $qst, $qen, $dst, $den, $evalue, $score) = @fields;
    my $qlen = length($$query_seqs{$qid}{seq});
    my $dlen = length(   $$db_seqs{$did}{seq});
    my $query_strand = "Plus";
    my $db_strand = "Plus";
    # no need for next 4 lines because qst always < qen in blast    
    #     if ($qst > $qen) {
    #         ($qst, $qen) = ($qen, $qst);
    #         $query_strand = "Minus";
    #     }
    if ($dst > $den) {
        ($dst, $den) = ($dlen - $dst + 1, $dlen - $den + 1);
        $db_strand = "Minus";
    }
    print "$qid\t$qlen\t$qst..$qen\t$did\t$dlen\t$dst..$den\t$evalue\t"
        . sprintf("%.2f", $perc_identity/100) . "\t$score\t$query_strand / $db_strand\n"; 
}

#############################################################################

sub fastafile2hash
{
    my $fastafile = shift @_;
    my %sequences;
    my $fh = &read_fh($fastafile);
    my $seqid;
    while (<$fh>)
    {
        if (/^>(\S+)(.*)/) {
            $seqid = $1;
            $sequences{$seqid}{desc} = $2;
            $sequences{$seqid}{order} = rand;
        }
        else {
            chomp($sequences{$seqid}{seq} .= $_)
        }
    }
    return \%sequences;
}

#############################################################################

sub read_fh {
    my $filename = shift @_;
    my $filehandle;
    if ($filename =~ /gz$/) {
        open $filehandle, "gunzip -dc $filename |" or die $!;
    }
    else {
        open $filehandle, "<$filename" or die $!;
    }
    return $filehandle;
}
