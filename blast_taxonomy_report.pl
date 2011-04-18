#!/usr/bin/perl

=head1 NAME

blast_taxonomy_report.pl - Takes a blast output (in tabular format) and the NCBI taxonomy nodes and names database text dumps and returns an NCBI style taxonomy report of the hits

=head1 SYNOPSIS

blast_taxonomy_report.pl -blast contigs_nr.1e-8.txt -nodes nodes.dmp -names names.dmp

=head1 DESCRIPTION

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2010.05.24

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($blastfile, $nodesfile, $namesfile, $gi_taxid_file) = ("nomito_auto_data_35_nr.1e-10.max5.gz","nodes.dmp.gz","names.dmp.gz", "gi_taxid_prot.dmp.gz");
my %tax_levels;

GetOptions (
    "blastfile=s"     => \$blastfile,
    "nodesfile=s"     => \$nodesfile,
    "namesfile=s"     => \$namesfile,
    "gi_taxid_file=s" => \$gi_taxid_file,
    "tax_levels=s"    => \%tax_levels,
);
%tax_levels = ("superkingdom",1,"phylum",1,"family",1,"species",1) unless keys %tax_levels >0;

die "Usage: blast_taxonomy_report.pl -blast contigs_nr.1e-8.txt -nodes nodes.dmp -names names.dmp -gi_taxid_file gi_taxid_prot.dmp.gz\n" unless -r $blastfile or $blastfile eq "-";

#---------------------------------

my (%taxid_has_parent, %taxid_has_taxlevel, %taxid_has_name, %gi_has_taxid);

&load_nodes_names;
&load_gi_taxids($blastfile, $gi_taxid_file);

#initialise curr query var and hits array
my $curr_query = "";
my @curr_query_hits;

# parse each line of blast tabular results file
my $blast_fh = &read_fh($blastfile);
while (my $line = <$blast_fh>) {
    my @blast_fields = split(/\t/,$line);

    # skip loop if it is not a blast tab result
    next if $line =~ /^#/ or scalar @blast_fields < 12;
    
    if ($blast_fields[0] eq $curr_query) {
        push @curr_query_hits, $blast_fields[1];
    }
    else {
        # process existing query and then start new query
        &print_taxonomy_report ($curr_query, @curr_query_hits);
        $curr_query = $blast_fields[0];
        @curr_query_hits = ($blast_fields[1]);
    }
}
# last one:
&print_taxonomy_report ($curr_query, @curr_query_hits);

close $blast_fh;

############################################################

sub print_taxonomy_report {
    my $query = shift @_;
    my @hits  = @_;
    my %hit_counts;
    print "$query\t";
    for my $hit_id (@hits) {
        next unless $hit_id =~ /^gi\|(\d+)\|/;
        my $hit_gi = $1;
        if (exists $gi_has_taxid{$hit_gi}) {
            my $hit_taxid = $gi_has_taxid{$hit_gi};
            my @parents = &get_parents($hit_taxid);
            # print "@parents\n";
            for my $parent (@parents) {
                if (exists $taxid_has_taxlevel{$parent} and exists $tax_levels{$taxid_has_taxlevel{$parent}}) {
                    $hit_counts{$taxid_has_taxlevel{$parent}}{$taxid_has_name{$parent}}++;
                }
            }
        }
    }
    for my $tax_level (keys %hit_counts) {
        for my $tax_name (keys %{$hit_counts{$tax_level}}) {
            print "$tax_level\t$tax_name\t$hit_counts{$tax_level}{$tax_name}\t";
        }
    }
    print "\n";
}

############################################################

sub get_parents {
    my @all = @_;
    my $current_id = $all[0];
    if (exists $taxid_has_parent{$current_id} and $current_id ne $taxid_has_parent{$current_id}) {
        unshift @all, $taxid_has_parent{$current_id};
        @all = &get_parents(@all);
    }
    return @all;
}

############################################################

sub load_nodes_names {
    my $fh;
    
    $fh = &read_fh($nodesfile);
    while (my $line = <$fh>) {
        # line in nodes.dmp should match the regexp below.
        # Change the regexp if NCBI changes their file format
        next if $line !~ /^(\d+)\s*\|\s*(\d+)\s*\|\s*(.+?)\s*\|/;
        $taxid_has_parent{$1} = $2;
        $taxid_has_taxlevel{$1} = $3;
    }
    close $fh;
    
    $fh = &read_fh($namesfile);
    while (my $line = <$fh>) {
        next unless $line =~ /^(\d+)\s*\|\s*(.+?)\s*\|.+scientific name/;
        $taxid_has_name{$1} = $2;
    }
}

############################################################

sub load_gi_taxids {

    # need blast file so that only those gis are loaded which are hits
    my $blastfile = shift @_;
    my $gi_taxid_file = shift @_;

    # initialize hash of gis of hits (to be filled in later with taxids)
    my $fh = &read_fh($blastfile);
    while (my $line = <$fh>) {
        next unless $line =~ /^\S+\tgi\|(\d+)\|/;
        $gi_has_taxid{$1} = 0;
    }
    close $fh;

    # now load gi_taxid mappings
    $fh = &read_fh($gi_taxid_file);
    while (my $line = <$fh>) {
        # Change the regexp if NCBI changes their file format for gi_taxid.dmp
        next if $line !~ /^(\d+)\s+(\d+)/;
        if (exists $gi_has_taxid{$1}) {
            $gi_has_taxid{$1} = $2;
        }
    }
    close $fh;
    
}

############################################################

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

############################################################

# use Proc::ProcessTable;
# 
# sub memory_usage {
#     my $t = new Proc::ProcessTable;
#     foreach my $got ( @{$t->table} ) {
#         next if not $got->pid eq $$;
#         return $got->size;
#     }
# }
# 
# print STDERR "\nmemory: ". &memory_usage()/1024/1024 ."\n";
# exit;
