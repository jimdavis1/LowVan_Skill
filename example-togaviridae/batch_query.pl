#! /usr/bin/env perl
# Batched BV-BRC SOLR query: same contract as query_PATRIC_bob.pl but sends
# identifiers N at a time via an "in" clause instead of one query per line.
use strict;
use P3DataAPI;
use Getopt::Long;
$| = 1;

my ($core, $input_type, $return, $batch) = ('genome_feature', 'genome_id', '', 200);
GetOptions('c=s' => \$core, 'i=s' => \$input_type, 'r=s' => \$return, 'b=i' => \$batch);
$return =~ s/ /,/g;
my @return = split(",", $return);

my $api = P3DataAPI->new();
my @buf;
while (<>) { chomp; next unless /\S/; push @buf, $_; flush_buf() if @buf >= $batch; }
flush_buf();

sub flush_buf {
    return unless @buf;
    my $list = join(",", map { my $v = $_; $v =~ s/([,()"\\])/\\$1/g; $v } @buf);
    $api->query_cb($core, \&handle_data, ["in", $input_type, "($list)"], ["select", $return]);
    print STDERR scalar(@buf), " sent\n";
    @buf = ();
}

sub handle_data {
    my ($data) = @_;
    for my $ent (@$data) { print join("\t", map { defined $ent->{$_} ? $ent->{$_} : '' } @return), "\n"; }
    return 1;
}
