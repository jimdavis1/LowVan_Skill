package FCP_Trunc_Utils;

#===============================================================================
#  Truncation QC.
#
#  Runs AFTER every cluster has been aligned with MAFFT but BEFORE any
#  alignment is corrected or turned into a PSSM.
#
#  The problem it solves: mmseqs clusters on bidirectional coverage, so a
#  fragment cannot join a full-length cluster -- but a set of similarly
#  truncated sequences happily clusters with ITSELF, covering each other
#  100%. That cluster then produces a perfectly clean-looking alignment and
#  a PSSM for a protein that does not exist. Neither the -c coverage rule,
#  the -fd dash rule, nor a "does it start with M" eyeball catches it,
#  because a C-terminally truncated protein starts with a good M and can sit
#  inside a literature-seeded length bound.
#
#  Method: take the consensus of each aligned cluster, blastp them all
#  against each other, and drop any cluster whose consensus is essentially a
#  SUBSTRING of a longer surviving cluster's consensus -- i.e. nearly all of
#  the short one aligns, at high identity, inside something materially
#  longer.
#
#  Nothing is deleted. Rejected alignments are moved to truncated_alis/ and
#  every decision is written to Truncation_Report.
#===============================================================================

use strict;
use warnings;
use File::Copy qw(move);
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = qw(qc_truncated_clusters);

use gjoseqlib;

#-------------------------------------------------------------------------------
#  Consensus of an aligned cluster: the modal residue of every column whose
#  modal character is not a gap. Ties broken alphabetically so the consensus
#  is deterministic.
#-------------------------------------------------------------------------------
sub alignment_consensus
{
    my (@ali) = @_;
    return "" unless @ali;
    my $len = length($ali[0][2]);
    my $cons = "";

    for my $c (0 .. $len - 1)
    {
        my %cnt;
        for my $s (@ali) { $cnt{ substr($s->[2], $c, 1) }++ }
        my ($best) = sort { $cnt{$b} <=> $cnt{$a} or $a cmp $b } keys %cnt;
        $cons .= $best unless $best eq '-';
    }
    return $cons;
}

#-------------------------------------------------------------------------------
#  qc_truncated_clusters(\%aligned, $options)
#
#  %aligned : cluster_label => arrayref of gjoseqlib triples (the MAFFT ali)
#
#  Returns a hashref: { keep => [labels], drop => { label => reason } }
#-------------------------------------------------------------------------------
sub qc_truncated_clusters
{
    my ($aligned, $options) = @_;

    my $pident    = defined $options->{qc_pident}   ? $options->{qc_pident}   : 90;
    my $min_qcov  = defined $options->{qc_qcov}     ? $options->{qc_qcov}     : 0.95;
    my $len_ratio = defined $options->{qc_lenratio} ? $options->{qc_lenratio} : 0.92;

    my @labels = sort keys %$aligned;
    my %cons;
    for my $lab (@labels)
    {
        my $c = alignment_consensus(@{ $aligned->{$lab} });
        $cons{$lab} = $c if length $c;
    }

    open(my $rep, ">Truncation_Report")
        or die "Cannot open Truncation_Report for writing: $!\n";
    print $rep "Cluster\tN_Seqs\tConsensus_Len\tVerdict\tSubstring_Of\tPercent_Id\tQuery_Coverage\tLen_Ratio\n";

    my %drop;
    my @keep;

    #  Fewer than two consensuses: nothing to compare against.
    if (keys(%cons) < 2)
    {
        @keep = @labels;
        for my $lab (@labels)
        {
            printf $rep "%s\t%d\t%d\tkeep\t\t\t\t\n",
                $lab, scalar @{$aligned->{$lab}}, length($cons{$lab} || "");
        }
        close $rep;
        print STDERR "Truncation QC: ", scalar(@keep),
                     " cluster(s) kept, 0 dropped (fewer than two consensuses to compare)\n";
        return { keep => \@keep, drop => \%drop };
    }

    #  All-vs-all blastp on the consensuses.
    my $qc_dir = "truncation_qc";
    mkdir($qc_dir) unless -d $qc_dir;
    my $faa = "$qc_dir/consensus.faa";
    open(my $fh, ">$faa") or die "Cannot open $faa for writing: $!\n";
    for my $lab (sort keys %cons) { print $fh ">$lab\n$cons{$lab}\n" }
    close $fh;

    my $db  = "$qc_dir/consensus";
    my $tab = "$qc_dir/consensus_vs_self.tsv";
    system("makeblastdb -in $faa -dbtype prot -out $db >/dev/null 2>&1") == 0
        or die "makeblastdb failed on $faa\n";
    my $n_seq = scalar keys %cons;
    my $cmd = "blastp -query $faa -db $db"
            . " -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send bitscore'"
            . " -evalue 1e-5 -max_target_seqs $n_seq -num_threads 4 > $tab 2>/dev/null";
    system($cmd) == 0 or die "blastp failed during truncation QC\n";

    #  Best single-HSP query coverage for each ordered pair.
    my (%cov, %pid);
    open(my $tfh, "<$tab") or die "Cannot open $tab: $!\n";
    while (<$tfh>)
    {
        chomp;
        my ($q, $s, $pi, undef, $ql, undef, $qs, $qe) = split /\t/;
        next if !defined $qe;
        next if $q eq $s;
        my $qcov = ($qe - $qs + 1) / $ql;
        if (!defined $cov{$q}{$s} || $qcov > $cov{$q}{$s})
        {
            $cov{$q}{$s} = $qcov;
            $pid{$q}{$s} = $pi;
        }
    }
    close $tfh;

    #  Longest consensus first. A cluster is dropped only if it is a
    #  substring of an already-ACCEPTED (longer) cluster, so a chain
    #  A within B within C keeps only C, and two clusters can never drop
    #  each other.
    my @by_len = sort { length($cons{$b}) <=> length($cons{$a}) or $a cmp $b }
                 keys %cons;

    my %verdict;
    for my $cand (@by_len)
    {
        my $hit;
        for my $ref (@keep)
        {
            next unless length($cons{$cand}) < $len_ratio * length($cons{$ref});
            my $c = $cov{$cand}{$ref};
            my $p = $pid{$cand}{$ref};
            next unless defined $c && defined $p;
            next unless $c >= $min_qcov && $p >= $pident;
            $hit = $ref;
            last;
        }

        if ($hit)
        {
            $drop{$cand} = $hit;
            $verdict{$cand} = [ "TRUNCATED", $hit,
                                sprintf("%.1f", $pid{$cand}{$hit}),
                                sprintf("%.3f", $cov{$cand}{$hit}),
                                sprintf("%.3f", length($cons{$cand}) / length($cons{$hit})) ];
        }
        else
        {
            push @keep, $cand;
            #  Report the closest near-miss anyway. A silent "keep" hides the
            #  case where a real truncation sat just under the threshold, which
            #  is exactly how a mis-set qc_pident goes unnoticed.
            my ($bp, $bref) = (-1, "");
            for my $ref (@keep)
            {
                next if $ref eq $cand;
                next unless length($cons{$cand}) < length($cons{$ref});
                my $c = $cov{$cand}{$ref};
                my $p = $pid{$cand}{$ref};
                next unless defined $c && defined $p && $c >= $min_qcov;
                if ($p > $bp) { $bp = $p; $bref = $ref }
            }
            if ($bp >= 0)
            {
                $verdict{$cand} = [ "keep", "(best: $bref)",
                                    sprintf("%.1f", $bp),
                                    sprintf("%.3f", $cov{$cand}{$bref}),
                                    sprintf("%.3f", length($cons{$cand}) / length($cons{$bref})) ];
            }
            else
            {
                $verdict{$cand} = [ "keep", "", "", "", "" ];
            }
        }
    }

    #  Any label that produced no consensus at all is kept and noted.
    for my $lab (@labels)
    {
        next if exists $verdict{$lab};
        push @keep, $lab;
        $verdict{$lab} = [ "keep (no consensus)", "", "", "", "" ];
    }

    #  Report, and move rejected alignments out of the way.
    mkdir("truncated_alis") unless -d "truncated_alis";
    for my $lab (sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] or $a cmp $b } @labels)
    {
        my $v = $verdict{$lab};
        printf $rep "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n",
            $lab, scalar @{$aligned->{$lab}}, length($cons{$lab} || ""),
            $v->[0], $v->[1], $v->[2], $v->[3], $v->[4];

        if (exists $drop{$lab} && -e "alis/$lab.fa")
        {
            move("alis/$lab.fa", "truncated_alis/$lab.fa")
                or die "Cannot move alis/$lab.fa to truncated_alis/: $!\n";
        }
    }
    close $rep;

    @keep = sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] or $a cmp $b } @keep;

    print STDERR "Truncation QC: ", scalar(@keep), " clusters kept, ",
                 scalar(keys %drop), " dropped as truncated";
    print STDERR " (", join(", ", map { "$_ within $drop{$_}" } sort keys %drop), ")" if %drop;
    print STDERR "\n";

    return { keep => \@keep, drop => \%drop };
}

1;
