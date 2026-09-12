package FCP_Ali_Utils;

use strict;
use warnings;
use Math::Round;
use Cwd;
use Exporter qw(import);
use gjoseqlib;
use FCP_PSSM_Utils qw(create_pssm);
use FCP_Trunc_Utils qw(qc_truncated_clusters);

our @EXPORT_OK = qw(process_cluster_alignments process_alignment);

# Main function to process all cluster alignments
sub process_cluster_alignments 
{
    my ($seq_data, $options, $base) = @_;
    my $tmp = $seq_data->{tmp};
    
    # Get cluster files
    opendir(my $dir_h, "clusters") or die "Cannot open clusters directory: $!\n";
    #  Sort numerically. readdir order is filesystem order, which made the
    #  Curation_Report rows come out scrambled (10, 12, 8, 13, ...).
    my @clusterF = sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] }
                   grep {$_ !~ /^\./} readdir($dir_h);
    closedir($dir_h);

    # Create report file
    open(my $report_fh, ">Curation_Report") or die "Cannot open Curation_Report for writing: $!\n";
    #  NB the last two labels used to be swapped relative to the values the
    #  writer emits: it prints end_dash_removed first and the cut-column
    #  list second. Labels corrected here rather than reordering the data,
    #  so existing Curation_Report files stay readable.
    print $report_fh "Alignment\tAli_Len\tFirst_AA\tOriginal Seqs\tFinal Seqs\tLow_Occ_Cols_Cut\tN_Term_Cut\tC_Term_Cut\tDash_Starts_Removed\tSeq_w_End_Dashes_Removed\tCut_Cols\tIdentical_Collapsed\tNterm_Stop_Col\tNterm_Stop_Res\tNterm_Flag\n";

    #---------------------------------------------------------------------------
    #  PASS 1 -- align every cluster with MAFFT and keep the alignment.
    #
    #  This used to be a single loop that aligned, corrected and made a PSSM
    #  one cluster at a time. The truncation QC has to compare every cluster
    #  consensus against every other one, so it cannot run until all of them
    #  are aligned -- hence the split into two passes with the QC between.
    #---------------------------------------------------------------------------
    my %aligned;
    my @order;
    my @nterm_flagged;

    foreach my $file (@clusterF) 
    {
        my $ali_file = $file;
        $ali_file =~ s/\.fasta//g;
        
        print STDERR "\n\nAligning $file\n";
        
        # Align the sequences with MAFFT
        open(my $in_fh, "mafft --thread 24 --quiet --reorder $base/$tmp/clusters/$file |") or die "Cannot run MAFFT: $!\n";
        my @ali = &gjoseqlib::read_fasta($in_fh);
        close($in_fh);

        # Print the original alignment
        open(my $out_fh, ">alis/$ali_file.fa") or die "Cannot open alis/$ali_file.fa for writing: $!\n";
        &gjoseqlib::print_alignment_as_fasta($out_fh, @ali);
        close($out_fh);

        $aligned{$ali_file} = \@ali;
        push @order, $ali_file;
    }

    #---------------------------------------------------------------------------
    #  TRUNCATION QC -- before any correction. Drops clusters whose consensus
    #  is a substring of a longer cluster's consensus.
    #---------------------------------------------------------------------------
    my %skip;
    unless ($options->{no_trunc_qc}) 
    {
        my $qc = qc_truncated_clusters(\%aligned, $options);
        %skip = %{ $qc->{drop} };
    }

    #---------------------------------------------------------------------------
    #  PASS 2 -- correct the survivors and build their PSSMs.
    #---------------------------------------------------------------------------
    foreach my $ali_file (@order) 
    {
        if (exists $skip{$ali_file}) 
        {
            print STDERR "Skipping $ali_file: consensus is a substring of $skip{$ali_file} (see Truncation_Report)\n";
            next;
        }

        my @ali = @{ $aligned{$ali_file} };

        # Process the alignment
        my $ali_data = process_alignment(\@ali, $options);
        
        #  A translated protein starts at Met. If the trim stopped somewhere
        #  else the alignment is FLAGGED for review, not trimmed further --
        #  non-Met starts are rare but real, and hunting onward for a
        #  conserved Met would delete good protein.
        my $st = $ali_data->{nterm_stop};
        if ($st && $st->{res} ne "M") 
        {
            push @nterm_flagged, [$ali_file, $st->{col}, $st->{res},
                                  $st->{occ}, $st->{id}, $st->{kept},
                                  $ali_data->{n_seqs_final}];
        }
        
        # Write to report file
        write_alignment_report($report_fh, $ali_file, $ali_data);
        
        # Write corrected alignment
        open(my $out_fh, ">corrected_alis/$ali_file.fa") or die "Cannot open corrected_alis/$ali_file.fa for writing: $!\n";
        &gjoseqlib::print_alignment_as_fasta($out_fh, @{$ali_data->{final_ali}});
        close($out_fh);
        
        #  Two floors, not one. -m decides which clusters are worth FORMING;
        #  this decides which formed-and-curated alignments are worth KEEPING,
        #  and it is lower on purpose. Gating both on -m throws away an
        #  alignment that cleared the bar at clustering and then lost a
        #  sequence or two to curation -- Alpharhabdovirinae M cluster 18 and
        #  Dichorhavirus M cluster 3 both went 6 -> 4 against -m 5 and were
        #  written to corrected_alis with no profile, leaving an orphan
        #  alignment and lost coverage for the sake of one sequence.
        my $keep_floor = defined($options->{min_keep}) ? $options->{min_keep} : 2;
        $keep_floor = $options->{min_seqs} if $options->{min_seqs} < $keep_floor;
        if ($ali_data->{n_seqs_final} >= $keep_floor) 
        {
            create_pssm_for_alignment($ali_file, $ali_data, $options, $tmp);
        }
        else
        {
            print STDERR "  $ali_file: $ali_data->{n_seqs_final} sequence(s) after ",
                         "curation, below the keep floor of $keep_floor; no profile\n";
        }
    }
    
    close($report_fh);

    if (@nterm_flagged) 
    {
        open(my $fh, ">FLAGGED_NTERM") or die "Cannot open FLAGGED_NTERM: $!\n";
        print $fh "Alignment\tStop_Col\tStop_Res\tOccupancy\tIdentity_Among_Occupied"
                . "\tSeqs_With_Residue\tFinal_Seqs\n";
        for my $r (sort { $a->[0] <=> $b->[0] } @nterm_flagged) 
        {
            printf $fh "%s\t%d\t%s\t%.3f\t%.3f\t%d\t%d\n", @$r;
        }
        close($fh);
        print STDERR "N-terminal review: ", scalar(@nterm_flagged),
                     " alignment(s) stopped on a non-Met column (see FLAGGED_NTERM)\n";
    }
    return 1;
}

# Process a single alignment
sub process_alignment 
{
    my ($ali, $options) = @_;
    
    # Extract options
    my $frac_dash = $options->{frac_dash};
    my $f_insert = $options->{f_insert};
    my $frac_n_term = $options->{frac_n_term};
    my $frac_c_term = $options->{frac_c_term};
    my $start_dash = $options->{start_dash};
    my $end_dash = $options->{end_dash};
    my $end_frac_occ = $options->{end_frac_occ};
    
    # Calculate original sequence count
    my $n_seqs_orig = scalar @$ali;
    
    # Remove sequences with too many dashes
    my @ali2;
    for my $i (0..$#$ali) 
    {
        my $len = length($ali->[$i][2]);
        my $dashes = $ali->[$i][2] =~ tr/-//;
        if (($dashes/$len) <= $frac_dash) 
        {
            push @ali2, $ali->[$i];
        }
    }
    
    my $n_seqs_internal_dash = scalar @ali2;
    
    # Pack the alignment (remove empty columns)
    my @ali3 = &gjoseqlib::pack_alignment(@ali2);
    
    # Process column occupancy
    my $n_seqs = scalar @ali3;
    my $min_n = round($f_insert * $n_seqs);
    my %col_sum;  # number of non-dash characters
    my $aa_sum = {}; # aa count
    
    # Count characters in each column
    for my $i (0..$#ali3) 
    {
        my $sequence = uc $ali3[$i][2];
        my @bases = split("", $sequence);
        for my $j (0..$#bases) 
        {
            if ($bases[$j] =~ /\w/) 
            {
                $aa_sum->{$j}->{$bases[$j]}++;  # only counts letters not dashes
                $col_sum{$j}++;
            }
        }
    }
    
    # Add low occupancy columns to exclusion hash
    my %exclude;
    foreach my $col (keys %col_sum) 
    {
        my $f_occ = ($col_sum{$col}/$n_seqs);
        if ($f_occ < $f_insert) 
        {
            print STDERR "Excluding $col\t$f_occ\n";
            $exclude{$col} = 1;
        }
    }
    
    my $low_occ = keys %exclude;  # count of low occupancy columns removed
    my $ali_len = length($ali3[0][2]);
    
    # Trim from the N-terminus
    my $n_occ = defined $options->{n_occ} ? $options->{n_occ} : 0.40;
    my $nterm_stop;
    my $n_term_trimmed = 0;
    my $c_term_trimmed = 0;
    
    for my $k (0..$ali_len) 
    {
        my $aa_resR = $aa_sum->{$k};
        next unless $aa_resR; # Skip if no data for this position
        
        my @sorted_res = sort {$aa_resR->{$b} <=> $aa_resR->{$a}} keys(%$aa_resR);
        next unless @sorted_res; # Skip if no residues
        
        my $mc_res = $sorted_res[0];

        #  Conservation is scored among the sequences that HAVE a residue here,
        #  not over every sequence in the alignment. $aa_sum counts letters
        #  only, so dividing it by the total sequence count conflates two
        #  different things: a column that is genuinely unconserved, and a
        #  column that is perfectly conserved but gappy because other
        #  sequences are truncated. When the truncated form is the majority it
        #  outvotes the real N-terminus and the column is deleted -- which is
        #  how Alpharhabdovirinae G cluster 2 lost 13 columns that were
        #  91-100% conserved among the 78 sequences carrying them.
        #
        #  The occupancy floor ($n_occ) is what stops a sparsely-populated
        #  leader halting the trim: a column has to be present in enough
        #  sequences to be worth keeping, not merely agree with itself.
        my $occ  = $col_sum{$k} || 0;
        my $cons = $occ ? ($aa_sum->{$k}->{$mc_res} / $occ) : 0;
        my $occf = $occ / $n_seqs_internal_dash;

        if ($occf < $n_occ || $cons < $frac_n_term) 
        {
            $exclude{$k} = 1;
            $n_term_trimmed++;
            printf STDERR "Trimmed %d %s occ=%.2f id=%.2f from N-terminal\n",
                          $k, $mc_res, $occf, $cons;
        } 
        else 
        {
            #  Stop here. A translated protein starts at Met, so if this column
            #  is not one, the alignment is flagged for review rather than
            #  trimmed further -- non-Met starts are rare but real, and hunting
            #  onward for a conserved Met would delete good protein.
            $nterm_stop = { col => $k, res => $mc_res,
                            occ => $occf, id => $cons, kept => $occ };
            last;
        }
    }
    
    # Trim from the C-terminus
    for (my $l = ($ali_len - 1); $l >= 0; $l--) 
    {
        my $aa_resR = $aa_sum->{$l};
        next unless $aa_resR; # Skip if no data for this position
        
        my @sorted_res = sort {$aa_resR->{$b} <=> $aa_resR->{$a}} keys(%$aa_resR);
        next unless @sorted_res; # Skip if no residues
        
        my $mc_res = $sorted_res[0];
        #  Same denominator fix as the N-terminus: score among the sequences
        #  that have a residue here. This trim never fires at the default
        #  frac_c_term of 0, but it carried the identical flaw.
        my $occ_c = $col_sum{$l} || 0;
        my $cons  = $occ_c ? ($aa_sum->{$l}->{$mc_res} / $occ_c) : 0;
        
        if ($cons < $frac_c_term) 
        {
            $exclude{$l} = 1;
            $c_term_trimmed++;
            print STDERR "Trimmed $l, $mc_res, $aa_sum->{$l}->{$mc_res}, $cons, from C-terminal\n";
        } 
        else 
        {
            last;
        }
    }
    
    #  Trim rather than delete when deletion would be expensive.
    #
    #  The N-terminal trim stops at the first column that clears -n_occ, and
    #  the loop below then deletes every sequence that does not reach that
    #  column. Those two steps can disagree about what a majority is: with the
    #  default -n_occ 0.40, a column present in 44% of sequences is kept and
    #  the other 56% are deleted.
    #
    #  Alpharhabdovirinae L cluster 11 was exactly that. Eight sequences began
    #  MKKTLTVIMDYSQEY..., ten began MDYSQEY... eight residues later -- both at
    #  a Met, a start-site annotation disagreement rather than truncation.
    #  Keeping the 8-residue extension cost ten of eighteen sequences.
    #
    #  So: if deleting the ragged starts would cost more than -max_start_frac
    #  of the alignment, trim the leading columns instead and keep everyone.
    #  Mirrors the C-terminal cap in the end-trim step below.
    my $max_start_frac = defined($options->{max_start_frac})
                       ? $options->{max_start_frac} : 0.35;
    unless ($start_dash)
    {
        my $n_all = scalar @ali3;
        my $would_lose = 0;
        my $need = 0;                      # columns to drop to rescue them
        for my $i (0..$#ali3)
        {
            my @b = split("", $ali3[$i][2]);
            my $j = 0;
            $j++ while ($j <= $#b && exists $exclude{$j});
            next if ($j > $#b || $b[$j] =~ /\w/);   # already starts on a residue
            $would_lose++;
            my $k = $j;
            $k++ while ($k <= $#b && ($b[$k] !~ /\w/ || exists $exclude{$k}));
            $need = $k if $k > $need;
        }
        if ($would_lose && $n_all && ($would_lose / $n_all) > $max_start_frac)
        {
            $exclude{$_} = 1 for (0 .. $need - 1);

            #  Trim to the COMMON MET START. Column $need is where the shorter
            #  sequences begin; when that is a Met -- the usual case, because
            #  the disagreement is over which of two in-frame Mets was called
            #  -- the trimmed alignment starts at a real initiation codon and
            #  nothing more is needed.
            #
            #  When it is not a Met, do NOT scan forward for one. That would
            #  delete real protein from every sequence in the cluster to
            #  satisfy the expectation, which is the same mistake the
            #  N-terminal trim was corrected for. Keep the trim, and flag the
            #  alignment so a curator looks at it.
            my %col0;
            for my $i (0..$#ali3)
            {
                my $ch = substr($ali3[$i][2], $need, 1);
                $col0{$ch}++ if $ch =~ /\w/;
            }
            my ($top) = sort { $col0{$b} <=> $col0{$a} || $a cmp $b } keys %col0;
            $top = "-" unless defined $top;
            my $occ_here = 0; $occ_here += $_ for values %col0;

            printf STDERR "  Start-trim: deleting the ragged starts would cost %d of "
                        . "%d sequences (%.0f%%), above the %.0f%% cap -- trimming %d "
                        . "leading column(s) to the common start '%s' instead%s\n",
                        $would_lose, $n_all, 100*$would_lose/$n_all,
                        100*$max_start_frac, $need, $top,
                        ($top eq "M" ? "" : "  [NOT a Met -- FLAGGED for review]");

            $nterm_stop = { col => $need, res => $top,
                            occ => ($n_all ? $occ_here / $n_all : 0),
                            id  => ($occ_here ? $col0{$top} / $occ_here : 0),
                            kept => $occ_here };
        }
    }

    # Apply exclusions and remove sequences starting with a dash if needed
    my @ali4;
    my $dash_starts_removed = 0;
    
    for my $i (0..$#ali3) 
    {
        my $id = $ali3[$i][0];
        my $anno = $ali3[$i][1];
        my @bases = split("", $ali3[$i][2]);
        my $string = "";
        
        for my $j (0..$#bases) 
        {
            unless (exists $exclude{$j}) 
            {
                $string .= $bases[$j];
            }
        }
        
        if ($start_dash) 
        {
            push @ali4, ([$id, $anno, $string]);
        } 
        elsif ($string !~ /^-/) 
        {
            push @ali4, ([$id, $anno, $string]);
        } 
        else 
        {
            $dash_starts_removed++;
        }
    }
    
    # Pack the alignment again
    my @ali5 = &gjoseqlib::pack_alignment(@ali4);
    
    # Get rid of identical sequences
    #  De-duplicate in place, PRESERVING the order MAFFT --reorder produced.
    #  The previous version built @ali6 by iterating `keys %$unique`, i.e. in
    #  Perl's randomized hash order. That threw away the similarity ordering
    #  that --reorder exists to create, and made the order differ on every
    #  run -- so corrected_alis never diffed clean, and because final_ali is
    #  what create_pssm consumes, the PSSMs were not reproducible either.
    #  Note: on a duplicate sequence this now keeps the FIRST occurrence
    #  (earliest in similarity order); the hash version kept the last.
    #  This is the point at which identical sequences are collapsed to a
    #  single member. It runs AFTER the column exclusions have been applied
    #  (ali4), so it catches sequences made identical by the N-terminal,
    #  C-terminal and low-occupancy trimming -- not just ones that arrived
    #  identical. Nothing further downstream can reintroduce a duplicate:
    #  ali7 only drops sequences, and the pack_alignment after it strips only
    #  columns that are all-gap in every survivor, which cannot merge two
    #  distinct sequences.
    #
    #  Order-preserving, first occurrence wins, so MAFFT's similarity order
    #  survives and the result is deterministic.
    my @ali6;
    my %seen;
    for my $i (0..$#ali5) 
    {
        my $seq = uc $ali5[$i][2];
        next if $seen{$seq}++;
        push @ali6, ([$ali5[$i][0], $ali5[$i][1], $seq]);
    }
    my $identical_collapsed = (scalar @ali5) - (scalar @ali6);

    #  A PSSM cannot be built from a single sequence -- PSI-BLAST has nothing
    #  to profile. If the collapse leaves exactly one member, put a second
    #  back so a PSSM can still be produced. The two are identical by
    #  definition, which is fine: the profile is then simply that sequence.
    if (@ali6 == 1 && @ali5 >= 2) 
    {
        push @ali6, ([$ali5[1][0], $ali5[1][1], uc $ali5[1][2]]);
        $identical_collapsed--;
        print STDERR "  Collapse left a single sequence; retained a second member so a PSSM can be built\n";
    }

    if ($identical_collapsed > 0) 
    {
        print STDERR "  Collapsed $identical_collapsed identical sequence(s) after trimming\n";
    }
    
    # Handle sequences with spurious C-terminal gaps
    #  Drop C-terminally truncated sequences, but never so many that the
    #  alignment is gutted.
    #
    #  The previous guard aborted the whole step unless every one of the last
    #  -e columns was at least (1 - end_frac_occ) occupied. That inverts under
    #  load: occupancy at the end is low precisely BECAUSE sequences are
    #  truncated, so the more raggedness there is, the less likely the step was
    #  to act on any of it. Alpharhabdovirinae L cluster 3 had 22 of 89
    #  sequences ending 10-38 residues early; the last columns were 75%
    #  occupied against an 85% threshold, so nothing was removed at all and the
    #  truncated proteins went into the profile.
    #
    #  The intent of the old guard was real -- a genuinely variable C-terminus
    #  should not be trimmed away -- so it is kept, expressed as a cap on how
    #  much may be removed rather than a veto on removing anything.
    my $len   = length($ali6[0][2]);
    my $nseqs = scalar @ali6;
    my $max_end_frac = defined($options->{max_end_frac})
                     ? $options->{max_end_frac} : 0.35;

    my @ragged = grep { substr($ali6[$_][2], -$end_dash) !~ /\w/ } (0..$#ali6);
    my $ragged_frac = $nseqs ? (scalar(@ragged) / $nseqs) : 0;

    my @ali7;
    if (!@ragged)
    {
        @ali7 = @ali6;
    }
    elsif ($ragged_frac > $max_end_frac)
    {
        #  Too many to be bad submissions: treat the short C-terminus as real.
        printf STDERR "  End-trim skipped: %d of %d sequences (%.0f%%) end early, "
                    . "above the %.0f%% cap -- treating the C-terminus as variable\n",
                    scalar(@ragged), $nseqs, 100*$ragged_frac, 100*$max_end_frac;
        @ali7 = @ali6;
    }
    elsif (($nseqs - scalar(@ragged)) < 2)
    {
        #  A profile needs two sequences; never trim below that.
        print STDERR "  End-trim skipped: removing the ragged sequences would "
                   . "leave fewer than 2\n";
        @ali7 = @ali6;
    }
    else
    {
        my %drop = map { $_ => 1 } @ragged;
        for my $i (0..$#ali6) { push @ali7, $ali6[$i] unless $drop{$i}; }
        printf STDERR "  End-trim: removed %d of %d sequences (%.0f%%) that end "
                    . "within the last %d columns\n",
                    scalar(@ragged), $nseqs, 100*$ragged_frac, $end_dash;
    }
    
    my @ali8 = &gjoseqlib::pack_alignment(@ali7);
    my $end_dash_removed = (scalar @ali6) - (scalar @ali7);

    
    # Final stats
    my $n_seqs_final = scalar @ali8;
    my $start_char = substr($ali8[0][2], 0, 1) if @ali8;
    my $final_ali_len = @ali8 ? length($ali8[0][2]) : 0;
    
    print STDERR "Original_Seqs = $n_seqs_orig, Final_Seqs = $n_seqs_final\n";
    print STDERR "Low Occ Cols Cut = $low_occ, N-term Cut = $n_term_trimmed, C-term Cut = $c_term_trimmed\t";
    print STDERR "Dash Starts Removed = $dash_starts_removed, Dash End Seqs Removed = $end_dash_removed\n\n";
    
    # Return all the data
    return {
        n_seqs_orig => $n_seqs_orig,
        n_seqs_final => $n_seqs_final,
        nterm_stop => $nterm_stop,
        identical_collapsed => $identical_collapsed,
        ali_len => $final_ali_len,
        start_char => $start_char,
        low_occ => $low_occ,
        n_term_trimmed => $n_term_trimmed,
        c_term_trimmed => $c_term_trimmed,
        dash_starts_removed => $dash_starts_removed,
        end_dash_removed => $end_dash_removed,
        exclude => \%exclude,
        final_ali => \@ali8
    };
}

# Write alignment processing report
sub write_alignment_report 
{
    my ($report_fh, $ali_file, $ali_data) = @_;
    
    print $report_fh join("\t", 
        $ali_file,
        $ali_data->{ali_len},
        $ali_data->{start_char},
        $ali_data->{n_seqs_orig},
        $ali_data->{n_seqs_final},
        $ali_data->{low_occ},
        $ali_data->{n_term_trimmed},
        $ali_data->{c_term_trimmed},
        $ali_data->{dash_starts_removed},
        $ali_data->{end_dash_removed},
        join(",", sort {$a <=> $b} keys %{$ali_data->{exclude}}),
        $ali_data->{identical_collapsed},
        ($ali_data->{nterm_stop} ? $ali_data->{nterm_stop}{col} : "-"),
        ($ali_data->{nterm_stop} ? $ali_data->{nterm_stop}{res} : "-"),
        ($ali_data->{nterm_stop} && $ali_data->{nterm_stop}{res} ne "M"
             ? "REVIEW-nonMet" : "ok")
    ), "\n";
    
    return 1;
}

# Create PSSM for an alignment
sub create_pssm_for_alignment 
{
    my ($ali_file, $ali_data, $options, $tmp) = @_;
    
    my $pssm_file = $ali_file;
    if ($options->{pssm_prefix}) 
    {
        $pssm_file = "$options->{pssm_prefix}.$ali_file";
    }
    
    # Set up PSSM options
    my $pssm_options = {
        outfile => "pssms/$tmp.pssm",
        final_file => "pssms/$pssm_file.pssm",
        process_titles => 1,
        cleanup => 1,
        tmp => $tmp,
        pssm_prefix => $options->{pssm_prefix}
    };
    
    # Create the PSSM using PSSMUtils
    create_pssm($ali_data->{final_ali}, $pssm_options);
    
    return 1;
}

1; # End of module