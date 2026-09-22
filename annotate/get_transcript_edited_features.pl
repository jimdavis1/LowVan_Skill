#! /usr/bin/env perl
use strict;
use Time::HiRes 'gettimeofday';
use GenomeTypeObject;
use File::Temp;
use JSON::XS;
use File::Slurp;
use IPC::Run qw(run);
use Cwd;
use gjoseqlib;
use Getopt::Long::Descriptive;

# Try to load version module; fall back to "dev" if not available
my $tool_version;
eval {
    require LowVanVersion;
    $tool_version = LowVanVersion::get_version();
};
if ($@ || !$tool_version) {
    $tool_version = "dev";
}


my $program_description = <<'END_DESCRIPTION';
This program performs feature calling for transcript edited proteins.  It reads and writes GTO files.
It works by using a set of hand-curated transcripts in --dir as the queries.

Both directions of edit are handled, and neither needs a flag -- the direction is a property of the
curated reference and is read off the alignment.  A reference carrying an EXTRA base relative to the
genome (a -1 ribosomal slip, or a non-templated insertion as in the paramyxovirus V/W proteins) gaps
the subject, and that base is inserted.  A reference that is one base SHORTER (a +1 slip, as in the
hepatitis C F/ARFP protein) gaps the query, and the corresponding genome base is dropped.  --gaps
counts both sides together.  --id, --gaps, and --cov refer 
to the strict inclusion criteria for enabling the mapping the nucleotides from the closest hand-curated transcript
onto the subject sequence.  When we enoucnter a blast match that is good, [defined by --lower_pid, --lower_pcov, and --eval], 
but not good enough to carry over the transcript-edited seqeunce, we call a partial_cds feature and the annotation becomes:
[Uncorrected . annotation string . encoding region]. It is considered a partial_cds feature because the translation would be 
be interrupted where the frame jump occurs, or shortly thereafter.

END_DESCRIPTION


my $default_data_dir = $ENV{LOWVAN_DATA_DIR} // "/home/jjdavis/bin/Viral_Annotation";

my ($help, $tmp); 
my($opt, $usage) = describe_options(
    "\n$program_description\nUsage: %c %o",  
				    ["input|i=s"            => "Input GTO"],
				    ["output|o=s"           => "Output GTO"],
				    ["cov|c=i"              => "Minimum BLASTn percent query coverage (D = 95)", { default => 95 }],
				    ["id|p=i"               => "Minimum BLASTn percent identity  (D = 95)", { default => 95 }],
				    ["gaps|g=i"             => "Maximum number of allowable gaps (D = 2)", { default => 2 }],
				    ["e_val|e=f"            => "Maximum BLASTn evalue for considering any HSP (D = 0.5)", { default => 0.5 }],
				    ["lower_pid|lpi=i"      => "Lower percent identity threshold for a feature call without transcript editing correction (D = 80)", {default => 80}],
				    ["lower_pcov|lpc=i"     => "Lower percent query coverage for for a feature call without transcript editing correction (D = 80)", {default => 80}],
				    ["threads|a=i"          => "Threads for the BLASTN (D = 24))", { default => 24 }],
				    ["json|j=s"             => "Full path to the JSON opts file", {default => "$default_data_dir/Viral_PSSM.json"}],
				    ["dir|d=s"              => "Full path to the directory hand curated transcripts", {default => "$default_data_dir/Transcript-Editing"}],
				    ["tmp|t=s"              => "Declare name for temp dir (D = randomly named in cwd)"], 
				    ["help|h"               => "Show this help message", { shortcircuit => 1 } ],
				    ["version|v"            => "Show version information", { shortcircuit => 1 } ],
				    ["debug|b"              => "Enable debugging"],
);

if ($opt->version) {
    print "get_transcript_edited_features.pl version $tool_version\n";
    exit 0;
}
print $usage->text and exit if $opt->help;
die($usage->text) if @ARGV != 0;

if ($opt->tmp){$tmp = $opt->tmp;}
#else{$tmp .= sprintf("%x", rand 16) for 1..20;}
else {$tmp = File::Temp->newdir(CLEANUP => ($opt->debug ? 0 : 1))}
print STDERR "Tempdir=$tmp\n" if $opt->debug;

my $dir = $opt->dir;

my $genome_in = GenomeTypeObject->create_from_file($opt->input);
$genome_in or die "Error reading and parsing input";
$genome_in->{features}->[0] or die "No features in GTO\n"; 

my $base = getcwd;


# Get the viral family from the GTO
my $fam = $genome_in->{viral_family};
$fam or die "GTO has no viral_family field (not annotated by LowVan?)\n";


# Next we read the JSON to see if there are any transcript edited features that we need to find
my $json      = decode_json(read_file($opt->json));
$genome_in or die "Error reading json protein feature data";

my @to_analyze;
foreach (keys %{$json->{$fam}->{features}})
{
	my $prot = $_;
	if ($json->{$fam}->{features}->{$prot}->{special} eq "transcript_edit")
	{
		my $anno = $json->{$fam}->{features}->{$prot}->{anno};
		my $symbol = $json->{$fam}->{features}->{$prot}->{gene_symbol};
		my $feature_type = $json->{$fam}->{features}->{$prot}->{feature_type};
		push @to_analyze, ([$prot, $anno, $feature_type, $symbol]);
	}
}

# If there are transcript edited features then we set up the blast
if (scalar @to_analyze)
{
	my %features;
	print STDERR "\n\nSearching for transcript-edited features\n------------------------------------\n"; 
	print STDERR "\n$fam detected from GTO\n\n"; 

	mkdir ($tmp); 
	chdir ($tmp);
	
	my $contigs = $genome_in->{id}."."."contigs"; 
	$genome_in->write_contigs_to_file($contigs);
	
	#make the blastn db in the temp dir.
	my $make_db = run("makeblastdb -dbtype nucl -in $contigs >/dev/null");

	if (!$make_db)
	{
  	 print STDERR "get_transcript_edited_features:  makeblastdb failed with rc=$?. Stdout:\n";
	}
	
	# create the GTO analysis event.
	my $event = {
 	   tool_name => "LowVan Transcript Edited Features",
 	   tool_version => $tool_version,
  	   execution_time => scalar gettimeofday,
	};
	my $event_id = $genome_in->add_analysis_event($event);

	
	
	#cycle through the transcript edited features and search for them one at a time with blastn.
	for my $i (0..$#to_analyze)
	{
		my $name   = $to_analyze[$i][0];
		my $anno   = $to_analyze[$i][1]; 
		my $ft     = $to_analyze[$i][2];
		my $symbol = $to_analyze[$i][3];

		
		print STDERR "\tAnalyzing $name\n\n"; 
		my $query = "$dir/$fam/$name.fasta"; 

		run ("cp $query ."); 
		

		my @blast_parms = (
		      "-query",         $query,
		      "-db",            $contigs, #### fix this in the original program
		      "-evalue",        $opt->e_val,
		      "-reward",          2,
		      "-penalty",        -3,
		      "-word_size",      28,
		      "-outfmt",         15,
		      "-soft_masking",   "false",
		      "-dust",           "no",
		      "-perc_identity",  $opt->lower_pid,   
		      "-qcov_hsp_perc",  $opt->lower_pcov,
		      "-num_threads",    $opt->threads);

		my $do_blast = run(["blastn", @blast_parms], ">", "$name.json", "2>", "$name.blastn.stderr.txt");


		# got stuck on blast chunking queries >10MB.
		# Had to merge the json blocks manually.
		open (IN, "<$name.json") or die "Cannot open JSON BLASTn output file $name.json\n";
		my $json_content = read_file(\*IN);
		close IN;
		
		# Split on "}{\n" pattern that separates root JSON objects
		my @json_strings = split(/\}\s*\{/, $json_content);
		
		# Re-add braces and parse each object
		my @results;
		for my $i (0..$#json_strings) 
		{
			my $json_str = $json_strings[$i];
			# Add opening brace if not first
			$json_str = '{' . $json_str unless $i == 0;
			# Add closing brace if not last
			$json_str = $json_str . '}' unless $i == $#json_strings;
			push @results, decode_json($json_str);
		}
		
		# Merge all BlastOutput2 arrays into one result structure
		my $merged = { BlastOutput2 => [] };
		foreach my $result (@results) 
		{
			push @{$merged->{BlastOutput2}}, @{$result->{BlastOutput2}};
		}
		
		my $results = $merged;
					
				
		#Gather in the best match.
		my $matches = best_blastn_match_by_loc($results);  #removed id, cov, and gap thresholds from here  
		
		foreach (keys %$matches)
		{
			my $sid = $_;
			foreach (keys %{$matches->{$sid}})
			{
				my $from    = $_; 
				my $to      = $matches->{$sid}->{$from}->{TO};
				my $sseq    = $matches->{$sid}->{$from}->{HSEQ};
				my $qseq    = $matches->{$sid}->{$from}->{QSEQ};
				my $iden    = $matches->{$sid}->{$from}->{IDEN};
				my $ali_len = $matches->{$sid}->{$from}->{ALI_LEN};				
				#  Gaps are counted on BOTH sides, because the two directions of
				#  frameshift put the gap on opposite sides of the alignment and
				#  both are real biology:
				#
				#    -1 slip / non-templated insertion (NS1', V, W, ORF1ab)
				#         the reference carries an EXTRA base, so the SUBJECT
				#         gaps and the fill inserts the reference base.
				#    +1 slip (HCV F / ARFP)
				#         the ribosome skips a base, so the reference is one base
				#         SHORTER, the QUERY gaps, and the fill must DROP the
				#         corresponding subject base.
				#
				#  Only the first was implemented, which meant a +1 product
				#  silently reproduced the genome's own reading frame instead of
				#  the edited one -- the loop copied the subject through and
				#  nothing signalled that anything was wrong. references/
				#  special-features.md already claimed one mechanism served both.
				my $sdash   = ($sseq =~ tr/-//);
				my $qdash   = ($qseq =~ tr/-//);
				my $dashes  = $sdash + $qdash;
				my $sruns   = (() = $sseq =~ /-+/g) || 0;
				my $qruns   = (() = $qseq =~ /-+/g) || 0;
				my $runs    = $sruns + $qruns;
				my $pid     =  (($iden/$ali_len) * 100);				
				my $qcov    =  (($ali_len/(length $qseq)) * 100); 				
				my $gaps    = $matches->{$sid}->{$from}->{GAPS};	
						
				
				#  $corrected records whether a trustworthy protein was actually produced.
				#  The four gates below describe the NUCLEOTIDE match; none of them looks
				#  at what the gap-fill translated to. Since the fill takes its bases from
				#  the SUBJECT, an N-masked target yields an X-bearing protein. Screening
				#  only the alignment let 31 of 293 Togaviridae TF calls ship with an X in
				#  them and one with an internal stop; no amount of reference curation can
				#  prevent that, because the defect is in the genome being annotated.
				my $corrected = 0;

				#If all inclusion critreria are met (%id, %Qcov, num gaps, runs of gaps)
				
				#  A reference may edit in ONE direction only. A single frameshift
				#  is either +1 or -1; a reference asking for both an insertion
				#  and a deletion is describing two events and is far more likely
				#  to be a mis-curated reference or a spurious alignment than a
				#  real double frameshift. $runs <= 1 already excludes it, since a
				#  subject gap and a query gap are necessarily separate runs, but
				#  that is emergent from summing the two and a later reader could
				#  undo it without noticing. State it.
				my $mixed = ($sdash > 0 && $qdash > 0);
				if ($mixed)
				{
					print STDERR "\t$name: reference asks for an insertion AND a deletion "
					           . "($sdash subject, $qdash query); refusing, an edit is one "
					           . "direction only\n";
				}

				if ( !$mixed && ($pid >= $opt->id) && ($qcov >= $opt->cov) && ($dashes <= $opt->gaps) && ($runs <= 1))
				{
					my @snts = split ("", $sseq);
					my @qnts = split ("", $qseq);
					my @mod_seq;

					#  Build the edited transcript from the SUBJECT, so that every
					#  base except the edit itself comes from the genome being
					#  annotated and each genome yields its own protein.
					#
					#    subject gap  -> the reference has a base the genome lacks;
					#                    insert it
					#    query gap    -> the genome has a base the reference lacks;
					#                    skip it
					#    otherwise    -> take the genome's base
					#
					#  No flag selects between these. The direction is a property
					#  of the reference that was curated, so the right behaviour is
					#  whichever the alignment shows.
					for my $i (0..$#snts)
					{
						if ($snts[$i] =~ /\-/)
						{
							push @mod_seq, $qnts[$i];
						}
						elsif ($qnts[$i] =~ /\-/)
						{
							next;
						}
						else
						{
							push @mod_seq, $snts[$i]; 
						}
					}

					my $mod = join ("", @mod_seq);
					my $mod_aa = &gjoseqlib::translate_seq( $mod );
			
					my ($len, $strand);
					if ($from < $to)
					{
						$strand = "+";
						$len = ($to - $from) + 1;
					}
					elsif ($from > $to){
						$strand = "-";
						$len = ($from - $to) + 1;	
					}
			
					my $feature = {
						type        => $ft,
						contig      => $sid,
						aa_sequence => $mod_aa,
						location    => ([[$sid, $from, $strand, $len]]),
						product     => $anno,
						symbol      => $symbol,
						pssm        => ([[$fam, $name, $anno, "LowVan Transcript Edited Feature"]]),
					};
					
					#  Screen the PROTEIN before emitting it. A failure here is not a
					#  dropped call: it falls through to the partial_cds branch below,
					#  which is exactly what that branch is for -- the region is located
					#  but no translation can be trusted across the frame jump.
					my $aa_body = $mod_aa;
					$aa_body =~ s/\*$//;          # a terminal stop is the ORF's own
					if ($aa_body =~ /X/i)
					{
						print STDERR "\t$name: gap-filled sequence translates with ambiguous residues; demoting to partial_cds\n";
					}
					elsif ($aa_body =~ /\*/)
					{
						print STDERR "\t$name: gap-filled sequence translates with an internal stop; demoting to partial_cds\n";
					}
					else
					{
						push(@{$features{$ft}}, $feature);
						$corrected = 1;
					}
				}
			
			
				# If the inclusion criteria were NOT met, or they were met but the
				# resulting protein did not survive the screen above, and there is
				# still a decent HSP from the blast, we add it as a partial CDS that
				# goes uncorrected.  No protein translation is given.
				if (! $corrected)
				{
					my $feature_type = "partial_cds";
					my ($len, $strand);
					if ($from < $to)
					{
						$strand = "+";
						$len = ($to - $from) + 1;
					}
					elsif ($from > $to){
						$strand = "-";
						$len = ($from - $to) + 1;	
					}

					my $lc_anno = lcfirst($anno);
				
					my $feature = {
						type        => $feature_type,
						contig      => $sid,
						location    => ([[$sid, $from, $strand, $len]]),
						product     => "Uncorrected "."$lc_anno"." encoding region",
						pssm        => ([[$fam, $name, $anno, "LowVan Transcript Edited Feature"]]),  #I don't actually use this but i left it there.
					};
					push(@{$features{$feature_type}}, $feature);
				}
			}
		}
	}	
	
	# Push features into the GTO
	if (%features)
	{
		foreach (keys %features)
		{
			my $type = $_; 
					
			foreach (@{$features{$type}})
			{
				my $data = $_;
	
				if ($type eq 'CDS' || $type eq 'mat_peptide')
				{
					my $p = {
						-id	                 => $genome_in->new_feature_id($type),
						-type 	             => $type,
						-location 	         => $data->{location},
						-analysis_event_id 	 => $event_id,
						-annotator           => 'LowVan Transcript Edited Feature',
						-protein_translation => $data->{aa_sequence},
						-function            => $data->{product},
						-family_assignments  => $data->{pssm},
						};				
					if (defined $data->{symbol} && $data->{symbol} ne '') 
					{
						$p->{-alias_pairs} = [[gene => $data->{symbol}]];
					}
					$genome_in->add_feature($p);
				}
				
				# Call a partial cds and do not add the AA seq if its a distant match.
				# No family assignment is generated
				elsif ($type eq 'partial_cds')
				{
					my $p = {
						-id	                 => $genome_in->new_feature_id($type),
						-type 	             => $type,
						-location 	         => $data->{location},
						-analysis_event_id 	 => $event_id,
						-annotator           => 'LowVan Transcript Edited Feature',
						-function            => $data->{product},
						};
					$genome_in->add_feature($p);
				} 	
			}
		}
		chdir ($base);
		$genome_in->destroy_to_file($opt->output);			
	}
 	else
 	{
		#handle condition where there should have been a blast match, but none was found.
		print STDERR "Transcript edited features curated for: $fam, but none were found.\n";
		chdir($base);
		$genome_in->destroy_to_file($opt->output);			
 	}
 }


 else 
 {
 	print STDERR "No proteins from transcript editing for: $fam\n";
	chdir($base);
	$genome_in->destroy_to_file($opt->output);			
 }
 
 







##########################sub best_blastn_match_by_loc##############
# Find the best blast match(s) for a blastn json
#
#   This will return the best blastn matched sequence by location.  Only one HSP within 
#   A distance defined as > abs(alignment length) is returned per contig.  
#   (i.e., you can have more than one "best" blast match per contig, but they
#   can't overlap.  The amount of overlap could be turned into a parameter, 
#   but i just used the ali length, which seemed reasonable. 
#
#   NOTE: this reads the -db formatted json output not the -subject [fasta] formatted version. 
#   The returned JSONs are slightly different 
#
#   Relevant paramaters such as %id %coverage should be set in the blast command line options,
#   or post-processed.
#
#   usage:
#   $hash = best_blastn_match_by_loc($blastn_json);
# 
#   The returned hash reference is in the following format:
#     
#    SubjectID->{HitFrom}->{TO}->{HitTo} 
#             ->{HitFrom}->{HSEQ}->{SubjectSeq}
#             ->{HitFrom}->{QSEQ}->{QuerySeq}
#             ->{HitFrom}->{BIT}->{BitScore} 
#             ->{HitFrom}->{IDEN}->{NumIdentities}         ### need to add this
#             ->{HitFrom}->{ALI_LEN}->{Alignment_Length}   ### need to add this
#             ->{HitFrom}->{GAPS}->{gaps}
#
# 
#----------------------------------------------------------
sub best_blastn_match_by_loc
{
	my ($blast)  =  @_;
	
	my $matches = {};
	
	my @output = @{$blast->{BlastOutput2}};
	for my $i (0..$#output)
	{
		my @hits = @{$blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}};
		for my $j (0..$#hits)
		{
			my @hsps = @{$blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}};
			for my $k (0..$#hsps)
			{
				my $ident     = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{identity};
				my $ali_len   = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{align_len};
				my $qlen      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{query_len};
				my $gaps      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{gaps};
				my $sid       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{description}->[0]->{title}; #changed from $k per go analysis
				my $hit_from  = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hit_from};
				my $bit       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{bit_score};
				my $hit_to    = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hit_to};
				my $qseq      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{qseq};
				my $hseq      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hseq};
								
				my ($pid, $qcov);
				if ($ident && $ali_len && $qlen) # this ensures that we got search results.
				{
					# Look for a match already recorded on this contig that is close
					# enough (proximity < this HSP's alignment length) to be treated as
					# the SAME location.  Keyed on proximity, NOT on "is the contig seen",
					# so a genuinely separate location on an already-seen contig is kept.
					my $near_loc;
					if (exists $matches->{$sid})
					{
						foreach my $loc (keys %{$matches->{$sid}})
						{
							if (abs($hit_from - $loc) < $ali_len)
							{
								$near_loc = $loc;
								last;
							}
						}
					}

					if (defined $near_loc)
					{
						# Overlapping location already recorded: keep the better bit score.
						if ($bit > $matches->{$sid}->{$near_loc}->{BIT})
						{
							delete $matches->{$sid}->{$near_loc};
							$matches->{$sid}->{$hit_from}->{BIT}      = $bit;
							$matches->{$sid}->{$hit_from}->{TO}       = $hit_to;
							$matches->{$sid}->{$hit_from}->{QSEQ}     = $qseq;
							$matches->{$sid}->{$hit_from}->{HSEQ}     = $hseq;
							$matches->{$sid}->{$hit_from}->{IDEN}     = $ident;
							$matches->{$sid}->{$hit_from}->{ALI_LEN}  = $ali_len;
							$matches->{$sid}->{$hit_from}->{GAPS}     = $gaps;
						}
						# else: existing match is better or equal -> discard this HSP
					}
					else
					{
						# No nearby match on this contig (contig unseen, or seen but this
						# HSP is at a distinct location): record it as a NEW entry.
						$matches->{$sid}->{$hit_from}->{BIT}      = $bit;
						$matches->{$sid}->{$hit_from}->{TO}       = $hit_to;
						$matches->{$sid}->{$hit_from}->{QSEQ}     = $qseq;
						$matches->{$sid}->{$hit_from}->{HSEQ}     = $hseq;
						$matches->{$sid}->{$hit_from}->{IDEN}     = $ident;
						$matches->{$sid}->{$hit_from}->{ALI_LEN}  = $ali_len;
						$matches->{$sid}->{$hit_from}->{GAPS}     = $gaps;
					}
				}
			}
		}
	}
	return $matches;
}

###########################################################













