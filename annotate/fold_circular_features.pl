#! /usr/bin/env perl
use strict;
use GenomeTypeObject;
use Getopt::Long::Descriptive;
use Time::HiRes 'gettimeofday';

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
Fold a circular genome back to one clean replicon.  It reads and writes GTO files and is the
LAST stage of the annotation pipeline, after the special-feature callers and before
viral_genome_quality.pl.

WHY THIS EXISTS.  A circular genome reaches GenBank as the loop cut at one arbitrary point, so
a CDS lying across that cut arrives split in two: part at the end of the contig, part at the
start.  Searched as deposited it is two partial hits, neither clearing its bit_cutoff, and the
feature is simply lost.  Measured on 400 near-full-length hepatitis B deposits, the polymerase
crosses the origin in 97.2 per cent of them and the large envelope protein in all of them, so for HBV this is
not an edge case, it is the normal state of the data.

annotate_by_viral_pssm.pl therefore DOUBLES a contig whose module declares
`replicon_geometry: circular` and searches the doubled copy, where a feature crossing the cut
appears whole and in order, and where every downstream step -- stop scanning, Met scanning,
subsequence extraction, translation -- works unchanged because nothing is circular any more.
annotate_by_viral_pssm-GTO.pl doubles the contig in the GTO to match, and the special-feature
callers inherit the doubled view for free, which is what lets a SPLICED product across the
origin be found at all.

This program undoes all of it.  It folds coordinates past the true length back into two-part
join locations, drops the second copy of everything, restores each contig to its deposited
length, and leaves a genome that the quality checker and every reader after it can treat as
ordinary.  Nothing doubled is allowed past this point.

WHAT IT DOES, in order, for each contig carrying `circular_true_length`:

  1. Folds EACH PIECE of each location, not each feature.  A two-exon gene cut partway through
     its first exon folds to three pieces, and the per-piece rule produces that without knowing
     anything about splicing.
  2. Drops any feature longer than the replicon.  On a doubled contig every feature appears
     twice, so a search can pair the 5' end of one copy with the 3' end of the next and report
     one feature spanning a whole genome length.  A real feature on a circle of length L cannot
     exceed L, so this single comparison removes that entire class of artifact -- including the
     spliced version of it, a false "intron" one genome long, which is why the test has to be
     here and not only in the annotator.
  3. Removes duplicates.  Everything is found once per copy; two features that fold to the same
     location are the same feature.  This is arithmetic, not a judgement about which call looks
     better, and a feature genuinely present twice survives because its two copies fold to
     different locations.
  4. Restores the contig to its deposited sequence and length.

WHAT IT DOES NOT DO.  It never reorders the pieces of a location.  On a circle cut at an
arbitrary point the second exon of a gene can have a LOWER coordinate than the first, so sorting
the pieces into ascending order would silently reverse the exons and yield a protein that is
wrong but still looks plausible.  The 5'->3' order produced upstream is preserved exactly.

A contig with no `circular_true_length` is left untouched, so running this on a linear genome,
or on a genome annotated before circular support existed, changes nothing.
END_DESCRIPTION

my($opt, $usage) = describe_options(
    "\n$program_description\nUsage: %c %o",
				    ["input|i=s"            => "Input GTO"],
				    ["output|o=s"           => "Output GTO"],
				    ["keep_doubled|k"       => "Fold the coordinates but leave the contigs doubled (debugging only)"],
				    ["help|h"               => "Show this help message", { shortcircuit => 1 } ],
				    ["version|v"            => "Show version information", { shortcircuit => 1 } ],
				    ["debug|b"              => "Enable debugging"],
);

print($usage->text), exit 0 if $opt->help;
if ($opt->version) { print "fold_circular_features.pl version $tool_version\n"; exit 0; }

my $genome_in = GenomeTypeObject->create_from_file($opt->input);
$genome_in or die "Error reading and parsing input";

#  contig id -> true (deposited) length, for the contigs that were doubled.
my %true_len;
for my $ctg (@{$genome_in->{contigs} || []})
{
	next unless $ctg->{circular_true_length};
	$true_len{$ctg->{id}} = $ctg->{circular_true_length};
}

unless (%true_len)
{
	#  Nothing was doubled. Pass the genome through untouched rather than
	#  rewriting it, so a linear family's output is byte-identical to what it
	#  was before this stage existed.
	$genome_in->destroy_to_file($opt->output);
	exit 0;
}

#  RECORD THE FOLD AS AN ANALYSIS EVENT.
#
#  Every other stage of this pipeline registers one, and this stage DELETES
#  features and REWRITES locations -- it is the only stage that removes a call
#  another stage made. Without an event there is nothing in the finished GTO to
#  say that happened, and a genome that was folded is indistinguishable from one
#  that was never doubled: the fold clears circular_true_length on purpose, so
#  that it is idempotent, and that is also the only field that recorded the
#  doubling. The event is what survives.
my $event_id = $genome_in->add_analysis_event({
	tool_name      => "LowVan Fold Circular Features",
	tool_version   => $tool_version,
	execution_time => scalar gettimeofday,
});

my (@to_delete, $n_folded, $n_dup, $n_over, $n_contained);

#  CHOOSE ONE INSTANCE PER FEATURE, BEFORE FOLDING ANYTHING.
#
#  On a doubled contig every feature is found twice, and the two copies are NOT
#  interchangeable. A copy that runs off the end of the doubled sequence comes
#  back truncated, and a truncated copy folds to a DIFFERENT, shorter location
#  than the whole one -- so comparing folded locations keeps both and emits a
#  short bogus feature beside the correct one. Measured on a rubella genome
#  rotated by 8000 nt: E2 was found whole at 9174-10019 (crossing the origin)
#  and again truncated at 18936-19523, cut off by the end of the doubled
#  sequence at 19524.
#
#  Which copy is the whole one depends on the strand -- a plus-strand feature
#  is complete in the first copy, a minus-strand one reads downward and is
#  complete in the second -- so the rule cannot be "keep the first". Instead,
#  group the copies by where the feature STARTS on the replicon, which is the
#  same for both, and keep the longest. That is one rule for both strands and
#  for truncation and duplication alike.
#
#  The group key carries the function as well, so two genuinely different
#  features that happen to share a start -- nested ORFs do this routinely; in
#  this same rubella genome p150 and p200 both start at 1803 -- are never
#  collapsed into one.
{
	my %group;
	for my $feat (@{$genome_in->{features} || []})
	{
		my $loc = $feat->{location} or next;
		next unless ref $loc eq 'ARRAY' && @$loc;
		my $L = $true_len{$loc->[0][0]} or next;

		my $start  = $loc->[0][1];
		my $anchor = (($start - 1) % $L) + 1;
		my $total  = 0;
		$total += $_->[3] for @$loc;

		my $key = join("|", $feat->{type} // '', $loc->[0][0], $loc->[0][2] // '',
		                    $anchor, $feat->{function} // '');
		push @{$group{$key}}, [$feat, $total];
	}

	for my $key (sort keys %group)
	{
		my @inst = @{$group{$key}};
		next if @inst < 2;
		my ($keep) = sort { $b->[1] <=> $a->[1] } @inst;   # longest wins
		for my $i (@inst)
		{
			next if $i->[0] == $keep->[0];
			push @to_delete, $i->[0]{id};
			$n_dup++;
		}
	}
}
my %dropped = map { $_ => 1 } @to_delete;

for my $feat (@{$genome_in->{features} || []})
{
	my $loc = $feat->{location} or next;
	next unless ref $loc eq 'ARRAY' && @$loc;
	next if $dropped{$feat->{id} // ''};

	my $contig = $loc->[0][0];
	my $L      = $true_len{$contig} or next;

	my @folded;
	my $wrapped = 0;
	for my $piece (@$loc)
	{
		my ($c, $start, $strand, $len) = @$piece;
		my @out = fold_piece($c, $start, $strand, $len, $L);
		$wrapped++ if @out > 1;
		push @folded, @out;
	}

	#  Total length is the sum of the pieces: a spliced feature's introns do
	#  not count toward what it occupies on the replicon, but no feature, split
	#  or not, can have more coding sequence than the circle has bases.
	my $total = 0;
	$total += $_->[3] for @folded;
	if ($total > $L)
	{
		push @to_delete, $feat->{id};
		$n_over++;
		warn "$feat->{id} ($feat->{function}) totals $total nt on a circle of $L nt -- ",
		     "impossible; dropped as a duplicate-copy artifact of the doubled contig\n";
		next;
	}

	$feat->{location} = \@folded;
	$n_folded++ if $wrapped;
}

#  CONTAINMENT: the partial call that the whole one swallows.
#
#  Doubling lets a gene across the origin be found whole, but it does not stop
#  the search ALSO finding the leftover piece on its own, and that piece is a
#  legitimate hit in its own right -- it clears bit_cutoff, it translates, and
#  nothing above this point has any reason to doubt it.  On the hepatitis B
#  reference NC_003977.2 the surface ORF produced both the correct
#  join(2850..3182,1..837) and a bare 1..837, which is the same protein's last
#  279 residues called a second time.  Grouping by where the feature starts
#  cannot catch it: the two copies start in different places, because one of
#  them is not a copy, it is a fragment.
#
#  So, last of all and on the FOLDED coordinates, drop any feature whose bases
#  are wholly inside another call of the same type and function on the same
#  strand.  Containment is the test rather than overlap, because overlap is
#  ordinary in this family -- every hepadnavirus gene overlaps another -- while
#  one call of a protein sitting entirely inside another call of the SAME
#  protein is never anything but a duplicate.
{
	my %by;
	for my $feat (@{$genome_in->{features} || []})
	{
		my $loc = $feat->{location} or next;
		next unless ref $loc eq 'ARRAY' && @$loc;
		next unless $true_len{$loc->[0][0]};
		next if $dropped{$feat->{id} // ''};

		my %pos;
		for my $p (@$loc)
		{
			my ($c, $st, $sd, $ln) = @$p;
			if (($sd // '+') eq '-') { $pos{$_} = 1 for ($st - $ln + 1) .. $st }
			else                     { $pos{$_} = 1 for $st .. ($st + $ln - 1) }
		}
		my $key = join("|", $feat->{type} // '', $loc->[0][0],
		                    $loc->[0][2] // '', $feat->{function} // '');
		push @{$by{$key}}, [$feat, \%pos, scalar keys %pos];
	}

	for my $key (sort keys %by)
	{
		my @g = @{$by{$key}};
		next if @g < 2;
		for my $i (0 .. $#g)
		{
			my $fid = $g[$i][0]{id} // '';
			next if $dropped{$fid};
			for my $j (0 .. $#g)
			{
				next if $i == $j || $dropped{$g[$j][0]{id} // ''};
				next unless $g[$j][2] > $g[$i][2];      # only a BIGGER call can contain
				my $inside = 1;
				for my $p (keys %{$g[$i][1]})
				{
					next if $g[$j][1]{$p};
					$inside = 0;
					last;
				}
				next unless $inside;
				push @to_delete, $fid;
				$dropped{$fid} = 1;
				$n_contained++;
				warn "$fid ($g[$i][0]{function}) lies wholly inside a longer call of ",
				     "the same feature; dropped as the leftover piece of a gene the ",
				     "doubled contig recovered whole\n";
				last;
			}
		}
	}
}

for my $fid (@to_delete)
{
	if ($genome_in->can('delete_feature')) { $genome_in->delete_feature($fid); }
	else { @{$genome_in->{features}} = grep { $_->{id} ne $fid } @{$genome_in->{features}}; }
}

#  Restore each doubled contig to the sequence that was deposited.
unless ($opt->keep_doubled)
{
	for my $ctg (@{$genome_in->{contigs} || []})
	{
		my $L = $true_len{$ctg->{id}} or next;
		$ctg->{dna}    = substr($ctg->{dna}, 0, $L);
		$ctg->{length} = $L if exists $ctg->{length};
		delete $ctg->{circular_true_length};
		#  replicon_geometry is left in place: it is true of the replicon, and
		#  viral_genome_quality.pl sets the same field from the module anyway.
	}
}

print STDERR "Circular fold: ", scalar(keys %true_len), " contig(s); ",
             ($n_folded // 0), " feature(s) folded across the origin; ",
             ($n_dup // 0), " duplicate(s) removed; ",
             ($n_over // 0), " over-length artifact(s) dropped; ",
             ($n_contained // 0), " contained partial(s) dropped\n";

$genome_in->destroy_to_file($opt->output);


##############sub fold_piece####################
#  One piece of a location, in DOUBLED coordinates, to one or two pieces in the
#  replicon's own coordinates.
#
#  GTO convention: start is the 5'-most base of the piece and len runs in the
#  reading direction, so on the minus strand start is the HIGHER coordinate and
#  the piece extends downward.
#
#    fold_piece($contig, $start, $strand, $len, $true_length)
#
#  Returns a list of [contig, start, strand, len] tuples in 5'->3' order.
#
sub fold_piece
{
	my ($c, $start, $strand, $len, $L) = @_;
	return ([$c, $start, $strand, $len]) unless $len;

	if ($strand eq '-')
	{
		my $b = $start;              # 5' end, higher coordinate
		my $e = $start - $len + 1;   # 3' end, lower coordinate

		#  Wholly in the second copy: shift both ends back by one genome.
		if ($e > $L) { $b -= $L; $e -= $L; }

		#  Straddles the seam: 5' end in the second copy, 3' end in the first.
		#  Runs down from ($b - $L) to 1, then continues from $L down to $e.
		if ($b > $L && $e <= $L)
		{
			return ([$c, $b - $L, '-', $b - $L],
			        [$c, $L,      '-', $L - $e + 1]);
		}
		return ([$c, $b, '-', $b - $e + 1]);
	}

	my $b = $start;                  # 5' end, lower coordinate
	my $e = $start + $len - 1;       # 3' end, higher coordinate

	if ($b > $L) { $b -= $L; $e -= $L; }

	#  Straddles the seam: runs from $b to $L, then continues from 1.
	if ($e > $L && $b <= $L)
	{
		return ([$c, $b, '+', $L - $b + 1],
		        [$c, 1,  '+', $e - $L]);
	}
	return ([$c, $b, '+', $e - $b + 1]);
}
