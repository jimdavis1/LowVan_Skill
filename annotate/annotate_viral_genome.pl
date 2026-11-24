#! /usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path remove_tree);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin;
use Getopt::Long::Descriptive;

my $tool_version;
eval {
    require LowVanVersion;
    $tool_version = LowVanVersion::get_version();
};
if ($@ || !$tool_version) {
    $tool_version = "dev";
}

my $program_description = <<'END_DESCRIPTION';
Annotate ONE viral genome, end to end, and leave a GTO that the quality checker
can score.  In goes a contig GTO, out comes a finished annotation GTO; with -q
it also writes the quality GTO.

    annotate_viral_genome.pl -i 10407.1002.contig.gto -o 10407.1002.anno.gto \
                             -q 10407.1002.qual.gto -x 10407.1002

WHY THIS EXISTS.  The annotation of a viral genome is four programs in a fixed
order, and until now nothing in this kit ran all four.  Callers assembled the
chain themselves, each one slightly differently, and every one of them stopped
after three:

    annotate_by_viral_pssm-GTO.pl | get_transcript_edited_features.pl
                                  | get_splice_variant_features.pl

That chain is correct and incomplete.  The missing fourth stage is
fold_circular_features.pl, and leaving it out is not a cosmetic omission -- it
produces a finished-looking GTO that is quantitatively wrong.

THE DUPLICATION EVENT, and what has to reconcile it.  A circular replicon
reaches a database as the loop cut at one arbitrary point, so a gene lying
across that cut arrives in two pieces and is found by neither.  On hepatitis B
that is not an edge case: the polymerase crosses the origin in 97.2 per cent of
near-full-length deposits and the large envelope protein in all of them.  So
annotate_by_viral_pssm.pl DOUBLES a contig whose module declares
`replicon_geometry: circular`, searches the doubled copy where every gene
appears whole and in order, and the GTO wrapper doubles the contig in the GTO to
match so that the special-feature callers inherit the same view.

Everything after that point is then looking at a genome that does not exist.
Each real feature has been found twice, and the contig is twice its deposited
length, and viral_genome_quality.pl -- which has no way to know any of this --
measures what it is given.  The result is not a warning, it is a verdict:

    Contig  Len   Min_Len  Max_len  Exceptions
    ..._1   6622  3130     3400     Too many HSPs for: RNA-dependent DNA
                                    polymerase;Too many HSPs for: Large surface
                                    antigen;...;Contig is too long

Eight of eight Orthohepadnavirus exemplars scored `Poor` that way, with every
feature flagged, on a module whose features were in fact all called correctly.
Measured over the full exemplar set it read `0 good, 589 poor`.  A module that
works perfectly and a module that is broken produce the same report, which makes
the quality score useless on exactly the taxa that need it most.

fold_circular_features.pl reconciles the duplication: it folds coordinates past
the true length into two-part join locations, drops the second copy of
everything, drops the leftover partial that the whole call swallows, and
restores each contig to its deposited sequence.  It must run BEFORE
viral_genome_quality.pl and AFTER the special-feature callers -- after, because
a spliced product crossing the origin can only be found on the doubled view.
That is a one-place-only slot in the order, which is the argument for this
script existing rather than for documenting the chain and hoping.

STAGES, in order.  Each one is skipped cleanly when it has nothing to do, so
this is the right driver for a linear taxon too:

    1. annotate_by_viral_pssm-GTO.pl      profile-based CDS/mat_peptide calling
    2. get_transcript_edited_features.pl  features needing an RNA edit
    3. get_splice_variant_features.pl     features called from spliced references
    4. fold_circular_features.pl          undo the doubling (no-op if not doubled)
    5. viral_genome_quality.pl            only with -q

WHY STAGED FILES AND NOT LITERAL PIPES.  Every stage reads its whole GTO through
GenomeTypeObject->create_from_file and writes it whole, so nothing streams and a
`|` buys no concurrency.  What it costs is diagnosis: in a shell pipeline the
exit status of an inner stage is invisible by default, a stage that dies mid-way
emits a truncated GTO that the next stage parses as best it can, and the failure
surfaces several stages later as missing features.  Writing each stage to its own
file in a work directory means a failure is attributed to the stage that caused
it, and --keep-work leaves every intermediate for inspection.

RESUMING.  -from fold takes an ALREADY ANNOTATED GTO and runs only the fold and
the quality step.  That is the repair path for a GTO produced by some other
driver that stopped at stage 3, and it is safe on a GTO that was already folded:
the fold clears circular_true_length when it restores a contig, so a second pass
is a no-op.
END_DESCRIPTION

my ($opt, $usage) = describe_options(
    "\n$program_description\nUsage: %c %o",
    ["input|i=s"      => "Input contig GTO (default: STDIN)"],
    ["output|o=s"     => "Output annotation GTO (default: STDOUT)"],
    ["quality|q=s"    => "Also write a quality GTO here (runs viral_genome_quality.pl)"],
    ["prefix|x=s"     => "Genome/file prefix, passed to the annotator and the quality step",
                         { default => "Viral_Anno" }],
    ["json|j=s"       => "Full path to the JSON opts file; passed to every stage that takes one"],
    ["threads|t=i"    => "Threads for the annotator", { default => 8 }],
    ["from|f=s"       => "First stage to run: annotate (default) or fold",
                         { default => "annotate" }],
    ["work|w=s"       => "Work directory for the intermediates (default: a temp dir)"],
    ["keep-work|k"    => "Keep the work directory and every intermediate GTO"],
    ["quiet"          => "Suppress the per-stage progress lines on STDERR"],
    ["help|h"         => "Show this help message", { shortcircuit => 1 }],
    ["version|v"      => "Show version information", { shortcircuit => 1 }],
);

print($usage->text), exit 0 if $opt->help;
if ($opt->version) { print "annotate_viral_genome.pl version $tool_version\n"; exit 0; }

my %valid_from = (annotate => 1, fold => 1);
$valid_from{$opt->from}
    or die "-from must be 'annotate' or 'fold', not '" . $opt->from . "'\n";

#  SIBLINGS, NOT $PATH.  The four stages are versioned together with this
#  script and a mixed set is a silent wrong answer, not an error -- an older
#  fold_circular_features.pl on $PATH would fold with the pre-containment rules
#  and leave the duplicate partial behind. So each stage is resolved next to
#  this file, and only falls back to $PATH if it is genuinely not there (a
#  deployment that installs the stages elsewhere).
sub stage_path
{
	my ($name) = @_;
	my $sib = File::Spec->catfile($FindBin::RealBin, $name);
	return $sib if -e $sib;
	return $name;
}

my $work = $opt->work;
my $work_is_temp = 0;
if ($work) { make_path($work) unless -d $work; }
else       { $work = tempdir("annotate_viral_genome.XXXXXX", TMPDIR => 1); $work_is_temp = 1; }

my $prefix = $opt->prefix;
my $step   = 0;

sub note
{
	return if $opt->quiet;
	print STDERR "[annotate_viral_genome] @_\n";
}

#  Run one stage, from $in to a fresh file in $work, and return that file.
#
#  A stage is only allowed to succeed if it exits zero AND leaves a non-empty
#  output. Both halves are needed: a stage that dies after opening its output
#  leaves a zero-byte GTO behind, and a stage whose perl dependencies are
#  missing exits non-zero having written nothing. Neither must be passed on as
#  if it were an annotation.
sub run_stage
{
	my ($script, $in, $label, @args) = @_;
	$step++;
	my $out = File::Spec->catfile($work, sprintf("%s.%d.%s.gto", $prefix, $step, $label));
	my @cmd = ($^X, stage_path($script), "-i", $in, "-o", $out, @args);
	note("$label: $script");
	my $rc = system(@cmd);
	if ($rc != 0)
	{
		my $why = $rc == -1 ? "could not be run: $!"
		        : ($rc & 127) ? sprintf("died on signal %d", $rc & 127)
		        : sprintf("exited %d", $rc >> 8);
		die "annotate_viral_genome: stage '$label' ($script) $why.\n"
		  . "Intermediates are in $work\n";
	}
	unless (-s $out)
	{
		die "annotate_viral_genome: stage '$label' ($script) exited 0 but wrote "
		  . "nothing to $out.\nIntermediates are in $work\n";
	}
	return $out;
}

#  STDIN has to become a file, because every stage is given -i explicitly and
#  because an intermediate that can be re-read is the whole point of staging.
my $cur = $opt->input;
unless (defined $cur)
{
	$cur = File::Spec->catfile($work, "$prefix.0.input.gto");
	open(my $IN, ">", $cur) or die "Cannot write $cur: $!\n";
	binmode $IN;
	local $/;
	my $all = <STDIN>;
	defined $all && length $all or die "annotate_viral_genome: no GTO on STDIN\n";
	print $IN $all;
	close $IN;
}
-s $cur or die "annotate_viral_genome: input GTO $cur is missing or empty\n";

my @json = $opt->json ? ("-j", $opt->json) : ();

if ($opt->from eq 'annotate')
{
	$cur = run_stage("annotate_by_viral_pssm-GTO.pl", $cur, "annotate",
	                 "-x", $prefix, "-t", $opt->threads, @json);
	$cur = run_stage("get_transcript_edited_features.pl", $cur, "edited", @json);
	$cur = run_stage("get_splice_variant_features.pl",   $cur, "spliced", @json);
}
else
{
	note("starting at the fold; the input is taken to be already annotated");
}

#  THE STAGE EVERY OTHER CALLER FORGOT. Unconditional: it reads
#  circular_true_length off the contigs and passes a genome that was never
#  doubled through untouched, so there is no geometry to test here and no way
#  for a caller to get the test wrong.
$cur = run_stage("fold_circular_features.pl", $cur, "folded");

if (defined $opt->output)
{
	open(my $SRC, "<", $cur)           or die "Cannot read $cur: $!\n";
	open(my $DST, ">", $opt->output)   or die "Cannot write " . $opt->output . ": $!\n";
	binmode $SRC; binmode $DST;
	print $DST $_ while <$SRC>;
	close $SRC; close $DST;
}
else
{
	open(my $SRC, "<", $cur) or die "Cannot read $cur: $!\n";
	binmode $SRC; binmode STDOUT;
	print STDOUT $_ while <$SRC>;
	close $SRC;
}

#  QUALITY LAST, AND ONLY ON THE FOLDED GENOME.
#
#  viral_genome_quality.pl also drops two tabular reports, $prefix.feature_quality
#  and $prefix.contig_quality, into its working directory. They are part of the
#  output, not debris, so the step is run beside the quality GTO and they land
#  next to it instead of in whatever directory the caller happened to be in.
if ($opt->quality)
{
	$step++;
	my $qdir = dirname(File::Spec->rel2abs($opt->quality));
	make_path($qdir) unless -d $qdir;
	my @cmd = ($^X, stage_path("viral_genome_quality.pl"),
	           "-i", File::Spec->rel2abs($cur),
	           "-o", File::Spec->rel2abs($opt->quality),
	           "-p", $prefix, @json);
	note("quality: viral_genome_quality.pl");
	my $cwd = `pwd`; chomp $cwd;
	chdir $qdir or die "Cannot chdir to $qdir: $!\n";
	my $rc = system(@cmd);
	chdir $cwd;
	$rc == 0 or die "annotate_viral_genome: viral_genome_quality.pl failed (rc=$rc).\n"
	              . "Intermediates are in $work\n";
	-s $opt->quality
	    or die "annotate_viral_genome: viral_genome_quality.pl wrote no quality GTO.\n"
	         . "Intermediates are in $work\n";
}

if ($opt->keep_work) { note("intermediates kept in $work"); }
elsif ($work_is_temp) { remove_tree($work); }

exit 0;
