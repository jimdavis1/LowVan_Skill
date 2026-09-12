#! /usr/bin/env perl
use strict;
use JSON::XS;
use File::Slurp;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Cwd;
use POSIX qw(:sys_wait_h);
use gjoseqlib;
use GenomeTypeObject;

# NOTE: Proc::ParallelLoop is no longer used.  Its pareach() dies outright
# ("unable to fork") the first time the OS refuses a fork, which happens on
# busy machines or under a low `ulimit -u` because every worker here spawns a
# shell plus a multi-program pipeline.  run_parallel() below is a drop-in
# replacement that retries the fork, throttles itself down when forks keep
# failing, and falls back to running the item in the parent as a last resort.


my $usage = 'annotate_viral_taxon.pl [options]

    Annotate viral genomes from BVBRC or from a local contigs directory.

    Modes:
      BVBRC mode    : requires -i Taxon_Name
      Contig dir mode: requires -contigs Contigs_Dir and -meta Metadata_File
                       -i is optional in this mode; if supplied, min/max genome
                       length will be read from the JSON file and used to filter
                       contig files before annotation.
      Download-only mode: -download-only with BVBRC mode.  Downloads the
                       contigs for every genome within [min, max] into the
                       contigs directory, writes a metadata file, and exits
                       without annotating.  If the taxon has no segments entry
                       in the JSON, supply -min/-max explicitly, or the
                       download runs unfiltered.

    General options:
      -h        Help
      -t        Threads (D = 24)
      -j        Full path to the JSON opts file
                (D = /home/jjdavis/bin/Viral_Annotation/Viral_PSSM.json)

    Output directory options:
      -c        Contig GTO folder (D = Contig_GTOs)
      -a        Annotation GTO folder (D = Anno_GTOs)
      -q        Quality GTO folder (D = Quality_GTOs)
      -d        Contigs download folder, BVBRC mode only (D = Contigs)

    BVBRC mode options:
      -i        Taxon name to query BVBRC (also usable in contig dir mode
                to pull min/max lengths from the JSON file)
      -exclude  File of genome IDs to exclude
      -min      Force a min genome length integer (must pair with -max)
      -max      Force a max genome length integer (must pair with -min)

      -download-only  Download the contigs for the genomes within [min, max]
                      into the contigs directory (-d), write a metadata file,
                      and stop.  No GTOs are built and nothing is annotated.
                      Aliases: -download, -dl
      -metaout  Name of the metadata file written by -download-only
                (D = <Taxon_Name>.metadata).  The file has the columns
                required by contig dir mode:
                   file_name TAB genome_name TAB taxon_id
                so a later run can be started with:
                   -contigs <contigs dir> -meta <metaout file>

      NOTE: BVBRC mode filters genomes by the length of a single BVBRC genome
      record.  Segmented viruses whose segments are stored as separate genome
      records in BVBRC will NOT be merged; each segment will be treated as an
      independent genome.  Use contig dir mode (-contigs/-meta) if you need to
      annotate pre-assembled multi-segment genomes as a single unit.

    Contig dir mode options:
      -contigs  Path to the directory of contig files
      -meta     Metadata file with columns: file_name TAB genome_name TAB taxon_id

    Archiving:
      -tar      Tar and gzip the contigs directory when the run finishes.
                Takes an optional file name; with no value the archive is
                named <contigs dir>.tar.gz.  In download-only mode the
                metadata, .processed and .rejected files are added to the
                archive as well.  Examples:
                   -tar
                   -tar Rhabdoviridae_contigs.tar.gz

    Parallelism / fork failures:
      -retries  Times to retry a failed fork before throttling further
                (D = 8).  Each retry waits a little longer than the last.
      Every worker spawns a shell plus a multi-program pipeline, so the real
      process count is several times -t.  If you see fork failures, lower -t
      or raise the shell process limit (ulimit -u).  The run no longer aborts
      on a failed fork: it reduces the number of concurrent workers and keeps
      going, and if even one worker cannot be forked the item is processed in
      the parent process.

    Length filtering (contig dir mode):
      If -min/-max are supplied explicitly, or if -i is supplied and the taxon
      has segments defined in the JSON, contig files whose total sequence length
      falls outside [min, max] are skipped and written to a rejected file.

    Output layout:
      Contig_GTOs/   - one GTO per genome created from raw contigs
      Anno_GTOs/     - annotated GTOs (PSSM + transcript edit + splice variants),
                       plus per-genome .stdout.txt and .stderr.txt from the PSSM step
      Quality_GTOs/  - quality-assessed GTOs plus .feature_quality and .contig_quality files

    Notes:
      The annotation pipeline pipes annotate_by_viral_pssm-GTO.pl into
      get_transcript_edited_features.pl and get_splice_variant_features.pl.
      Both downstream programs skip their processing automatically when the
      viral family has no transcript-edited or splice-variant features.

';


my $contig_gto  = "Contig_GTOs";
my $anno_gto    = "Anno_GTOs";
my $quality_gto = "Quality_GTOs";
my $contigD     = "Contigs";
my $threads     = 24;
my $retries     = 8;
my $json        = "/home/jjdavis/bin/Viral_Annotation/Viral_PSSM.json";
my $base        = getcwd;

my ($help, $taxon, $min, $max, $exF, $contigDir, $metaF, $download_only, $metaOut, $tar);

my $opts = GetOptions(
	'h'         => \$help,
	'i=s'       => \$taxon,
	't=i'       => \$threads,
	'j=s'       => \$json,
	'a=s'       => \$anno_gto,
	'q=s'       => \$quality_gto,
	'c=s'       => \$contig_gto,
	'd=s'       => \$contigD,
	'contigs=s' => \$contigDir,
	'meta=s'    => \$metaF,
	'min=i'     => \$min,
	'max=i'     => \$max,
	'exclude=s' => \$exF,
	'download-only|download|dl' => \$download_only,
	'metaout=s' => \$metaOut,
	'tar:s'     => \$tar,
	'retries=i' => \$retries,
);

if ($help)
{
	die "$usage\n";
}

# Validate mode
if ($contigDir || $metaF)
{
	unless ($contigDir && $metaF)
	{
		die "Contig dir mode requires both -contigs and -meta\n\n$usage\n";
	}
}
else
{
	unless ($taxon)
	{
		die "Must supply -i Taxon_Name for BVBRC mode, or use -contigs/-meta for contig dir mode\n\n$usage\n";
	}
}

if (($min && !$max) || ($max && !$min))
{
	die "Must declare both -min and -max\n\n$usage\n";
}

# Download-only is a BVBRC-mode operation; there is nothing to download in
# contig dir mode.
if ($download_only && $contigDir)
{
	die "-download-only cannot be used with -contigs/-meta (contig dir mode)\n\n$usage\n";
}

$threads = 1 if $threads < 1;
$retries = 0 if $retries < 0;


# ============================================================================
# Fork-tolerant parallel loop
#
#   run_parallel(\@items, sub { ... }, $max_workers)
#
# Replaces Proc::ParallelLoop::pareach.  Differences that matter here:
#
#   * A fork failure is retried with an increasing delay instead of killing
#     the run.  Finished children are reaped first, since the usual cause is
#     simply that every process slot is in use.
#   * After repeated failures the worker cap is lowered, so a machine that
#     cannot sustain 24 workers settles at whatever it can sustain rather
#     than failing over and over.
#   * If a fork still cannot be had, the item runs in the parent process.
#     Slower, but the item is not silently skipped.
#   * Children exit through POSIX::_exit so they never flush a copy of the
#     parent's buffered output or run its cleanup a second time.
# ============================================================================

my %KIDS;   # pid => item, for children of the loop currently running


sub reap_kids
{
	my ($blocking) = @_;

	my $reaped = 0;

	while (1)
	{
		my $pid = waitpid(-1, $blocking ? 0 : WNOHANG);
		last if $pid <= 0;

		my $status = $?;
		my $item   = delete $KIDS{$pid};

		if ($status && defined $item)
		{
			warn "Worker for '$item' exited with status " . ($status >> 8) . "\n";
		}

		$reaped++;
		last if $blocking;
	}

	return $reaped;
}


sub run_parallel
{
	my ($items, $code, $workers) = @_;

	$workers = 1 if !$workers || $workers < 1;

	my $cap     = $workers;   # what the user asked for
	my $streak  = 0;          # consecutive successful forks since a failure

	%KIDS = ();

	foreach my $item (@$items)
	{
		# Wait for a free slot
		while (scalar(keys %KIDS) >= $workers)
		{
			reap_kids(1);
		}

		my $pid;
		my $delay = 1;

		for (my $try = 0; $try <= $retries; $try++)
		{
			$pid = fork();
			last if defined $pid;

			# Fork refused.  Reap anything that has finished, give the
			# machine a moment, and shrink the pool.
			warn "fork failed ($!)"
			   . ($try < $retries ? "; retrying in ${delay}s" : "")
			   . " [$workers workers, " . scalar(keys %KIDS) . " running]\n";

			reap_kids(0);

			if (scalar(keys %KIDS))
			{
				reap_kids(1);
			}

			$streak = 0;

			if ($workers > 1)
			{
				$workers--;
				warn "Reducing concurrency to $workers workers\n";
			}

			last if $try == $retries;

			sleep $delay;
			$delay *= 2;
			$delay = 30 if $delay > 30;
		}

		unless (defined $pid)
		{
			# Could not fork at all: do the work here so nothing is lost
			warn "Could not fork for '$item'; running it in the parent process\n";
			eval { $code->($item); 1 } or warn "Error processing '$item': $@";
			next;
		}

		if ($pid == 0)
		{
			# Child
			eval { $code->($item); 1 } or warn "Error processing '$item': $@";
			POSIX::_exit(0);
		}

		$KIDS{$pid} = $item;

		# The machine is behaving again: creep back up toward the requested
		# number of workers rather than staying throttled for the whole run.
		if (++$streak >= 20 && $workers < $cap)
		{
			$workers++;
			$streak = 0;
			warn "Raising concurrency back to $workers workers\n";
		}
	}

	# Drain
	while (scalar(keys %KIDS))
	{
		reap_kids(1);
	}
}


# ----------------------------------------------------------------------------
# Tar up a contigs directory, optionally with a few loose files alongside it.
# Does nothing unless -tar was given.
# ----------------------------------------------------------------------------
sub tar_contigs
{
	my ($dir, @extra) = @_;

	return unless defined $tar;

	unless (-d $dir)
	{
		warn "Nothing to tar: $dir is not a directory\n";
		return;
	}

	my $parent = dirname($dir);
	my $name   = basename($dir);

	my $archive = ($tar ne "") ? $tar : "$name.tar.gz";
	$archive .= ".tar.gz" unless $archive =~ /\.(?:tar\.gz|tgz|tar)$/;

	my $cmd = "tar -czf \"$archive\" -C \"$parent\" \"$name\"";

	foreach my $f (@extra)
	{
		next unless defined $f && -f $f;
		$cmd .= " -C \"" . dirname($f) . "\" \"" . basename($f) . "\"";
	}

	my $rc = system $cmd;

	if ($rc == 0)
	{
		print STDERR "Wrote $archive\n";
	}
	else
	{
		warn "tar failed (exit " . ($rc >> 8) . "): $cmd\n";
	}
}


# Load the exclude hash if provided
my %exclude;
if ($exF)
{
	open(IN, "<$exF") or die "Could not open exclude file: $exF\n";
	while (<IN>)
	{
		chomp;
		$exclude{$_}++;
	}
	close IN;
}


# Read the JSON options file
open(IN, "<$json") or die "Cannot find JSON options file: $json\n";
my $options = decode_json(scalar read_file(\*IN));
close IN;


# Determine min/max genome length from JSON when not supplied on the command line.
# This applies to both modes when -i is provided.
unless ($min && $max)
{
	if ($taxon)
	{
		if (exists $options->{$taxon}->{segments})
		{
			$max = 0;
			$min = 0;
			foreach my $segment (values %{$options->{$taxon}->{"segments"}})
			{
				$max += $segment->{"max_len"};
				$min += $segment->{"min_len"};
			}
		}
		elsif (exists $options->{$taxon}->{min_len} && exists $options->{$taxon}->{max_len})
		{
			# Non-segmented taxon that carries its lengths at the top level
			$min = $options->{$taxon}->{min_len};
			$max = $options->{$taxon}->{max_len};
		}
		elsif ($download_only)
		{
			# Nothing to derive lengths from.  Downloading is harmless without
			# a filter, so warn and pull everything in the taxon.
			warn "No segments entry for '$taxon' in JSON; " .
			     "downloading all genomes in the taxon (use -min/-max to filter)\n";
		}
		else
		{
			die "No segments entry for '$taxon' in JSON\n";
		}
	}
}


# Create output directories.
# In download-only mode no GTOs are produced, so these are not made.
unless ($download_only)
{
	mkdir $contig_gto;
	mkdir $anno_gto;
	mkdir $quality_gto;
}


# Collect the list of items to process and their metadata.
# In BVBRC mode  : keys are genome IDs, metadata comes from the BVBRC query.
# In contig mode : keys are filenames, metadata comes from the -meta file.

my @items;          # genome IDs (BVBRC) or contig filenames (contig dir)
my %genome_meta;    # item => { NAME => ..., TAX => ... }
my ($goodF, $badF); # names of the .processed / .rejected files


if ($contigDir)
{
	# -------------------------------------------------------------------------
	# Contig directory mode: read metadata file, then optionally filter by
	# total contig length using gjoseqlib::read_fasta.
	# -------------------------------------------------------------------------

	my $label = $taxon // "contigs";
	$goodF = "$label.processed";
	$badF  = "$label.rejected";

	open(GOOD, ">$goodF") or die "Cannot open $goodF\n";
	open(BAD,  ">$badF")  or die "Cannot open $badF\n";

	open(IN, "<$metaF") or die "Could not open metadata file: $metaF\n";
	while (<IN>)
	{
		chomp;
		my ($file, $name, $tax) = split /\t/;

		next if exists $exclude{$file};

		# Sum the lengths of all contigs in the file
		my @seqs     = gjoseqlib::read_fasta("$contigDir/$file");
		my $total_len = 0;
		foreach my $seq (@seqs)
		{
			$total_len += length($seq->[2]);
		}

		# Apply length filter if min/max are defined
		if ($min && $max && (($total_len < $min) || ($total_len > $max)))
		{
			print BAD "$file\t$total_len\t$name\n";
			next;
		}

		push @items, $file;
		$genome_meta{$file}->{NAME} = $name;
		$genome_meta{$file}->{TAX}  = $tax;
		print GOOD "$file\t$total_len\t$name\n";
	}
	close IN;
	close GOOD;
	close BAD;
}
else
{
	# -------------------------------------------------------------------------
	# BVBRC download mode: query BVBRC for genomes of the right length
	# -------------------------------------------------------------------------
	mkdir $contigD;

	$goodF = "$taxon.processed";
	$badF  = "$taxon.rejected";

	open(IN,   "echo $taxon | query_PATRIC_bob.pl -c genome -i taxon_lineage_names -r \"genome_id genome_length genome_name\" | ");
	open(GOOD, ">$goodF") or die "Cannot open $goodF\n";
	open(BAD,  ">$badF")  or die "Cannot open $badF\n";

	while (<IN>)
	{
		chomp;
		my ($id, $len, $name) = split /\t/;
		print "$_\n";

		unless (exists $exclude{$id})
		{
			# With no min/max defined (download-only on a taxon with no
			# lengths in the JSON) every genome passes.
			if (!($min && $max) || (($len < $max) && ($len > $min)))
			{
				push @items, $id;
				$genome_meta{$id}->{NAME} = $name;
				$genome_meta{$id}->{TAX}  = $id;   # taxon derived from genome_id below
				print GOOD "$id\t$len\t$name\n";
			}
			else
			{
				print BAD "$id\t$len\t$name\n";
			}
		}
	}
	close IN;
	close GOOD;
	close BAD;
}


# Keep the parent's buffers from being duplicated into the children
$| = 1;


# ============================================================================
# Download-only mode
#
# Downloads the contigs for every genome that passed the [min, max] filter
# into $contigD, writes a metadata file that can be fed straight back in
# with -contigs/-meta, and exits before any GTO or annotation step.
# ============================================================================

if ($download_only)
{
	$metaOut //= "$taxon.metadata";

	run_parallel(
		\@items,

		sub
		{
			my $item = shift @_;
			system "echo $item | BVBRC_clean_contigs.pl -d $base/$contigD";
		},

		$threads
	);

	open(META, ">$metaOut") or die "Cannot open $metaOut\n";

	my $got     = 0;
	my $missing = 0;

	foreach my $item (@items)
	{
		my $file = "$item.contigs";

		unless (-s "$base/$contigD/$file")
		{
			warn "No contigs downloaded for $item\n";
			$missing++;
			next;
		}

		my $tax = $item;
		$tax =~ s/\..+//g;   # strip version suffix to get integer taxon ID

		print META join("\t", $file, $genome_meta{$item}->{NAME}, $tax), "\n";
		$got++;
	}
	close META;

	print STDERR "Downloaded $got genomes into $contigD\n";
	print STDERR "Failed to download $missing genomes\n" if $missing;
	print STDERR "Metadata written to $metaOut\n";
	print STDERR "To annotate: $0 -contigs $contigD -meta $metaOut" .
	             ($taxon ? " -i \"$taxon\"" : "") . "\n";

	# Archive the contigs directory plus the run's bookkeeping files
	tar_contigs($contigD, $metaOut, $goodF, $badF);

	exit 0;
}


# ============================================================================
# Main parallel processing loop
# ============================================================================

run_parallel(
	\@items,

	sub
	{
		my $item = shift @_;

		my $name = $genome_meta{$item}->{NAME};
		my $tax  = $genome_meta{$item}->{TAX};
		$tax =~ s/\..+//g;   # strip version suffix to get integer taxon ID

		# Derive a clean ID: strip common contig-file extensions if present
		my $id = $item;
		$id =~ s/\.(?:dna|fasta|fa|fna|contigs?)$//;


		# ------------------------------------------------------------------
		# Step 1: Build the contig GTO from raw sequences
		# ------------------------------------------------------------------
		if ($contigDir)
		{
			system "rast-create-genome "
			     . "--ncbi-taxonomy-id $tax "
			     . "--scientific-name \"$name\" "
			     . "--domain Viruses "
			     . "--genetic-code 11 "
			     . "--contigs $base/$contigDir/$item "
			     . "> $base/$contig_gto/$id.contig.gto";
		}
		else
		{
			system "echo $item | BVBRC_clean_contigs.pl -d $base/$contigD";
			system "rast-create-genome "
			     . "--ncbi-taxonomy-id $tax "
			     . "--scientific-name \"$name\" "
			     . "--domain Viruses "
			     . "--genetic-code 11 "
			     . "--contigs $base/$contigD/$item.contigs "
			     . "> $base/$contig_gto/$id.contig.gto";
		}


		# ------------------------------------------------------------------
		# Step 2: Annotate and post-process
		#
		# Pipeline:
		#   annotate_by_viral_pssm-GTO.pl  (PSSM-based feature calls)
		#     -> get_transcript_edited_features.pl  (skips if family has none)
		#     -> get_splice_variant_features.pl     (skips if family has none)
		#
		# The final GTO from this pipeline is the Anno GTO.
		# annotate_by_viral_pssm-GTO.pl writes $id.stdout.txt and
		# $id.stderr.txt to the working directory; we move those into
		# Anno_GTOs for traceability.
		# ------------------------------------------------------------------

		system "annotate_by_viral_pssm-GTO.pl -x $id "
		     . "< $base/$contig_gto/$id.contig.gto "
		     . "| get_transcript_edited_features.pl "
		     . "| get_splice_variant_features.pl "
		     . "> $base/$anno_gto/$id.anno.gto";

		system "mv $base/$id.stdout.txt $base/$anno_gto/" if -f "$base/$id.stdout.txt";
		system "mv $base/$id.stderr.txt $base/$anno_gto/" if -f "$base/$id.stderr.txt";


		# ------------------------------------------------------------------
		# Step 3: Quality assessment
		#
		# Reads the finished Anno GTO and writes a Quality GTO plus two
		# tabular report files that we move into Quality_GTOs.
		# ------------------------------------------------------------------

		system "viral_genome_quality.pl "
		     . "-p $id "
		     . "-i $base/$anno_gto/$id.anno.gto "
		     . "-o $base/$quality_gto/$id.qual.gto";

		system "mv $base/$id.feature_quality $base/$quality_gto/" if -f "$base/$id.feature_quality";
		system "mv $base/$id.contig_quality  $base/$quality_gto/" if -f "$base/$id.contig_quality";
	},

	$threads
);


# Archive the contigs directory that was used for this run
tar_contigs($contigDir ? $contigDir : $contigD);
