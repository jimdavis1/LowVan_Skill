#!/bin/bash
# Verify a clone of this repository can actually run, before you start a build.
#
# Everything from BV-BRC and SEED is vendored in vendor/bv-brc/, so this checks
# the two things that are not: external binaries, and the handful of CPAN
# modules that cannot be bundled.
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PERL=${LOWVAN_PERL:-perl}
export PERL5LIB="$here/vendor/bv-brc/lib:$here/build${PERL5LIB:+:$PERL5LIB}"
fail=0

echo "perl:   $($PERL -e 'print $^V') ($PERL)"
echo

echo "== binaries =="
for b in blastn tblastn psiblast makeblastdb mmseqs mafft; do
  if command -v "$b" >/dev/null 2>&1; then printf "  ok      %s\n" "$b"
  else printf "  MISSING %s\n" "$b"; fail=1; fi
done
if command -v python3 >/dev/null 2>&1; then
  printf "  ok      python3 (%s)\n" "$(python3 -c 'import sys;print("%d.%d"%sys.version_info[:2])')"
else printf "  MISSING python3\n"; fail=1; fi

echo
echo "== perl modules that must come from CPAN =="
for m in File::Slurp JSON::XS Getopt::Long::Descriptive Class::Accessor IPC::Run Math::Round; do
  if $PERL -M"$m" -e1 >/dev/null 2>&1; then printf "  ok      %s\n" "$m"
  else printf "  MISSING %s\n" "$m"; fail=1; fi
done
if $PERL -MData::UUID -e1 >/dev/null 2>&1 || $PERL -MUUID -e1 >/dev/null 2>&1; then
  printf "  ok      Data::UUID or UUID\n"
else
  printf "  MISSING Data::UUID (or UUID) -- needed to mint feature ids\n"; fail=1
fi

echo
echo "== vendored BV-BRC / SEED modules =="
for m in GenomeTypeObject IDclient P3DataAPI BlastInterface gjoseqlib SeedUtils; do
  if $PERL -M"$m" -e1 >/dev/null 2>&1; then printf "  ok      %s\n" "$m"
  else printf "  BROKEN  %s (should be in vendor/bv-brc/lib)\n" "$m"; fail=1; fi
done

echo
echo "== entry points compile =="
for s in annotate/annotate_by_viral_pssm.pl annotate/annotate_by_viral_pssm-GTO.pl \
         annotate/viral_genome_quality.pl annotate/get_transcript_edited_features.pl \
         annotate/get_splice_variant_features.pl build/fasta-cluster-pssm-2.pl \
         evaluate/New-annotate-viral-taxon.pl; do
  if $PERL -c "$here/$s" >/dev/null 2>&1; then printf "  ok      %s\n" "$s"
  else printf "  FAIL    %s\n" "$s"; fail=1; fi
done

echo
if [ $fail -eq 0 ]; then
  echo "All checks passed."
  echo
  echo "  export PERL5LIB=\"$here/vendor/bv-brc/lib:$here/build:\$PERL5LIB\""
  echo "  export PATH=\"$here/annotate:$here/vendor/bv-brc/bin:\$PATH\""
else
  echo "Some checks failed. For missing CPAN modules:"
  echo "  cpanm File::Slurp Data::UUID JSON::XS Getopt::Long::Descriptive \\"
  echo "        Class::Accessor IPC::Run Math::Round"
  echo "For binaries, conda is easiest:"
  echo "  conda install -c bioconda blast mmseqs2 mafft"
fi
exit $fail
