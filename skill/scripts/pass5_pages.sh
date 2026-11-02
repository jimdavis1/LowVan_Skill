#!/bin/bash
#  Everything downstream of the coverage run, for one module or for all nine.
#  Usage: bash pass5_pages.sh [Module ...]
#
#  Order matters and is not arbitrary:
#    copy_num must be applied BEFORE the module JSON is regenerated, because
#    the JSON carries it; the vocabulary TSV must be written AFTER that,
#    because its "Used for Genome Quality" column is derived from copy_num
#    and reads 0 for everything until then.
#
#  env.sh already provides mmseqs, blast and the BV-BRC perl, in that order.
#  Do not prepend the lowvan conda env to PATH: its perl has no
#  GenomeTypeObject.pm and every genome then fails to annotate while the
#  runner still exits 0.
set -uo pipefail
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Retroviridae/env.sh
P=$HOME/lowvan_scratch/parvo
S=$LOWVAN_KIT/skill/scripts
REPO=$HOME/lowvan_scratch/VA_spina
#  ROUTE can be overridden to a module-only reference repo built by
#  make_module_repo.sh. The rescore reuses the existing .ann.gto, so the
#  reference set barely matters here -- but evaluate_coverage re-annotates
#  its exemplars and does route them, so point it at the same repo the
#  coverage run used or the two measurements disagree.
ROUTE=${ROUTE:-$HOME/lowvan_scratch/VA_parvo_route}
A=$P/artifacts
MODS="${@:-Dependoparvovirus Protoparvovirus Bocaparvovirus Erythroparvovirus Amdoparvovirus TetraCopiparvovirus Hamaparvovirinae DensovirinaeAmbi DensovirinaeMono}"
mkdir -p "$P/coverage" "$A/rare" "$A/notes" "$A/generators"

perl -MGenomeTypeObject -e1 2>/dev/null || {
  echo "FATAL: perl on PATH has no GenomeTypeObject.pm -- $(command -v perl)"; exit 1; }

disp () {  # display name, per the convention the string-collapse pages set
  case $1 in
    TetraCopiparvovirus) echo "TetraCopiparvovirus — Tetraparvovirus and Copiparvovirus";;
    DensovirinaeAmbi)    echo "DensovirinaeAmbi — ambisense densoviruses";;
    DensovirinaeMono)    echo "DensovirinaeMono — monosense densoviruses";;
    *) echo "$1";;
  esac
}

for m in $MODS; do
  echo "=============== $m"
  c=$P/cov/$m
  n=$(ls $c/tbl/*.feature.tbl 2>/dev/null | wc -l | tr -d ' ')
  [ "$n" -gt 0 ] || { echo "  no coverage tables yet, skipped"; continue; }
  echo "  $n coverage tables"

  #  1. copy_num / essentiality, from the .qual.gto the coverage run wrote.
  #     Retroviridae ENV is a curator override elsewhere; nothing here needs
  #     --require, and spliced features are never essential.
  #     apply_copy_num rewrites the MODULE JSON (it is keyed by module), not
  #     features_all.json -- the first run of this script passed the latter and
  #     died on KeyError. Do NOT regenerate the module JSON afterwards: it is
  #     built from features_all.json, which carries no copy_num, so a
  #     regeneration silently strips what this step just wrote.
  python3 "$S/apply_copy_num.py" --eval "$c/anno" --json "$P/$m/${m}_Viral_PSSM.json" \
          --module "$m" --write > "$P/logs/$m.copynum.log" 2>&1
  echo "  copy_num exit=$? :: $(grep -c 'copy_num' $P/logs/$m.copynum.log 2>/dev/null) line(s)"

  #  2. reinstall (VA_parvo_route symlinks Viral_PSSM.json to $REPO), then
  #     RE-SCORE quality: the .qual.gto files were scored with no feature
  #     essential, so report.tsv -- the coverage audit's panel -- would not
  #     reflect copy_num. Annotations are reused; run_gto_eval skips the
  #     special passes on a reused .ann.gto, and qc_duplicate_calls checks it.
  python3 "$S/install_module.py" --workdir "$P/$m" --repo "$REPO" \
          > "$P/logs/$m.install5.log" 2>&1
  echo "  install exit=$?"
  rm -f "$c"/anno/*.qual.gto
  python3 "$S/run_gto_eval.py" --gto-dir "$c/gto" --repo "$ROUTE" \
          --out "$c/anno" --report "$c/report.tsv" --tbl-dir "$c/tbl" \
          --jobs 4 > "$P/logs/$m.rescore.log" 2>&1
  echo "  rescore exit=$? :: $(grep -E 'good \(' $P/logs/$m.rescore.log | xargs)"
  python3 "$S/qc_duplicate_calls.py" --eval "$c/anno" --json "$P/$m/${m}_Viral_PSSM.json" \
          --module "$m" > "$P/logs/$m.dupcalls.log" 2>&1
  #  grep -c EXACT also matches the "EXACT -- none" header, so a clean module
  #  reported "1 EXACT line(s)" and looked like it had a duplicate. Count the
  #  rows UNDER that header instead.
  echo "  duplicate-call qc exit=$? :: $(awk '/^  EXACT/{f=1;next} /^$/{f=0} f&&/[A-Za-z]/{n++} END{print n+0}' $P/logs/$m.dupcalls.log) EXACT pair(s)"

  #  3. registry again (cutoffs and copy_num both show on the page).
  python3 "$S/collect_registry.py" --workdir "$P/$m" --out "$P/$m/registry.json" \
          > "$P/logs/$m.registry5.log" 2>&1
  cp "$P/$m/registry.json" "$A/registry/$m.json"
  echo "  registry exit=$?"

  #  4. vocabulary growth: what the module calls now versus what was deposited.
  #     The deposited side is joined back through P.id_md5 -- uniq.id_ann holds
  #     one row per UNIQUE protein, so reading it directly would give each
  #     genome only its representatives' vocabulary, which is the very thing
  #     being measured.
  python3 "$S/annotation_rarefaction.py" --new "$c/tbl" --old "$A/rare/$m.old.tsv" \
          --out "$P/$m/rarefaction.json" > "$P/logs/$m.rare.log" 2>&1
  echo "  rarefaction exit=$?"

  #  5. coverage JSON -- exemplar-based, a far smaller run than the 2,200-genome
  #     sweep, and it must NOT be given --no-special where specials exist.
  python3 "$S/evaluate_coverage.py" --contigs "$c/contigs" --repo "$ROUTE" \
          --work "$P/coverage/_wd_$m" --json "$P/$m/${m}_Viral_PSSM.json" \
          --taxon "$m" --out "$P/coverage/$m.json" \
          > "$P/logs/$m.evalcov.log" 2>&1
  echo "  evaluate_coverage exit=$? -> $(ls -l $P/coverage/$m.json 2>/dev/null | awk '{print $5}') bytes"

  #  6. the two pages that needed coverage.
  s=$(echo "$m" | tr 'A-Z' 'a-z')
  python3 "$S/gen_rarefaction.py" --rarefaction "$P/$m/rarefaction.json" \
          --taxon "$m" --display "$(disp $m)" \
          --out "$A/$s-vocabulary-saturation.html" > "$P/logs/$m.satpage.log" 2>&1
  echo "  saturation page exit=$?"

  if [ -f "$A/generators/${s}_prose.json" ]; then
    #  --meta is NOT optional here. Without it gen_facts falls back to
    #  Picorna.full.tsv, which does not exist outside the picornavirus build,
    #  and the audit page dies with FileNotFoundError while every other step
    #  reports success -- the module looks finished and is missing a page.
    python3 "$S/gen_facts.py" "$m" "$A/generators/${s}_prose.json" \
            "$A/generators/${s}_audit.json" --scratch "$P" --workdir "$P/$m" \
            --coverage "$P/coverage/$m.json" --contigs "$c/contigs" \
            --panel "$c/report.tsv" --family Parvoviridae \
            --meta "$P/$m/G.meta.tsv" \
            > "$P/logs/$m.facts.log" 2>&1 \
      && python3 "$S/gen_coverage_audit.py" --facts "$A/generators/${s}_audit.json" \
            --out "$A/$s-coverage-audit.html" > "$P/logs/$m.auditpage.log" 2>&1
    echo "  coverage audit exit=$?"
  else
    echo "  coverage audit SKIPPED -- no $A/generators/${s}_prose.json"
  fi
done

echo "=============== pages"
python3 "$S/check_artifact_html.py" $A/*.html 2>&1 | tail -6
ls "$A"/*.html | wc -l
