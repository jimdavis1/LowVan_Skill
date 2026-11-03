#!/bin/bash
#  Adenoviridae dump, 3 Nov 2026. Second DNA family in the project, and an
#  order of magnitude larger per genome than Parvoviridae: 46-74 CDS against
#  5-8.
#
#  --min-len 20000, MEASURED. BV-BRC holds 43,636 Adenoviridae records and
#  only 5,351 reach 5 kb -- ~88% are hexon/fiber typing submissions, the same
#  shape as the parvo capsid-gene holding. Above 5 kb the distribution is
#  nearly flat:
#
#      >=5000  5,351   >=20000 5,015   >=25000 4,975   >=28000 4,898
#
#  25000 looks like the obvious number (a complete adenovirus genome is
#  26-48 kb) and would cut Siadenovirus 66 -> 54. Siadenovirus records run
#  23,811-27,209 nt, so that band is the genus's own size range, not
#  partials. Going down to 20000 costs 40 records family-wide (0.8%) and
#  keeps every genus whole. Per-genus length QC at triage handles the rest.
#
#  Genus holdings at >=20 kb:
#      Mastadenovirus 4,362   Aviadenovirus 450   (blank) 80
#      Siadenovirus 66        Barthadenovirus 55  Ichtadenovirus 2
#  Barthadenovirus is BV-BRC's label for ICTV Atadenovirus (duck/ovine/
#  bovine atadenovirus, snake and lizard adenoviruses). Testadenovirus has
#  219 records and none above 20 kb.
set -eo pipefail
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Retroviridae/env.sh
cd "$HOME/lowvan_scratch/adeno"
python3 "$LOWVAN_KIT/skill/scripts/dump_api.py" \
    --family Adenoviridae --min-len 20000 --prefix A --workdir .
echo "DUMP_DONE $(date)"
