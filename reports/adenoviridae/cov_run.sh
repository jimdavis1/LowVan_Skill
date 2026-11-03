#!/bin/bash
#  Family-wide coverage for Adenoviridae. STRICTLY SERIAL and strictly after
#  the panel run: three kernel panics on Parvoviridae (28 and 30 Oct) traced
#  to IPC voucher exhaustion, ~515 vouchers per genome when routing ran one
#  blastn per reference against 1,271 reference files. VA_adeno holds 150, so
#  the per-genome cost is ~8x lower, but two heavy runs at once is how the
#  ceiling was reached the first time. The vouchers RECLAIM when runs stop --
#  a panic does not need a reboot to recover from, it needs the runs stopped.
#
#  --min-complete 20000 is the dump's own floor (set 3 Nov after measuring
#  that 25000 would cut Siadenovirus 66 -> 54, and that band is the genus's
#  real size range). Below it a record is a partial submission, which is
#  flagged poor for missing essential features -- correct, and uninformative
#  about the module.
#
#  Specials are NOT skipped. A coverage audit that omits a module's spliced or
#  transcript-edited features understates it and says nothing about having
#  done so.
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Retroviridae/env.sh
P=$HOME/lowvan_scratch/adeno; S=$LOWVAN_KIT/skill/scripts
REPO=$HOME/lowvan_scratch/VA_adeno
export LOWVAN_DATA_DIR=$REPO
mkdir -p $P/cov
python3 $S/evaluate_coverage.py --contigs $P/Contigs --repo $REPO \
        --work $P/cov/work --min-complete 20000 --threads ${THREADS:-4} \
        --taxon Adenoviridae \
        --out $P/cov/coverage.json > $P/logs/coverage.log 2>&1
echo "coverage exit=$?"
