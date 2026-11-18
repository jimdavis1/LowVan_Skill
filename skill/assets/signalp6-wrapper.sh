#!/bin/bash
#  SignalP 6.0 on PATH.
#
#  SignalP lives in its own conda environment (python 3.10 + torch) and is not
#  otherwise reachable: ~/anaconda3/envs/signalp/bin is on nobody's PATH, so
#  two LowVan builds concluded it was not installed and skipped signal-peptide
#  prediction. This wrapper is the fix -- it execs the env's own interpreter
#  directly, so it needs no `conda activate` and works from any shell, cron
#  job or subprocess.
#
#  The model weights live in the package's model_weights/ directory and are NOT
#  installed by pip; they come from the licensed tarball, which is on this
#  machine at
#    .../CEPI/LowVan_Skill_Development/Rhabdoviridae/box_folder/Claude/signalp-6.0i.fast.tar.gz
#  Restored 18 Nov 2026 by copying models/distilled_model_signalp6.pt into
#    $SP_ENV/lib/python3.10/site-packages/signalp/model_weights/
#  Only the FAST model is in that tarball, so --mode fast works and --mode slow
#  does not.
#
#  Installed 18 Nov 2026.
SP_ENV="${SIGNALP_ENV:-$HOME/anaconda3/envs/signalp}"
if [ ! -x "$SP_ENV/bin/signalp6" ]; then
    echo "signalp6: environment not found at $SP_ENV" >&2
    echo "  set SIGNALP_ENV to the conda env that holds it." >&2
    exit 127
fi
W="$SP_ENV/lib/python3.10/site-packages/signalp/model_weights"
if [ ! -e "$W/distilled_model_signalp6.pt" ]; then
    echo "signalp6: the package is installed but its model weights are missing." >&2
    echo "  Expected: $W/distilled_model_signalp6.pt" >&2
    echo "  Restore from the licensed tarball already on this machine:" >&2
    echo "    tar -xzf <...>/Rhabdoviridae/box_folder/Claude/signalp-6.0i.fast.tar.gz -C /tmp/sp" >&2
    echo "    cp /tmp/sp/signalp6_fast/signalp-6-package/models/*.pt '$W/'" >&2
    exit 127
fi
exec "$SP_ENV/bin/signalp6" "$@"
