#!/usr/bin/env bash
# create_envs.sh: build the small, separate conda environments (one per pipeline stage).
# Separate environments avoid the dependency conflicts that a single large environment causes.
#
# Usage: bash create_envs.sh [env-name ...]     (no arguments = build all)
# A failed environment does not stop the others; a summary is printed at the end.

set -uo pipefail

conda config --set channel_priority strict

ALL=(chol-qc chol-asm chol-checkm2 chol-amr chol-phylo chol-r)
TARGETS=("${@:-${ALL[@]}}")

OK=(); FAIL=()
for e in "${TARGETS[@]}"; do
  echo "=============== building $e ==============="
  if conda env list | awk '{print $1}' | grep -qx "$e"; then
    echo "$e already exists, skipping (remove with: conda env remove -n $e)"
    OK+=("$e (existing)"); continue
  fi
  if conda env create -f "envs/${e}.yml"; then OK+=("$e"); else FAIL+=("$e"); fi
done

echo
echo "Built:  ${OK[*]:-none}"
echo "Failed: ${FAIL[*]:-none}"
[[ ${#FAIL[@]} -gt 0 ]] && echo "Send me the error text for the failed environment(s)."
exit 0
