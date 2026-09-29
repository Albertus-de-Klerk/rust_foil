#!/usr/bin/env bash
# Summarises tests/golden/: convergence per case and the three bit-identity checks.
set -euo pipefail
G="$(cd "$(dirname "$0")/.." && pwd)/tests/golden"

echo "## Convergence per case (fresh process per alpha, ITER 100)"
echo
echo "| case | Re | points | converged | unconverged alphas |"
echo "|---|---|---|---|---|"
for f in "$G"/polars/*.tsv; do
  b="$(basename "$f" .tsv)"; c="${b%_re*}"; re="${b##*_re}"
  awk -F'\t' -v c="$c" -v re="$re" 'NR>1{n++; if($2==1) k++; else u=u (u?", ":"") $1}
    END{printf "| %s | %s | %d | %d | %s |\n", c, re, n, k, (u?u:"none")}' "$f"
done
echo
echo "## Bit-identity checks (checks.tsv)"
echo
awk -F'\t' 'NR>1{n++; d+=$7; a+=$8; s+=$9; if($6>0) nan++}
  END{printf "* points: %d\n* plain == dump-build stdout: %d/%d\n* dump == -fautomatic final record: %d/%d\n* plain == shipped QFoil stdout: %d/%d\n* points with NaN in stdout: %d\n", n, d, n, a, n, s, n, nan}' "$G/checks.tsv"
echo
echo "Mismatches (case re alpha dump auto shipped):"
awk -F'\t' 'NR>1 && ($7!=1 || $8!=1 || $9!=1){print "  " $1, $2, $3, $7, $8, $9}' "$G/checks.tsv"
