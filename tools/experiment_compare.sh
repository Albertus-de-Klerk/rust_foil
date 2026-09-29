#!/usr/bin/env bash
# Runs the golden point matrix with an experimental reference variant and compares it to
# tests/golden/polars/*.tsv (the unmodified reference).
#
# Usage: tools/experiment_compare.sh <variant> [case-filter-regex]
#   <variant>  a dump-patched build in reference-build/<variant>/ (see build_reference.sh)
#
# Reports per case: converged points (reference -> variant), points gained/lost, and the
# max |dCL|, |dCM| and max relative dCD over points converged in both.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VAR="$1"; FILTER="${2:-.}"
BIN="$ROOT/reference-build/$VAR/bin/qfoil"
[[ -x "$BIN" ]] || { echo "build it first: tools/build_reference.sh $VAR" >&2; exit 1; }
WORK="$ROOT/reference-build/experiment-$VAR"
G="$ROOT/tests/golden"

run_point() {  # case re alpha
  local id="$1" re="$2" a="$3" dir="$WORK/$1/re$2/a$3"
  rm -rf "${dir:?}"; mkdir -p "$dir/d"
  local load
  case "$id" in
    naca*) load="NACA ${id#naca}" ;;
    *)     cp "$G/airfoils/$id.dat" "$dir/foil.dat"; load="LOAD foil.dat" ;;
  esac
  printf '%s\n' PLOP "G F" "" "$load" OPER VPAR "N 9" "" "ITER 100" \
    "VISC $re" "ALFA $a" "" QUIT > "$dir/cmd.txt"
  ( cd "$dir" && ulimit -s unlimited && QFOIL_DUMP="$dir/d" QFOIL_DUMP_TAGS=",final," \
      timeout 300 "$BIN" < cmd.txt > /dev/null 2>&1 ) || true
}
export -f run_point; export WORK G BIN

if [[ "${EXP_SKIP_RUN:-}" != 1 ]]; then
  rm -rf "${WORK:?}"
  for f in "$G"/polars/*.tsv; do
    b="$(basename "$f" .tsv)"; [[ "$b" =~ $FILTER ]] || continue
    awk -F'\t' -v c="${b%_re*}" -v re="${b##*_re}" 'NR>1{print c, re, $1}' "$f"
  done | xargs -P "$(nproc)" -L 1 bash -c 'run_point "$@"' _
fi

printf '| case | Re | conv ref → var | gained | lost | max|ΔCL| | max|ΔCM| | max rel ΔCD |\n|---|---|---|---|---|---|---|---|\n'
for f in "$G"/polars/*.tsv; do
  b="$(basename "$f" .tsv)"; [[ "$b" =~ $FILTER ]] || continue
  c="${b%_re*}"; re="${b##*_re}"
  while IFS=$'\t' read -r a conv _ cl cd _ cm _; do
    [[ "$a" == alpha ]] && continue
    ff="$(ls "$WORK/$c/re$re/a$a"/d/*_final.txt 2>/dev/null | head -1 || true)"
    if [[ -n "$ff" ]]; then
      awk -v a="$a" -v rc="$conv" -v rcl="$cl" -v rcd="$cd" -v rcm="$cm" '
        $1=="#I"&&$2=="CONVERGED"{v=$3} $1=="#S"&&$2=="CL"{cl=$3} $1=="#S"&&$2=="CD"{cd=$3} $1=="#S"&&$2=="CM"{cm=$3}
        END{print a, rc, v, rcl, cl, rcd, cd, rcm, cm}' "$ff"
    else
      echo "$a $conv 0 NA NA NA NA NA NA"
    fi
  done < "$f" | awk -v c="$c" -v re="$re" '
    function abs(x){return x<0?-x:x}
    {n++; r+=$2; v+=$3; if($2==0&&$3==1) g++; if($2==1&&$3==0) l++;
     if($2==1&&$3==1){d=abs($5-$4); if(d>dcl)dcl=d; d=abs($9-$8); if(d>dcm)dcm=d; d=abs($7-$6)/abs($6); if(d>dcd)dcd=d}}
    END{printf "| %s | %s | %d → %d | %d | %d | %.2e | %.2e | %.2e |\n", c, re, r, v, g, l, dcl, dcm, dcd}'
done
