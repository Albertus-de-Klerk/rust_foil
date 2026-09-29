#!/usr/bin/env bash
# Generates golden reference data under tests/golden/ from the reference QFoil builds.
#
# Usage:  tools/gen_golden.sh            full matrix
#         tools/gen_golden.sh --quick    NACA 0012, Re 1e6, alpha 0 5 10 only (smoke test)
#
# Prerequisite: tools/build_reference.sh (builds qfoil, qfoil-dump, qfoil-auto).
#
# Every operating point is a FRESH QFoil process with a single ALFA command. This is how
# QBlade drives QFoil (see docs/PORTING_PLAN.md, S10). In-process sequential ALFA/ASEQ
# sweeps diverge to NaN after the first point in QFoil 0.9.
#
# For every point the script runs:
#   plain   reference-build/qfoil         -> polar row (original format), Cp/BL dumps, stdout
#   dump    reference-build/qfoil-dump    -> 'final' record (full-precision CL, CD, ...)
#   auto    reference-build/qfoil-auto    -> 'final' record; must equal dump bit-for-bit (S9 probe)
#   shipped reference/bin/QFoil           -> stdout; must equal plain (checks our build = QBlade's)
# and records whether plain/dump stdout are identical (dump instrumentation is side-effect free).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN_PLAIN="$ROOT/reference-build/qfoil/bin/qfoil"
BIN_DUMP="$ROOT/reference-build/qfoil-dump/bin/qfoil"
BIN_AUTO="$ROOT/reference-build/qfoil-auto/bin/qfoil"
BIN_SHIP_SRC="$ROOT/reference/bin/QFoil"
BIN_SHIP="$ROOT/reference-build/shipped/QFoil"
WORK="$ROOT/reference-build/golden-work"
OUT="$ROOT/tests/golden"
AIRFOILS="$OUT/airfoils"

ITMAX=100          # ITER sent to QFoil (QBlade sends a user-set ITER; QFoil default is 20)
NCRIT=9
TIMEOUT=300
KEEP_ALPHAS=" 0.0 5.0 10.0 15.0 "   # per-station Cp / BL kept for these

# case id | load command
CASES=(
  "naca0012|NACA 0012"
  "naca0020|NACA 0020"
  "naca4412|NACA 4412"
  "du91w2250|LOAD du91w2250.dat"
  "e387|LOAD e387.dat"
)
RES=(1e5 1e6 5e6)

# (case re alpha [vm]) for which ALL intermediate dumps are kept
FULL_DUMPS=(
  "naca0012 1e6 0.0"
  "naca0012 1e6 5.0"
  "naca0012 1e6 10.0"
  "naca0012 1e6 15.0"
  "naca4412 1e6 5.0"
  "naca4412 1e6 15.0"
  "e387 1e5 5.0"
  "du91w2250 1e6 10.0"
  "naca0012 1e6 5.0 vm"
)

load_cmd() {
  local id="$1" c
  for c in "${CASES[@]}"; do
    [[ "${c%%|*}" == "$id" ]] && { echo "${c#*|}"; return; }
  done
  echo "unknown case $id" >&2; exit 2
}

write_cmd() {  # $1=dir $2=load-cmd $3=re $4=alpha
  local load="$2"
  if [[ "$load" == LOAD* ]]; then
    # QFoil's command parser truncates long paths: copy the file in, load by short name
    cp "$AIRFOILS/${load#LOAD }" "$1/foil.dat"
    load="LOAD foil.dat"
  fi
  printf '%s\n' PLOP "G F" "" "$load" OPER VPAR "N $NCRIT" "" "ITER $ITMAX" \
    "VISC $3" PACC pol.txt "" "ALFA $4" "CPWR cp.txt" "DUMP bl.txt" "" QUIT > "$1/cmd.txt"
}

run_bin() {  # $1=bin $2=dir ; extra env from caller ; stdout -> out.log
  ( cd "$2" && ulimit -s unlimited && timeout "$TIMEOUT" "$1" < cmd.txt > out.log 2>&1 ) || echo $? > "$2/exitcode"
}

# ---------------------------------------------------------------- one operating point
point() {
  local id="$1" re="$2" a="$3"
  local load; load="$(load_cmd "$id")"
  local base="$WORK/$id/re$re/a$a"
  rm -rf "${base:?}"; mkdir -p "$base"/{plain,dump,auto,ship}
  local v
  for v in plain dump auto ship; do write_cmd "$base/$v" "$load" "$re" "$a"; done

  run_bin "$BIN_PLAIN" "$base/plain"
  mkdir -p "$base/dump/d" "$base/auto/d"
  QFOIL_DUMP="$base/dump/d" QFOIL_DUMP_TAGS=",final,viscal_end," run_bin "$BIN_DUMP" "$base/dump"
  QFOIL_DUMP="$base/auto/d" QFOIL_DUMP_TAGS=",final," run_bin "$BIN_AUTO" "$base/auto"
  run_bin "$BIN_SHIP" "$base/ship"

  local fdump fauto
  fdump="$(ls "$base"/dump/d/*_final.txt 2>/dev/null | head -1 || true)"
  fauto="$(ls "$base"/auto/d/*_final.txt 2>/dev/null | head -1 || true)"

  local same_dump same_auto same_ship
  cmp -s "$base/plain/out.log" "$base/dump/out.log" && same_dump=1 || same_dump=0
  if [[ -n "$fdump" && -n "$fauto" ]] && cmp -s <(grep -v '^#TAG' "$fdump") <(grep -v '^#TAG' "$fauto"); then
    same_auto=1; else same_auto=0; fi
  # the shipped binary's runtime prints an IEEE-flags note at exit; ignore that line only
  local ieee='^Note: The following floating-point exceptions'
  cmp -s <(grep -v "$ieee" "$base/plain/out.log") <(grep -v "$ieee" "$base/ship/out.log") \
    && same_ship=1 || same_ship=0

  local niter failed nan
  niter=$(grep -c '^ *a = ' "$base/plain/out.log" || true)
  failed=$(grep -c 'VISCAL:  Convergence failed' "$base/plain/out.log" || true)
  nan=$(grep -c 'NaN' "$base/plain/out.log" || true)

  # summary line: id re alpha niter viscal_failed nan_lines same_dump same_auto same_ship final_file
  echo -e "$id\t$re\t$a\t$niter\t$failed\t$nan\t$same_dump\t$same_auto\t$same_ship\t${fdump:-none}" > "$base/summary.tsv"

  # prune big stdout of non-kept points (keep the convergence history tail)
  tail -n 400 "$base/plain/out.log" > "$base/plain/out.tail.log"
  if [[ -z "${GOLDEN_KEEP_LOGS:-}" ]]; then
    rm -f "$base"/{plain,dump,auto,ship}/out.log "$base"/{auto,ship}/{pol.txt,cp.txt,bl.txt}
  fi
}

# ---------------------------------------------------------------- full intermediate dump
full_dump() {
  local id="$1" re="$2" a="$3" vm="${4:-}"
  local load; load="$(load_cmd "$id")"
  local name="${id}_re${re}_a${a}${vm:+_vm}"
  local dir="$WORK/fulldump/$name"
  rm -rf "${dir:?}"; mkdir -p "$dir/d"
  write_cmd "$dir" "$load" "$re" "$a"
  if [[ -n "$vm" ]]; then
    QFOIL_DUMP="$dir/d" QFOIL_DUMP_VM=1 QFOIL_DUMP_TAGS=",setbl," run_bin "$BIN_DUMP" "$dir"
  else
    QFOIL_DUMP="$dir/d" run_bin "$BIN_DUMP" "$dir"
  fi
  rm -rf "${OUT:?}/dumps/$name"; mkdir -p "$OUT/dumps/$name"
  cp "$dir"/d/*.txt "$OUT/dumps/$name/"
  cp "$dir/cmd.txt" "$OUT/dumps/$name/cmd.txt"
}

# ---------------------------------------------------------------- committed fixtures
# Trimmed copies of selected full dumps: every single-stage record plus the Newton
# iterations 1, 2 and last (setbl/blsolv/update/iter). The default `cargo test` uses these.
# Full dumps (tests/golden/dumps/, gitignored) serve the `--ignored` tests.
FIXTURES=(naca0012_re1e6_a5.0 naca4412_re1e6_a5.0 e387_re1e5_a5.0)
make_fixtures() {
  local p src dst f tag k n
  for p in "${FIXTURES[@]}"; do
    src="$OUT/dumps/$p"; dst="$OUT/fixtures/$p"
    [[ -d "$src" ]] || continue
    rm -rf "${dst:?}"; mkdir -p "$dst"
    cp "$src/cmd.txt" "$dst/"
    n=$(ls "$src"/*_setbl.txt | wc -l)
    for tag in setbl blsolv update iter; do
      k=0
      for f in $(ls "$src"/*_"$tag".txt); do
        k=$((k+1))
        [[ $k == 1 || $k == 2 || $k == "$n" ]] && cp "$f" "$dst/"
      done
    done
    for f in "$src"/*.txt; do
      case "$f" in *_setbl.txt|*_blsolv.txt|*_update.txt|*_iter.txt|*/cmd.txt) ;; *) cp "$f" "$dst/" ;; esac
    done
  done
}

# ---------------------------------------------------------------- assemble
val() {  # $1=file $2=record-name  -> scalar value
  awk -v n="$2" '($1=="#S"||$1=="#I") && $2==n {print $3; exit}' "$1"
}
vec() {  # $1=file $2=name $3=index(1-based)
  awk -v n="$2" -v k="$3" '$1=="#R"&&$2==n{f=1;i=0;next} f{i++; if(i==k){print $1; exit}}' "$1"
}

assemble() {
  mkdir -p "$OUT/polars" "$OUT/points"
  cat "$WORK"/*/re*/a*/summary.tsv | sort -k1,1 -k2,2 -k3,3g > "$OUT/checks.tsv.tmp"
  { echo -e "case\tre\talpha\tniter_printed\tviscal_failed\tnan_lines\tplain_eq_dump_stdout\tdump_eq_auto_final\tplain_eq_shipped_stdout"
    cut -f1-9 "$OUT/checks.tsv.tmp"; } > "$OUT/checks.tsv"

  local c id re
  for c in "${CASES[@]}"; do
    id="${c%%|*}"
    for re in "${RES[@]}"; do
      [[ -d "$WORK/$id/re$re" ]] || continue
      local tsv="$OUT/polars/${id}_re${re}.tsv" pol="$OUT/polars/${id}_re${re}.pol"
      echo -e "alpha\tconverged\tniter\tCL\tCD\tCDp\tCM\tCDf\tXtr_top\tXtr_bot\trmsbl" > "$tsv"
      local header_done=0 rows=""
      local a
      for a in $(ls "$WORK/$id/re$re" | sed 's/^a//' | sort -g); do
        local base="$WORK/$id/re$re/a$a" f
        f="$(ls "$base"/dump/d/*_final.txt 2>/dev/null | head -1 || true)"
        if [[ -n "$f" ]]; then
          local conv it
          conv=$(val "$f" CONVERGED); it=$(val "$f" ITER)
          [[ "$conv" == 1 ]] || it=$ITMAX
          echo -e "$a\t$conv\t$it\t$(val "$f" CL)\t$(val "$f" CD)\t$(val "$f" CDP)\t$(val "$f" CM)\t$(val "$f" CDF)\t$(vec "$f" XOCTR 1)\t$(vec "$f" XOCTR 2)\t$(val "$f" RMSBL)" >> "$tsv"
        else
          echo -e "$a\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> "$tsv"
        fi
        # original-format polar: header from first run, converged rows in alpha order
        if [[ -f "$base/plain/pol.txt" ]]; then
          if [[ $header_done == 0 ]]; then
            sed -n '1,/^ *------/p' "$base/plain/pol.txt" > "$pol"; header_done=1
          fi
          sed -n '/^ *------/,$p' "$base/plain/pol.txt" | tail -n +2 >> "$pol"
        fi
        if [[ "$KEEP_ALPHAS" == *" $a "* ]]; then
          local pd="$OUT/points/${id}_re${re}/a$a"
          mkdir -p "$pd"
          for x in cp.txt bl.txt; do [[ -f "$base/plain/$x" ]] && cp "$base/plain/$x" "$pd/$x"; done
          local ve; ve="$(ls "$base"/dump/d/*_viscal_end.txt 2>/dev/null | head -1 || true)"
          [[ -n "$ve" ]] && cp "$ve" "$pd/viscal_end.txt"
          cp "$base/plain/out.tail.log" "$pd/stdout_tail.log"
        fi
      done
    done
  done
  rm -f "$OUT/checks.tsv.tmp"
}

# ---------------------------------------------------------------- main
if [[ "${1:-}" == "--point" ]]; then shift; point "$@"; exit 0; fi
if [[ "${1:-}" == "--fulldump" ]]; then shift; full_dump "$@"; exit 0; fi
if [[ "${1:-}" == "--fixtures" ]]; then make_fixtures; exit 0; fi

for b in "$BIN_PLAIN" "$BIN_DUMP" "$BIN_AUTO"; do
  [[ -x "$b" ]] || { echo "missing $b: run tools/build_reference.sh" >&2; exit 1; }
done
mkdir -p "$(dirname "$BIN_SHIP")"; install -m 0755 "$BIN_SHIP_SRC" "$BIN_SHIP"

ALPHAS="$(awk 'BEGIN{for(i=-20;i<=40;i++) printf "%.1f\n", i/2}')"
QUICK=0
if [[ "${1:-}" == "--quick" ]]; then
  QUICK=1; CASES=("naca0012|NACA 0012"); RES=(1e6); ALPHAS=$'0.0\n5.0\n10.0'; FULL_DUMPS=()
fi

rm -rf "${WORK:?}"; mkdir -p "$WORK"
jobs=()
for c in "${CASES[@]}"; do
  for re in "${RES[@]}"; do
    for a in $ALPHAS; do jobs+=("${c%%|*} $re $a"); done
  done
done
echo "running ${#jobs[@]} points x 4 binaries on $(nproc) cores ..."
printf '%s\n' "${jobs[@]}" | xargs -P "$(nproc)" -L 1 bash "$0" --point
for fd in "${FULL_DUMPS[@]}"; do bash "$0" --fulldump $fd; done

assemble
make_fixtures
{
  echo "generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "quick: $QUICK"
  echo "itmax: $ITMAX  ncrit: $NCRIT"
  for v in qfoil qfoil-dump qfoil-auto; do echo "--- $v"; cat "$ROOT/reference-build/$v/BUILD_INFO.txt"; done
  echo "--- shipped: $(sha256sum "$BIN_SHIP_SRC" | cut -d' ' -f1)"
  echo "--- patch: $(sha256sum "$ROOT/tools/patches/dump/0001-golden-dumps.patch" | cut -d' ' -f1)"
} > "$OUT/GENERATION_INFO.txt"
echo "done -> $OUT"
