#!/usr/bin/env bash
# Builds the reference QFoil binary out-of-tree in reference-build/.
#
#   reference-build/qfoil/bin/qfoil        unmodified QFoil 0.9 source, original compile flags
#   reference-build/qfoil-dump/bin/qfoil   same + tools/patches/dump/*.patch (intermediate dumps)
#   reference-build/qfoil-auto/bin/qfoil   dump patch, -fautomatic (SAVE-dependence probe, S9)
#
# reference/ is never modified. Deviations from the original build, all build-only:
#   * plotlib/libPlt.a in the source zip is a MinGW (Windows) object; it is rebuilt
#     from plotlib/dummy_plot.f (empty stubs). That file is not valid fixed-form as
#     shipped (line 1 unindented, line 94 > 72 cols), so it is compiled with
#     -ffixed-line-length-none after indenting line 1.
#   * LDFLAGS: '-static -s' dropped. Arch's gcc-fortran ships no static libgfortran.
#     Numerics are unaffected (same objects, same glibc libm).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/reference/qfoil"
OUT="$ROOT/reference-build"
FFLAGS_ORIG="-O2 -std=legacy -fno-automatic -finit-local-zero -fno-align-commons -fdefault-real-8 -mtune=generic -flto"

build_variant() {
  local name="$1" fflags="$2"; shift 2
  local patchdirs=("$@")
  local dir="$OUT/$name"
  rm -rf "$dir"
  cp -r "$SRC" "$dir"
  chmod -R u+w "$dir"

  local pd p
  for pd in "${patchdirs[@]}"; do
    for p in "$pd"/*.patch; do
      [[ -e "$p" ]] || continue
      patch -d "$dir" -p1 --quiet < "$p"
    done
  done

  # plot stubs (see header)
  ( cd "$dir/plotlib"
    rm -f libPlt.a dummy_plot.o
    sed -i '1s/^SUBROUTINE/      SUBROUTINE/' dummy_plot.f
    gfortran -c -O2 -std=legacy -fdefault-real-8 -ffixed-line-length-none dummy_plot.f
    ar rcs libPlt.a dummy_plot.o )

  ( cd "$dir/bin"
    make clean >/dev/null 2>&1 || true
    make qfoil FFLAGS="$fflags" FFLOPT="$fflags" LDFLAGS="-flto" > build.log 2>&1 \
      || { tail -30 build.log; exit 1; } )

  {
    echo "variant: $name"
    gfortran --version | head -1
    echo "FFLAGS: $fflags"
    echo "LDFLAGS: -flto"
    echo "patches: ${patchdirs[*]:-none}"
  } > "$dir/BUILD_INFO.txt"
  echo "built $dir/bin/qfoil"
}

mkdir -p "$OUT"
variants="${*:-qfoil qfoil-dump qfoil-auto}"
for v in $variants; do
  case "$v" in
    qfoil)      build_variant qfoil "$FFLAGS_ORIG" ;;
    qfoil-dump) build_variant qfoil-dump "$FFLAGS_ORIG" "$ROOT/tools/patches/dump" ;;
    # SAVE-dependence probe (PORTING_PLAN S9): locals on the stack, uninitialised
    # reals are signalling NaN. Carries the dump patch so results compare bitwise.
    qfoil-auto) build_variant qfoil-auto \
                  "-O2 -std=legacy -fautomatic -finit-real=snan -finit-integer=-99999999 -fno-align-commons -fdefault-real-8 -mtune=generic -flto" \
                  "$ROOT/tools/patches/dump" ;;
    # experiment: PORTING_PLAN S1, MRCHDU relaxation capped at 0.7 (not the reference)
    qfoil-s1cap) build_variant qfoil-s1cap "$FFLAGS_ORIG" \
                  "$ROOT/tools/patches/dump" "$ROOT/tools/patches/experiments/s1" ;;
    *) echo "unknown variant $v" >&2; exit 2 ;;
  esac
done
