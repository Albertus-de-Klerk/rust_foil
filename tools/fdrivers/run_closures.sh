#!/usr/bin/env bash
# Builds tools/fdrivers/closures.f against the reference xblsys.f/xbl.f (original compile
# flags) and writes tests/golden/closures.txt.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SRC="$ROOT/reference/qfoil/src"
WORK="$ROOT/reference-build/fdrivers"
FFLAGS="-O2 -std=legacy -fno-automatic -finit-local-zero -fno-align-commons -fdefault-real-8 -mtune=generic"
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
# xblsys.f holds the closures; BLPINI (xbl.f) sets the BLPAR constants. Only BLPINI is
# taken from xbl.f, so the driver needs no XFOIL.INC state.
awk '/^      SUBROUTINE BLPINI/,/^      END/' "$SRC/xbl.f" > blpini.f
gfortran $FFLAGS -I"$SRC" -c "$SRC/xblsys.f" blpini.f "$ROOT/tools/fdrivers/closures.f"
gfortran -o closures closures.o xblsys.o blpini.o 2> link.log || {
  # xblsys.f references GAUSS etc. only from routines the driver never calls
  gfortran $FFLAGS -c "$SRC/xsolve.f"
  gfortran -o closures closures.o xblsys.o blpini.o xsolve.o
}
./closures
cp closures.txt "$ROOT/tests/golden/closures.txt"
wc -l "$ROOT/tests/golden/closures.txt"
