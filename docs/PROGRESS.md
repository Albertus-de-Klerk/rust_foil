# Progress

## Done

### Phase 0: survey
* `docs/PORTING_PLAN.md`: routine table, COMMON → struct mapping, call graph, and the 14 QFoil deviations (D1–D14, file:line, each confirmed by diffing against the bundled XFOIL 6.99 and against upstream MIT 6.99). Suspected issues S1–S12.
* Reference sources unpacked read-only: `reference/qfoil/` (QFOIL 0.9, GPL v2), `reference/xfoil-6.99/`, and the shipped binaries in `reference/bin/`.
* `LICENSE` (GPL v2 text from QFoil) and `NOTICE` (provenance) added.

### Phase 1: golden harness
* `tools/build_reference.sh` builds three variants in `reference-build/`:
  `qfoil` (original flags), `qfoil-dump` (+ `tools/patches/dump/0001-golden-dumps.patch`), and `qfoil-auto` (dump patch + `-fautomatic`, S9 probe).
  Build-only deviations: plot stub rebuilt for Linux, dynamic link (S11).
* `tools/gen_golden.sh` runs one fresh process per point (QBlade's usage, S10), 915 points × 4 binaries. Summary: `tools/golden_summary.sh`.
* Output in `tests/golden/`. Layout and formats are in `tests/golden/README.md`.

#### Harness self-checks (all points)

| check | result |
|---|---|
| plain build == dump build (stdout, byte-identical) | 915/915 |
| dump build == `-fautomatic`/sNaN build (final record, bit-identical) | 915/915, so **no SAVE-variable dependence**: the port may use plain locals |
| our build == QBlade's shipped `QFoil` (stdout, byte-identical except the runtime IEEE note) | 915/915 |
| NaN anywhere in stdout | 0 points |

#### Reference convergence (ITER 100)

| case | 1e5 | 1e6 | 5e6 |
|---|---|---|---|
| naca0012 | 49/61 | 61/61 | 61/61 |
| naca0020 | 53/61 | 59/61 | 61/61 |
| naca4412 | 53/61 | 55/61 | 59/61 |
| du91w2250 | 23/61 | 61/61 | 57/61 |
| e387 (160 nodes) | 34/61 | 42/61 | 52/61 |

Unconverged α are listed per case by `tools/golden_summary.sh`.

## Next

Phase 2, step 1: create the Cargo workspace (`qfoil-core`, `qfoil-cli`) and port `spline` +
`geometry` (+ NACA4/5, .dat parsing), tested against `dumps/*/…_pangen.txt` (XB, YB, SB, X, Y, S, XP, YP).

## Decisions (2026-09-29)

* **Q1.** Independent α points (QBlade semantics). There's no warm start, and the sweep parallelises per point.
* **Q2.** MRCHDU relaxation cap tested and **not beneficial** (740 → 741 converged, scattered +22/−21). Ported as in the reference.
* **Q3.** E387 golden case now uses a 160-node re-panelled file: convergence 88 → 128 of 183 points.
* **Q4.** Full dumps gitignored. A trimmed fixture set (~9 MB) is committed for the default test run, and `--ignored` tests use the full dumps.
* **Q5.** `OperatingPoint` struct approved.
* **Q6.** Git repository initialised.

## Known deviations of the Rust port from the reference

None yet (no Rust code).
