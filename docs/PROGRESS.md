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

### Phase 2 (in progress)

Workspace: `crates/qfoil-core` (library, no I/O), `crates/qfoil-cli` (binary `qfoil`, stub),
`crates/qfoil-golden` (test support: dump/polar readers in `BTreeMap`s, ULP comparison and the
`assert_golden!` macro). Lints: `unsafe_code = forbid`, `missing_docs = warn`, clippy `-D warnings`.

| step | module(s) | routines | parity |
|---|---|---|---|
| 1 | `spline`, `geometry` (`dat`, `naca`) | SPLINE, SPLIND, SEGSPL(D), TRISOL, SEVAL, DEVAL, D2VAL, CURV(S), SINVRT, SCALC; AREAD/GETFLT, LOAD (orientation), NACA/NACA4/NACA5, LEFIND | **bit-identical** (0 ULP): buffer XB, YB, SB, XBP, YBP for NACA 0012, 4412, E387, DU 91-W2-250 |
| 2 | `paneling` | PANGEN, ABCOPY, TECALC (geometry), NCALC, APCALC | **bit-identical** (0 ULP): X, Y, S, XP, YP, NX, NY, APANEL, SLE, LE, TE, CHORD, ANTE, ASTE, DSTE |
| CLI | `qfoil-cli` (binary `qfoil`, clap) | POLWRIT (polar file writer, Fortran `Fw.d` formatting), TINDEX (Itr columns, in SETBL) | **all 15 golden polar files byte-identical** to QFoil's `PACC` output (`cargo test --release -p qfoil-cli -- --ignored`) |
| 7–8 | `viscous` (VISCAL loop, CDCALC◆ D10), `polar` (public API `analyse_polar`, `PreparedAirfoil`) | VISCAL (QVFUE, GAMQV, STMOVE, CLCALC, CDCALC with GWAKE), ALFA driver | **all 915 golden points bit-identical, 0 convergence-flag mismatches** (`cargo test --release --test golden_polars -- --ignored`), incl. non-converged and deep-stall points. Every dumped Newton iteration of 8 full-dump runs bit-identical |
| 6 | `newton` (NewtonSystem), `bl::coupling` | SETBL, BLSOLV, UPDATE◆ (D6–D8), UESET; STFIND, IBLPAN, XICALC (incl. WGAP), IBLSYS, UICALC, QVFUE, STMOVE; QFoil VISCAL initialisation (D9) | **bit-identical** (0 ULP): BL set-up from scratch vs `viscal_init` (IST, SST, pointers, ξ, WGAP, UINV, clamped UEDG); first Newton iteration: SETBL VA/VB/VDEL, full VM (~10⁵ entries), derived arrays; BLSOLV solution; UPDATE BL arrays, CL, RLX, RMSBL, RMXBL, on 7 runs |
| 5 | `bl::station`, `bl::equations`, `bl::march`, `bl` (BoundaryLayer, Side) | BLPRV, BLKIN, BLVAR, BLMID, BLDIF◆ (D1–D3: tanh Kc(Hk), SCC_HKA Jacobian, wake ALD), TRCHEK2, TRDIF, BLSYS, TESYS, MRCHUE, MRCHDU◆ (D5), XIFSET, DSLIM | **bit-identical** (0 ULP): all BL arrays after MRCHUE (θ, δ*, Cτ/N, Ue, mass, τ, D, Cτeq, δ, θ*, ITRAN, XSSITR), primary arrays after MRCHDU, on 7 runs incl. transition, separation (inverse mode) and wake. Station Jacobians VS1/VS2 match finite differences for laminar, wake and turbulent intervals, except the upstream S15 approximation, which is isolated and verified |
| 4 | `bl::closure`, `bl::transition`, `settings::BlParams` | HKIN, HCT, HSL, HST◆ (QFoil HSMIN/DHSINF), CFL, CFT, DIL, DILW, DIT; DAMPL, DAMPL2, AXSET; BLPINI | **bit-identical** (0 ULP) on all 6,331 rows of the Fortran reference table (`tools/fdrivers/closures.f`, every branch); every analytic partial matches central differences (rel 1e-5), except DILW ∂/∂HK, an upstream sign error (S14) |
| 3 | `linalg`, `inviscid` (`influence`, `forces`), `wake`, `settings`, `operating` | LUDCMP, BAKSUB, GAUSS; PSILIN, PSWLIN; GGCALC (incl. sharp-TE branch, E387), SPECAL, QISET, TECALC (strengths), MRCL, COMSET, CPCALC, CLCALC; XYWAKE, SETEXP, ATANC, QWCALC, QDCALC | **bit-identical** (0 ULP): AIJ, BIJ, LU factors and pivots, GAMU, GAM, QINV, CPI, CL, CM, CDP, wake x/y/s/normals/angles, QINVU, DIJ, for NACA 0012, 4412, E387 (fixtures) and DU 91-W2-250, NACA 0012 α 0/15, NACA 4412 α 15 (`--ignored`) |

Not ported (no effect on results): GEOPAR, NORM (`LNORM` is off by default), SPLINA, SPLNXY, CANG;
PSILIN's `GEOLIN` branch (inverse design) and ground-effect images (`LIMAGE`, not reachable from QFoil's CLI).

## Next

* Phase 2 wrap-up: idiomatic refactor commit (behaviour-preserving, guarded by the bit-parity tests).
* Phase 3: `docs/VALIDATION.md` with the per-case table and `plotters` plots. Parity is already
  exact, so every acceptance criterion holds with zero error.
* Phase 4: rayon over α (points are independent), criterion benchmarks against the reference binary.

## Decisions (2026-09-29)

* **Q1.** Independent α points (QBlade semantics). There's no warm start, and the sweep parallelises per point.
* **Q2.** MRCHDU relaxation cap tested and **not beneficial** (740 → 741 converged, scattered +22/−21). Ported as in the reference.
* **Q3.** E387 golden case now uses a 160-node re-panelled file: convergence 88 → 128 of 183 points.
* **Q4.** Full dumps gitignored. A trimmed fixture set (~9 MB) is committed for the default test run, and `--ignored` tests use the full dumps.
* **Q5.** `OperatingPoint` struct approved.
* **Q6.** Git repository initialised.

## Open questions (new)

* **S14 / S15.** Two upstream XFOIL Jacobian defects were found by the finite-difference checks: a sign error in DILW ∂/∂Hk, and a missing Rθ path in the turbulent shear-lag row. Both are ported faithfully. After parity: fix them (opt-in or default), or leave them as is? Fixing changes the convergence behaviour but not converged answers.

## Known deviations of the Rust port from the reference

None in results: all 915 golden points are bit-identical. Structural differences that do not
affect results: no fixed array sizes except where QFoil's behaviour depends on them (NACA point
count, input limits); LUDCMP's 499-node limit and ABCOPY/PANGEN limits are errors instead of
`STOP`; `GEOPAR`, `NORM`, plotting, inverse design and `TINDEX` are not ported.

### Parity pitfalls found (for future work)
1. Fortran `a*b**n` must be written `a * powi(b, n)`: gfortran evaluates the power first
   (a wake inverse-mode formula in MRCHUE differed by 12–38 ULP before the fix).
2. LLVM rewrites `pow(x, 0.5)` to `sqrt` without fast-math; glibc `pow` can differ by 1 ULP.
   All real-exponent powers go through `fortran::pow`, which hides the exponent with
   `std::hint::black_box`.
3. `a += b + c` is `a + (b + c)`; Fortran `A = A + B + C` is `(A + B) + C`.
