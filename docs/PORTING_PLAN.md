# QFoil → Rust porting plan

Status: Phase 0 (survey) and Phase 1 (golden harness) complete. Open questions in docs/PROGRESS.md.

## 1. Source survey

| Item | Finding |
|---|---|
| Reference | `reference/qfoil/` = **QFOIL 0.9** (`Qfoil_0.9_minimal_src.zip`, shipped inside `QBladeCE_2.0.9.7_unix.zip` → `Binaries/`) |
| Baseline | `reference/xfoil-6.99/` = `Xfoil_6.9.9_minimal_src.zip` from the same package |
| Upstream check | Baseline was diffed against MIT `xfoil6.99.tgz` (Dec 2013, fetched to scratch only, not kept). Baseline = upstream **plus** (a) raised array limits in `XFOIL.INC` and (b) the `.bl` debug file writer in `VISCAL` commented out. No other differences in `src/` or `osrc/`. |
| Licence | QFoil and XFOIL: **GPL v2** (`LICENSE_Qfoil`, `LICENSE_Xfoil`, source headers "version 2 … or (at your option) any later version"). QBlade itself is under the Academic Public License. That licence does not apply to QFoil. |
| Copyright | XFOIL © 2000 Mark Drela, Harold Youngren. QFOIL modifications © 2026 David Marten. |
| Language | Fortran 77, fixed form, legacy extensions (tab characters, `!` inline comments, `DO … ENDDO`), plus one C file (`osrc/getosfile.c`, unused by the solver path). |
| Build | GNU Make, `reference/qfoil/bin/Makefile`, gfortran. Flags: `-O2 -std=legacy -fno-automatic -finit-local-zero -fno-align-commons -fdefault-real-8 -flto`, static link. Plot library is a prebuilt dummy (`plotlib/libPlt.a` from `dummy_plot.f`). |
| Precision | `-fdefault-real-8`: every `REAL` and every real literal is 64-bit, so **`f64` is the correct target**. (The Makefile comment says single precision. The comment is wrong and the flag is what counts; see Suspected issues S6.) |
| Prebuilt binaries | `reference/bin/QFoil`, `reference/bin/XFoil` (x86-64 ELF, stripped). Used as a cross-check for our own reference build. |

### Physical/numerical scope that is actually ported

Only the analysis path: load or generate the airfoil → panel → inviscid solve → viscous Newton solve → polar sweep. Not ported: GDES/QDES/MDES (geometry and inverse design), plotting (`*plot*.f`, `xplots.f`, `plutil.f`, `gui.f`), `pplot`/`pxplot` programs, `profil.f`/`blu.f`, `polfit.f`, `modify*.f`, `xtcam.f`, `osmap.f` (Orr–Sommerfeld map, never called by the solver).

## 2. Call graph: OPER / polar entry → Newton solve

```
PROGRAM XFOIL (xfoil.f:21)
 ├─ INIT (xfoil.f:402) ── BLPINI (xbl.f:1593), MRCL, COMSET
 ├─ GETDEF 'xfoil.def' (xfoil.f:884)      ← harness must run with no xfoil.def present
 ├─ LOAD (xfoil.f:1257) ── AREAD, NORM, SCALC, SEGSPL, GEOPAR
 │   or NACA (xfoil.f:1613) ── NACA4/NACA5, SCALC, SEGSPL, GEOPAR, PANGEN
 ├─ PANGEN (xfoil.f:1668) ── SCALC, SEGSPL, LEFIND, TRISOL, CURV/SEVAL/DEVAL, TECALC, NCALC, APCALC
 └─ OPER (xoper.f:21)
     ├─ VISC/RE/ITER/PACC/VPAR(GW, N, XTR …)
     └─ ASEQ (xoper.f:545)  ── for IPOINT = 1..NPOINT:
         ├─ SPECAL (xoper.f:2754)
         │   ├─ GGCALC (xpanel.f:988) ── PSILIN, LUDCMP, BAKSUB     [only if !LGAMU || !LQAIJ]
         │   ├─ TECALC, QISET, MRCL, COMSET, CLCALC, CPCALC
         ├─ VISCAL(ITMAX+5) (xoper.f:2908)
         │   ├─ XYWAKE (xpanel.f:1270) ── PSILIN, SETEXP           [if !LWAKE; resets LWDIJ]
         │   ├─ QWCALC (xpanel.f:1127) ── PSILIN
         │   ├─ QISET
         │   ├─ [if !LIPAN] STFIND, IBLPAN, XICALC, IBLSYS
         │   ├─ UICALC
         │   ├─ ◆ QFOIL: UEDG ← min-clamp(UINV, 2.0); CTAU ← 0.01; LBLINI ← F
         │   ├─ [if LVCONV] QVFUE, CPCALC, GAMQV, CLCALC, CDCALC
         │   ├─ QDCALC (xpanel.f:1149) ── BAKSUB, PSILIN, PSWLIN   [if !LWDIJ || !LADIJ]
         │   └─ Newton loop ITER = 1..NITER:
         │       ├─ SETBL (xbl.f:21)
         │       │   ├─ MRCL, COMSET
         │       │   ├─ [!LBLINI] MRCHUE (xbl.f:542)   ◆ runs every alpha in QFoil
         │       │   │    └─ XIFSET, BLPRV, BLKIN, TRCHEK, BLSYS/TESYS, GAUSS, HKIN, DSLIM, BLVAR, BLMID
         │       │   ├─ MRCHDU (xbl.f:875)             ◆ RLX 0.7
         │       │   │    └─ (as MRCHUE)
         │       │   ├─ UESET (xpanel.f:1758)
         │       │   └─ per station: BLPRV, BLKIN, TRCHEK, BLSYS/TESYS, BLVAR, BLMID → VA/VB/VM/VDEL/VZ
         │       ├─ BLSOLV (xsolve.f:283)
         │       ├─ UPDATE (xbl.f:1256) ── DSLIM             ◆ RLX 0.7, DCL ±0.15, DHI/DLO
         │       ├─ MRCL, COMSET   (LALFA)   |  QISET, UICALC (!LALFA)
         │       ├─ QVFUE, GAMQV, STMOVE (xpanel.f:1628) ── STFIND, XICALC, IBLPAN, UICALC, IBLSYS
         │       ├─ CLCALC (xfoil.f:1105), CDCALC (xfoil.f:1189)  ◆ GWAKE drag
         │       └─ converged if RMSBL < 1e-4
         ├─ FCPMIN
         ├─ PLRADD/PLXADD → POLWRIT (iopol.f:432)   [if PACC and converged]
         └─ halt sweep after NSEQEX (=4) consecutive unconverged points

BL kernel (xblsys.f):
 BLSYS ── BLVAR ── HKIN, HSL, HST◆, CFL, CFT, DIL, DILW, HCT
       ── BLMID ── CFL, CFT
       ── BLDIF◆ ── AXSET ── DAMPL / DAMPL2
       ── TRDIF ── BLKIN, BLVAR, BLMID, BLDIF
 TRCHEK ── TRCHEK2 ── AXSET, BLKIN
 TESYS ── BLVAR
```
◆ = contains a QFoil modification.

## 3. Routine table (ported routines)

Line numbers refer to `reference/qfoil/src/`. The Rust module paths are all inside `qfoil-core` unless marked CLI.

### Step 1: splines and geometry

| File:line | Routine | Purpose | Calls | Rust module |
|---|---|---|---|---|
| spline.f:21 | SPLINE | cubic spline, zero-3rd-deriv ends | TRISOL | `spline` |
| spline.f:63 | SPLIND | spline with specified end derivatives | TRISOL | `spline` |
| spline.f:141 | SPLINA | Akima-ish spline (not in analysis path, port for completeness only if needed) | — | `spline` (optional) |
| spline.f:185 | TRISOL | tridiagonal solve | — | `spline` |
| spline.f:218/245/273 | SEVAL / DEVAL / D2VAL | spline value, 1st and 2nd derivative | — | `spline` |
| spline.f:302/349 | CURV / CURVS | curvature | — | `spline` |
| spline.f:404 | SINVRT | invert spline for s given x | SEVAL, DEVAL | `spline` |
| spline.f:436 | SCALC | arc length | — | `spline` |
| spline.f:533/562 | SEGSPL / SEGSPLD | segmented spline (corners) | SPLIND | `spline` |
| spline.f:452 | SPLNXY | x,y spline convenience | SCALC, SEGSPL | `spline` |
| xgeom.f:21 | LEFIND | LE location (s where tangent ⟂ chord) | SEVAL, DEVAL | `geometry` |
| xgeom.f:296 | NORM | normalise to unit chord | SCALC, SEGSPL, LEFIND | `geometry` |
| xgeom.f:323 | GEOPAR | area, thickness, camber etc. | LEFIND, AECALC, TCCALC | `geometry` |
| xgeom.f:387/520/226 | AECALC / TCCALC / SOPPS | geometric integrals, opposite-point | — | `geometry` |
| xgeom.f:1111 | CANG | panel corner angles (diagnostic) | — | `geometry` |
| naca.f:21/85 | NACA4 / NACA5 | NACA generator | — | `geometry::naca` |
| aread.f:2 | AREAD | .dat reader (plain/labeled/ISES/MSES) | — | parse `&str` in `geometry::dat` (pure); file I/O in CLI |
| xfoil.f:1257 | LOAD | load pipeline | AREAD, NORM, SCALC, SEGSPL, GEOPAR | `geometry::Geometry::from_coords` |
| xgdes.f:961 | **ABCOPY** | buffer → panel nodes **without re-panelling** (strips doubled points pairwise). This is the path LOAD takes, and so the path QBlade uses. | SCALC, SEGSPL, NCALC, LEFIND, SEVAL, TECALC, APCALC | `paneling::Paneling::from_buffer` |
| xutils.f:4/68 | SETEXP / ATANC | geometric-stretching distribution / continuous atan | — | `util` |

### Step 2: panelling

| File:line | Routine | Purpose | Calls | Rust module |
|---|---|---|---|---|
| xfoil.f:1668 | PANGEN | curvature-based node distribution (Newton on spacing) | SCALC, SEGSPL, LEFIND, TRISOL, CURV, TECALC, NCALC, APCALC | `paneling` |
| xfoil.f:2332 | TECALC | TE gap / sharpness (ANTE, ASTE, DSTE, SHARP) | — | `paneling` |
| xpanel.f:51 | NCALC | node normals | SEGSPL | `paneling` |
| xpanel.f:22 | APCALC | panel angles | — | `paneling` |

### Step 3: inviscid influence matrices and linear solve

| File:line | Routine | Purpose | Calls | Rust module |
|---|---|---|---|---|
| xpanel.f:99 | PSILIN | ψ and dψ/dn at a point from vortex+source panels (+ derivatives) | — | `inviscid::influence` |
| xpanel.f:803 | PSWLIN | ψ from wake source panels | — | `inviscid::influence` |
| xpanel.f:988 | GGCALC | build and LU-factor AIJ, solve α=0°, 90° γ | PSILIN, LUDCMP, BAKSUB | `inviscid` |
| xpanel.f:1127 | QWCALC | wake-point velocities for α=0°, 90° | PSILIN | `inviscid` |
| xpanel.f:1149 | QDCALC | source influence DIJ (airfoil + wake) | BAKSUB, PSILIN, PSWLIN | `inviscid` |
| xpanel.f:1597 | QISET | QINV = cosα·QINVU₁ + sinα·QINVU₂ | — | `inviscid` |
| xoper.f:2754 | SPECAL | converge inviscid at given α (Mach iteration) | GGCALC, TECALC, QISET, MRCL, COMSET, CLCALC, CPCALC | `inviscid` |
| xfoil.f:1105 | CLCALC | CL, CM, CDP with Karman–Tsien | — | `inviscid::forces` |
| xfoil.f:1074 | CPCALC | Cp | — | `inviscid::forces` |
| xfoil.f:1046 | COMSET | compressibility parameters | — | `inviscid::forces` |
| xfoil.f:801 | MRCL | Mach, Re from CL (RETYP/MATYP) | — | `settings` / `inviscid` |
| xsolve.f:174/250 | LUDCMP / BAKSUB | Crout LU with partial pivoting | — | `linalg` (hand-ported for FP-order parity; `faer`/`nalgebra` only after parity, behind a test) |
| xsolve.f:22 | GAUSS | small dense solve (4×4 station systems) | — | `linalg` |
| xoper.f:1768 | FCPMIN | min Cp | — | `inviscid::forces` |

### Step 4: BL closure relations and derivatives

| File:line | Routine | Purpose | Rust module |
|---|---|---|---|
| xblsys.f:2325 | HKIN | kinematic shape factor Hk(H, M²) | `bl::closure` |
| xblsys.f:2376 | HSL | laminar H* | `bl::closure` |
| xblsys.f:2437 | **HST** ◆ | turbulent H* (HSMIN, DHSINF changed) | `bl::closure` |
| xblsys.f:2403 | CFL | laminar Cf | `bl::closure` |
| xblsys.f:2535 | CFT | turbulent Cf | `bl::closure` |
| xblsys.f:2339 | DIL | laminar dissipation | `bl::closure` |
| xblsys.f:2357 | DILW | laminar wake dissipation | `bl::closure` |
| xblsys.f:2424 | DIT | turbulent dissipation (unused by BLVAR, but port) | `bl::closure` |
| xblsys.f:2566 | HCT | density shape factor H** | `bl::closure` |
| xblsys.f:2030 | DAMPL | e^N amplification rate (envelope) | `bl::transition` |
| xblsys.f:2148 | DAMPL2 | alternative rate (IDAMP=1) | `bl::transition` |
| xblsys.f:35 | AXSET | averaged amplification rate over interval | `bl::transition` |

Each gets an FD-vs-analytic derivative test.

### Step 5: BL march and transition

| File:line | Routine | Purpose | Rust module |
|---|---|---|---|
| xblsys.f:701 | BLPRV | set primary "2" station variables | `bl::station` |
| xblsys.f:725 | BLKIN | kinematic secondary vars + derivatives | `bl::station` |
| xblsys.f:784 | BLVAR | all secondary vars (closures) by ITYP | `bl::station` |
| xblsys.f:1124 | BLMID | midpoint Cf | `bl::station` |
| xblsys.f:1552 | **BLDIF** ◆ | discretised BL equations and Jacobian (lam/turb/wake) | `bl::equations` |
| xblsys.f:1195 | TRDIF | transition-interval equations | `bl::equations` |
| xblsys.f:583 | BLSYS | dispatch BLVAR/BLMID/BLDIF/TRDIF | `bl::equations` |
| xblsys.f:664 | TESYS | TE → wake junction system | `bl::equations` |
| xblsys.f:22/231 | TRCHEK / TRCHEK2 | transition check, implicit xtr | `bl::transition` |
| xbl.f:542 | MRCHUE | direct/inverse march with prescribed Ue | `bl::march` |
| xbl.f:875 | **MRCHDU** ◆ | march with Ue–Hk characteristic (quasi-inverse) | `bl::march` |
| xbl.f:1199 | XIFSET | forced-transition ξ | `bl::march` |
| xbl.f:1579 | DSLIM | δ* limiter (Hk ≥ HKLIM) | `bl::march` |

### Step 6: Newton system assembly and block solve

| File:line | Routine | Purpose | Rust module |
|---|---|---|---|
| xbl.f:21 | SETBL | assemble global Newton system (VA, VB, VM, VDEL, VZ) | `newton::assemble` |
| xbl.f:519 | IBLSYS | BL station → system row map | `newton` |
| xsolve.f:283 | BLSOLV | custom block-bidiagonal + dense mass-column solver (VACCEL drop tolerance) | `newton::solve` |
| xbl.f:1256 | **UPDATE** ◆ | relaxed update, RMS/max change, CL/α update | `newton::update` |

### Step 7: wake handling and viscous–inviscid coupling

| File:line | Routine | Purpose | Rust module |
|---|---|---|---|
| xpanel.f:1270 | XYWAKE | wake trajectory (streamline march), WGAP | `wake` (**needed by step 3 tests; ported with step 3, see §6**) |
| xpanel.f:1357 | STFIND | stagnation point panel IST, SST | `coupling` |
| xpanel.f:1395 | IBLPAN | BL ↔ panel pointers, NBL, IBLTE | `coupling` |
| xpanel.f:1455 | XICALC | BL arc length ξ (incl. wake), WGAP | `coupling` |
| xpanel.f:1542 | UICALC | UINV, UINV_A from QINV | `coupling` |
| xpanel.f:1562 | UECALC | UEDG from QVIS | `coupling` |
| xpanel.f:1580 | QVFUE | QVIS from UEDG | `coupling` |
| xpanel.f:1616 | GAMQV | GAM from QVIS | `coupling` |
| xpanel.f:1628 | STMOVE | shift BL arrays when the stagnation point moves | `coupling` |
| xpanel.f:1758 | UESET | Ue from mass defect via DIJ | `coupling` |
| xpanel.f:1786 | DSSET | δ* from mass | `coupling` |
| xfoil.f:1189 | **CDCALC** ◆ | Squire–Young + GWAKE correction, friction drag | `wake::drag` |

### Step 8: polar sweep driver

| File:line | Routine | Purpose | Rust module |
|---|---|---|---|
| xoper.f:2908 | **VISCAL** ◆ | viscous point driver | `solver::viscal` |
| xoper.f:545 | OPER/ASEQ | α sweep loop, NSEQEX halt, ITMAX+5 | `polar::sweep` |
| xoper.f:416 | OPER/ALFA | single point (ITMAX) | `polar::single` |
| xfoil.f:402, xbl.f:1593 | INIT, BLPINI | defaults | `settings::Settings::default()` |
| xoper.f:2490 | VPAR | viscous parameter edits (incl. ◆GW) | `settings` / CLI flags |
| xpol.f:550 | PLRADD | polar row (CDTOT, CDP, CDV, xtr …) | `polar` |
| iopol.f:432 | POLWRIT | polar file writer (format 9100) | **CLI** `qfoil-cli::polar_file` |
| xoper.f:1892/2265 | BLDUMP / CPDUMP | BL and Cp dumps | CLI (for golden comparison) |

## 4. Global and COMMON state → Rust structs

The Fortran has **no local state that is meant to persist**, but `-fno-automatic` makes every local variable static and `-finit-local-zero` zeroes it at program start. Phase 1 will check whether any routine reads a stale local (§7, S9).

| COMMON (file) | Contents | Rust owner |
|---|---|---|
| `/CR14/` (XFOIL.INC) | buffer airfoil XB, YB, XBP, YBP, SB, geometry params (THICKB, CAMBRB, …), camber/thickness arrays | **`Geometry`** |
| `/CR05/` | panel nodes X, Y, XP, YP, S (airfoil **and** wake), SLE, XLE, YLE, XTE, YTE, CHORD, WGAP, WAKLEN | **`Paneling`** (airfoil nodes, splines, LE/TE); wake nodes in `Paneling::wake` |
| `/CR06/` (part) | NX, NY, APANEL, SHARP, DSTE, ANTE, ASTE | **`Paneling`** |
| `/CR12/` | CVPAR, CTERAT, CTRRAT, XSREF1/2, XPREF1/2 | **`Settings::paneling`** |
| `/CI04/` (part) | N, NB, NW, NPAN, IST | `Paneling` (N, NW), `Geometry` (NB), `Settings` (NPAN), `BoundaryLayer` (IST) |
| `/CR06/` (part) | GAM, GAMU, GAM_A, SIG, GAMTE, SIGTE, …_A, SST, SST_GO, SST_GP | **`InviscidSolution`** (SST* → `BoundaryLayer`) |
| `/CR03/` + AIJPIV | AIJ (LU-factored), DIJ | **`InviscidSolution`** (`aij_lu`, `dij`) |
| BIJ, CIJ (EQUIVALENCEd into VM) | dγ/dσ, dQ/dγ scratch | local temporaries of `qdcalc`/`ggcalc`. Checked: no read-after-clobber in the analysis path. |
| `/QMAT/` | Q, DQ (scratch), DZDG, DZDN, DZDM, DQDG, DQDM, QTAN1/2, Z_* | PSILIN output struct `PsiResult` (returned, not global); Q/DQ scratch = locals |
| W1…W8 (EQUIVALENCEd into Q) | work arrays | local `Vec`s |
| `/CR04/` | QINV, QVIS, CPI, CPV, QINVU, QINV_A | **`InviscidSolution`** (QVIS, CPV → `BoundaryLayer` view) |
| `/CR09/` | ADEG, ALFA, AWAKE, AVISC, MVISC, XCMREF, YCMREF, CL, CM, CD, CDP, CDF, CL_ALF, CL_MSQ, PSIO, CIRC, COSA, SINA, QINF, GAMMA, GAMM1, MINF1, MINF, MINF_CL, TKLAM, TKL_MSQ, CPSTAR, QSTAR, CPMN… | inputs → **`Settings`**; α/CL/Mach state + results → **`OperatingPoint`** (**proposed extra struct**, see below) |
| `/CR11/` | PI, HOPI, QOPI, DTOR | `const`s |
| `/CL01/` (solver flags) | LVISC, LALFA, LWAKE, LBLINI, LIPAN, LQAIJ, LADIJ, LWDIJ, LGAMU, LVCONV, SHARP | cache validity → `Option<…>` fields / explicit `valid` flags on the owning struct; LALFA → `OperatingPoint::mode` |
| `/CL01/` (plot/UI flags) | ~50 plotting flags | **not ported** |
| `/CR15/` | XSSI, UEDG, UINV, MASS, THET, DSTR, CTAU, DELT, TSTR, USLP, GUXQ, GUXD, TAU, DIS, CTQ, VTI, UINV_A (IVX×2), REINF1, REINF, REINF_CL, ACRIT, XSTRIP, XOCTR, YOCTR, XSSITR, TINDEX | **`BoundaryLayer`** (per-side `[Side; 2]`); REINF1, ACRIT, XSTRIP → `Settings` |
| `/CI05/`, `/CL02/` | IXBLP, IBLTE, NBL, IPAN, ISYS, NSYS, ITRAN, IDAMP; TFORCE | **`BoundaryLayer`** (IDAMP → `Settings`) |
| `/VMAT/` | VA, VB, VDEL (3×2×IZX), VM (3×IZX×IZX), VZ | **`NewtonSystem`** |
| `/CR17/`, `/CI06/`, `/CC03/` | RMSBL, RMXBL, RLX, VACCEL; IMXBL, ISMXBL; VMXBL | **`NewtonSystem`** (VACCEL → `Settings`) |
| `/BLPAR/` (BLPAR.INC) | SCCON, GACON, GBCON, GCCON, DLCON, CTRCON, CTRCEX, DUXCON, CTCON, CFFAC, **GWAKE**◆ | **`Settings::bl`** (`BlParams`) |
| `/V_VAR1/`, `/V_VAR2/` + EQUIVALENCE COM1/COM2 (XBL.INC, NCOM=73) | station "1" / "2" primary + secondary vars + all partials | `bl::Station` (plain `Copy` struct). `COM1 = COM2` → `s1 = s2` |
| `/V_SAV/` | C1SAV, C2SAV | `Station` copies inside TRDIF |
| `/V_VARA/` | CFM + partials, XT + partials | `bl::Kernel::mid`, `bl::Kernel::xt` |
| `/V_VAR/` | DWTE, QINFBL, TKBL, RSTBL, HSTINV, REYBL (+partials), GAMBL, GM1BL, HVRAT, BULE, XIFORC, AMCRIT | `bl::KernelParams` (built once per SETBL) |
| `/V_INT/` | SIMI, TRAN, TURB, WAKE, TRFORC, TRFREE, IDAMPV | `bl::Kernel` flags |
| `/V_SYS/` | VS1, VS2, VSREZ, VSR, VSM, VSX | `bl::LocalSystem` returned by BLSYS/TESYS |
| `/CI01/`, `/CR07/`, `/CR10/`, `/CI03/`, `/CC01/`, `/CC02/`, `/CR13/`, `/CR18/`, `/CR19/`, `/CI19/`, `/CR01/` | inverse design, polar storage/plots, names, plot layout | not ported. The polar result is `Polar { points: Vec<PolarPoint> }` |
| `/CPI01/…` (CIRCLE.INC), `/WORK/`, `/COM_GUI/`, PPLOT/PXPLOT, `/AICOM_*/` | MDES, GUI, plot programs, OSMAP | not ported |

`bl::Kernel` (the XBL.INC state) is owned by `BoundaryLayer` as a scratch workspace and passed `&mut` explicitly. There is no global state.

**Approved addition:** a seventh struct, `OperatingPoint`, for the /CR09/ state that changes during a solve (α, CL, CM, CD, CDP, CDF, M∞(CL), Re(CL), Karman–Tsien terms, CL_ALF, CL_MSQ). It belongs in neither `Settings` (inputs) nor `InviscidSolution` (which α=0°/90° caching makes α-independent).

## 5. QFoil deviations from XFOIL 6.99

All are marked `C---QFOIL MOD---` in the source. Each was confirmed by `diff -u reference/xfoil-6.99 reference/qfoil`. ✔doc = documented in `version_notes_qfoil.txt`. ✔th = listed on the QBlade theory page.

| # | File:lines (qfoil) | Routine | Change | Doc | Rust location |
|---|---|---|---|---|---|
| D1 | xblsys.f:1766–1800 | BLDIF (turbulent shear-lag eq.) | `SCC = Kc(Hk)·1.333/(1+Us)`, `Kc = 4.65 − 0.95·tanh(0.5·(Hka − 3.5))` replaces `SCCON·1.333/(1+Us)`. New `SCC_HKA = −0.95·0.5·sech²(·)·1.333/(1+Us)`. | ✔doc ✔th | `bl::equations::bldif` |
| D2 | xblsys.f:1825–1831 | BLDIF | Jacobian: `Z_HKA += SCC_HKA·(CQA − SA·ALD)·DXI`. Flows into VS1/VS2 via Z_HK1/Z_HK2 and Z_UPW. | ✔doc ✔th | same; FD test required |
| D3 | xblsys.f:1699–1713 | BLDIF | wake (ITYP=3): `ALD = DLCON/4.0` (was `DLCON`) | ✔doc ✔th | same |
| D4 | xblsys.f:2442–2445 | HST | `HSMIN, DHSINF = 1.505, 0.04` (was 1.500, 0.015) | ✗doc ✔th | `bl::closure::hst` |
| D5 | xbl.f:1072–1075 | MRCHDU | station Newton `RLX = 0.7` (was 1.0), then `IF(DMAX>0.3) RLX=0.3/DMAX` (see S1) | ✔doc (as "RLX 0.7") | `bl::march::mrchdu` |
| D6 | xbl.f:1384–1387 | UPDATE | global Newton `RLX = 0.7` initial (was 1.0) | ✔doc ✔th | `newton::update` |
| D7 | xbl.f:1279–1285 | UPDATE | `DCLMAX/DCLMIN = ±0.15` (was ±0.5). **Removed** `IF(MATYP≠1) DCLMIN = MAX(−0.5, −0.9·CL)` | ✔doc (limits) ✗doc (guard removal) | `newton::update` |
| D8 | xbl.f:1414–1419 | UPDATE | `DHI = 1.0, DLO = −0.4` (was 1.5, −0.5) | ✗doc ✔th | `newton::update` |
| D9 | xoper.f:2951–2991 | VISCAL | always: `UEDG = (UINV > 2 ? 2 : UINV)`, `CTAU = 0.01` for IBL=1..NBL, `LBLINI = .FALSE.` (was: UEDG=UINV only if !LBLINI) | ✔doc ✔th | `solver::viscal` |
| D10 | xfoil.f:1191–1193, 1208–1231 | CDCALC | if `u = Ue_wake/Q∞ ≤ 1`: `H1 = 3.15 + 1.72/(H−1)`, `θ_tot = θ·[1 + (1−u)(u(GWAKE·H1 − 1) − 1)]`, `CD = 2θ_tot·u^((5+H)/2)`; else original Squire–Young | ✔doc ✔th | `wake::drag::cdcalc` |
| D11 | BLPAR.INC:12–16; xbl.f:1601–1603 | /BLPAR/, BLPINI | new `GWAKE`, default 0.40 | ✔doc ✔th | `BlParams::gwake` |
| D12 | xoper.f:2516–2535, 2567–2569, 2721–2730 | VPAR | display GWAKE; `GW r` command | ✔doc | `Settings` + CLI `--gwake` |
| D13 | XFOIL.INC:23–38 | — | `IQX=1400, NAX=1200` (upstream 370/800). Also present in bundled baseline. **Numerical side effect:** `NACA` uses `NSIDE = IQX/3` points per side, so QFoil's NACA buffer has 931 points (stock XFOIL: 245), and the PANGEN input and thus the panels differ from stock XFOIL. | ✔th | Rust uses `Vec`s; `limits::IQX` reproduces NSIDE and the input limits (ABCOPY ≤ IQX−5, PANGEN ≤ IQX−1) as errors |
| D14 | xfoil.f:40–63; xoper.f:3109–3145 | banner; VISCAL `.bl` dump | cosmetic / disabled debug output | ✔doc | not ported |

**Behavioural consequence of D9 (important for Phase 4).** Each α starts cold: `MRCHUE` re-initialises θ, δ* and Ctau from Thwaites at the stagnation point. What still carries over from the previous α in a sweep:
1. `LIPAN` stays `.TRUE.`, so `STFIND`/`IBLPAN`/`XICALC`/`IBLSYS` are **not** re-run at the start of `VISCAL`. The stagnation panel `IST`, `SST`, the `IPAN`/`ISYS` maps, `NBL`, `IBLTE` and `XSSI` are the **previous α's converged (`STMOVE`d) values**. `UINV`/`UEDG` initial values are indexed through them.
2. The sweep-halting counter (`NSEQEX`).
3. Any SAVEd locals (S9).

So in principle results depend on sweep order through the stagnation-point indexing. **Phase 1 showed this dependency is fatal rather than weak (S10):** in-session sweeps go NaN from the second point. In QBlade's usage (a fresh process per α) every point is independent. The Phase 4 requirement "preserve the warm start from the previous α within each sweep direction" therefore has nothing to preserve in QFoil-parity mode, and α points parallelise trivially. **Confirmed:** independent points; Phase 4 parallelises over α with no warm start.

## 6. Porting order

Bottom-up, as specified. One adjustment: **XYWAKE, QWCALC and QDCALC move into step 3**, because the source-influence matrix DIJ covers wake nodes and the step-3 golden dumps (DIJ) need them. Step 7 keeps the wake *BL* and coupling pieces: TESYS/wake closures parity, STFIND…STMOVE, UESET, CDCALC and GWAKE.

1. Splines and geometry: `spline`, `geometry`, `geometry::naca`, `geometry::dat`
2. Panelling: `paneling`
3. Inviscid: `linalg`, `inviscid::influence`, `inviscid` (GGCALC, QISET, SPECAL, CLCALC, CPCALC), `wake::trajectory` (XYWAKE), QWCALC, QDCALC
4. BL closures and derivatives: `bl::closure`, `bl::transition::{dampl, dampl2, axset}` (+ FD derivative tests)
5. BL march and transition: `bl::station`, `bl::equations`, `bl::transition::trchek`, `bl::march`
6. Newton system: `newton::{assemble, solve, update}`
7. Wake/coupling: `coupling`, `wake::drag`
8. Driver: `solver::viscal`, `polar::sweep`, `settings`, then `qfoil-cli`

### Parity strategy

* Floating-point operation order is kept, including the order of sums in the influence loops. gfortran at `-O2` without `-ffast-math` or `-march` does not reassociate and emits no FMA on generic x86-64, so for arithmetic-only routines Rust `f64` can be **bit-identical**.
* Transcendentals: gfortran and Rust `std` both call glibc `libm` on Linux (`exp`, `log`, `tanh`, `sqrt`, `pow`), so they should match. Fortran `X**2`/`X**3` with integer exponents is expanded by gfortran to multiplications, so it is ported as explicit `x*x`, `x*x*x` rather than `powi`. `X**0.5` stays `powf(0.5)` rather than `sqrt`. Real exponents stay `powf`.
* Golden dumps are written with `ES24.16E3` (17 significant digits) so bit-exactness can be checked, not just closeness.

## 7. Suspected issues (ported faithfully; your decision needed where marked)

| # | Location | Observation | Proposed handling |
|---|---|---|---|
| S1 | xbl.f:1074–1076 (MRCHDU) | `RLX = 0.7` then `IF(DMAX.GT.0.3) RLX = 0.3/DMAX`. For 0.3 < DMAX < 0.4286 this **raises** RLX above 0.7 (e.g. DMAX=0.35 → 0.857). Probably meant `MIN(0.7, 0.3/DMAX)`. | **Decided: port as is.** Capping was tested (`tools/patches/experiments/s1/`, `tools/experiment_compare.sh qfoil-s1cap`). Converged points go 740 → 741 over the 915-point matrix (+22/−21, scattered), and max \|ΔCL\| is about 2e-5 on points converged in both, so no benefit. |
| S2 | xblsys.f:1770–1776 | Comment says `Kc = 4.65 − 0.95·tanh(0.275·Hk − 3.5)`; code uses `0.5·(Hk − 3.5)`. The theory page and version notes agree with the code. | Port the code. Comment is stale; noted in the rustdoc. |
| S3 | xoper.f:2982–2984 | Comment: "Reset amplification factor … 0.01". `CTAU` holds N (laminar) / Ctau (turbulent). Since `LBLINI=F` forces `MRCHUE`, which overwrites `CTAU(2..NBL)` before any read, the reset has **no effect** on results. | Port it anyway (harmless) with a comment. |
| S4 | xbl.f:1284 (UPDATE) | Removal of the `MATYP≠1` DCLMIN guard is undocumented. It has no effect for fixed-Mach polars (MATYP=1, used here). | Port. Note for Type 2 polars. |
| S5 | xoper.f + BLPAR | `SCCON` is now unused (D1 hard-codes 4.65/0.95). `VPAR LAG` still edits it and the VPAR display still shows it as "Klag". | Keep `sccon` in `BlParams` for format parity, document as unused. |
| S6 | bin/Makefile comment | The comment says `-fdefault-real-8` was removed (single precision), but the flag **is present**. The binary is double precision. Without `-fdefault-double-8`, the `D0` literals in VISCAL (D9) are promoted to REAL(16). 2.0 and 0.01 round back to the same f64, so there is no effect. | Target f64. Verify by bit-comparing Ue after VISCAL init. |
| S7 | xfoil.f:1218 (CDCALC) | `H1 = 3.15 + 1.72/(H−1)` (Head entrainment, truncated form) is singular as H→1. The wake-end H ≥ 1.00005 via DSLIM, so there's no overflow, but CD becomes very sensitive for H near 1. Continuous at u=1 (factor → 0). | Port. Monitor in validation. |
| S8 | xbl.f:1419 | Closing marker reads `C---QFOIL MOD---` instead of `…END…`. Cosmetic. | None. |
| S10 | xoper.f:2951–2991 (D9) × `LIPAN` | **In-process α sweeps fail.** After the first converged point, every later `ALFA`/`ASEQ` point in the same session diverges to NaN (MRCHUE similarity station, side 2). D9 forces a cold MRCHUE start, but `LIPAN` stays true, so `IST`/`IPAN` are the previous α's. At the new α the first station next to the old stagnation point sees UINV ≤ 0, and Thwaites' `SQRT` goes NaN. Verified: α=1 converges in a fresh process and after `INIT` (which clears `LIPAN`), and fails straight after α=0. It reproduces with QBlade's shipped `QFoil` binary. QBlade is unaffected because it launches **one QFoil process per α** (its command template, recovered from `libQBladeCE`, has a single `ALFA` and no `ASEQ`/`PACC`). | Golden data = one fresh process per α. The Rust default `analyse_polar` matches QBlade (independent points, recompute stagnation indexing per α). **Decided:** independent points; no session-sweep mode. |
| S11 | xfoil.f / Makefile | The shipped `libPlt.a` is a MinGW (Windows) object and `dummy_plot.f` is not valid fixed-form (line 1 unindented, line 94 > 72 columns), so the Linux build needs the build-only fixes in `tools/build_reference.sh`. `ASEQ` in the shipped Linux `QFoil` segfaults (plot set-up path). | Build-only, no numerical effect. |
| S12 | aread/LOAD path | QFoil truncates long `LOAD` path arguments (an absolute path became `/`). | Harness copies airfoils into the run directory. Not relevant to the port. |
| S13 | xsolve.f:174 (LUDCMP), spline.f:23/65 | `NVX=500` in LUDCMP and `NMAX=1000` in SPLINE/SPLIND cap QFoil at **499 panel nodes** and 1000 points per spline segment, although `IQX=1400` suggests ~1395 (D13). QFoil `STOP`s beyond that. | Rust returns `LinalgError::TooLarge` for N+1 > 500. The spline has no limit (it can't change results for valid input). |
| S14 | xblsys.f:2365 (DILW), **upstream XFOIL 6.99** (identical in the MIT release) | Sign error in the analytic derivative of the laminar wake dissipation: `RCD_HK = -1.10*(1-1/HK)*2/HK**3 - RCD/HK`, but d/dHK[1.1(1−1/HK)²/HK] = **+**2.2(1−1/HK)/HK³ − RCD/HK. At HK = 1.3 the Jacobian entry is −0.27 against a true +0.20. Only the Newton Jacobian is affected (laminar wake via BLVAR, xblsys.f:1072), so it can slow convergence but cannot change a converged solution. | Ported faithfully (bit parity); `dilw_hk_derivative_sign_error` pins it. **Ask:** fix after parity, as an opt-in or as the default? |
| S15 | xblsys.f BLDIF shear-lag row, **upstream XFOIL 6.99** | The analytic Jacobian of the shear-lag equation omits the Rθ dependence of `HKC = Hk − 1 − GCCON/Rθ`: `UQ_T1, UQ_U1, …` are formed but never used, so the θ and Ue entries are about 0.6% off at a typical turbulent station. The Jacobian is exact when GCCON = 0 (verified). Only the convergence rate is affected. QFoil's added `SCC_HKA` term (D2) is verified correct. | Ported faithfully; `station_jacobians_match_finite_differences` pins it. Same decision as S14. |
| S16 | SETBL/BLPRV `HVRAT`, **upstream** | The Sutherland ratio `HVRAT` is set only by plotting routines (`blplot.f`, `dplot.f`), so it is 0 in every analysis. Effect only for M > 0. | Ported as 0 (`KernelParams::hvrat`). No action. |
| S9 | whole build | `-fno-automatic -finit-local-zero`: any routine that reads a local before assigning it in the current call gets the value from its previous call (or 0 on the first call). A faithful Rust port with fresh locals would then differ. | Probed in Phase 1 with a `-fautomatic -finit-real=snan -finit-integer=-99999999` build. Result is in `tests/golden/checks.tsv` (`dump_eq_auto_final`) and summarised in `docs/PROGRESS.md`. |

## 8. Golden data (Phase 1)

Implemented. See `tests/golden/README.md` for layout and formats. The bullets below are the original plan, amended.

* `reference-build/qfoil/`: out-of-tree build of the unmodified QFoil source (copied, never edited in `reference/`). Plus `reference-build/qfoil-dump/`: a patched copy with dump statements. Patches are kept as `tools/patches/*.patch` so they are reviewable.
* Driver: `tools/gen_golden.sh`, which feeds QFoil command scripts via stdin in a clean temporary directory (no `xfoil.def`).
* ~~Sweeps from 0° in each direction~~. **Amended (S10):** every α ∈ [−10, 20], step 0.5, is a fresh QFoil process with one `ALFA` command, which is how QBlade drives QFoil. Airfoils from files are `LOAD`ed, so their raw coordinates are the panel nodes (ABCOPY, as in QBlade). NACA sections use the built-in `NACA` command (PANGEN, 160 nodes). `ITER 100`, Ncrit 9, M = 0.
* Cases: NACA 0012, NACA 0020, NACA 4412 (built-in `NACA` command), DU 91-W2-250 (**no local coordinates, see open question**), and a cambered low-Re section (proposal: **E387**, which ships in `reference/qfoil/runs/e387.dat`) × Re {1e5, 1e6, 5e6} × Ncrit 9, M = 0.
* Captured: polar (α, CL, CD, CDp, CM, Xtr top/bot), convergence flag per α (from stdout, since the polar file only holds converged points), Cp and BL dumps (`CPWR`, `DUMP`; Ctau via patched dump) at α ∈ {0, 5, 10, 15}. Intermediate dumps (AIJ, γ, DIJ, BL arrays per Newton iteration) come from the patched build.
