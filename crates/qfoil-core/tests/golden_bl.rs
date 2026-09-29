//! Step 5 parity: BL march (MRCHUE, MRCHDU) with the full station kernel (BLPRV, BLKIN,
//! BLVAR, BLMID, BLDIF incl. QFoil shear lag, TRCHEK2, TRDIF, BLSYS, TESYS).

mod common;

use qfoil_core::bl::march::MarchInputs;
use qfoil_core::bl::station::{Kernel, KernelParams};
use qfoil_core::bl::transition::AmplificationModel;
use qfoil_core::bl::{BoundaryLayer, Side};
use qfoil_core::inviscid::forces::Compressibility;
use qfoil_core::settings::{BlParams, FlowConditions};
use qfoil_golden::{Dump, DumpSet, assert_golden, fixture, full_dump};

/// BL state as dumped by the reference (`viscal_init`, `mrchue`, `setbl`, `update`, ...).
fn bl_from_dump(d: &Dump) -> BoundaryLayer {
    let nbl = d.ints("NBL");
    let iblte = d.ints("IBLTE");
    let itran = d.ints("ITRAN");
    let side = |k: &str| {
        let r = |name: &str| d.reals(&format!("{name}{k}")).to_vec();
        // station 1 (the dummy) has IPAN = 0 in the Fortran; keep it at 0
        let idx = |name: &str| {
            d.ints(&format!("{name}{k}"))
                .iter()
                .map(|&v| (v - 1).max(0) as usize)
                .collect()
        };
        Side {
            xssi: r("XSSI"),
            uedg: r("UEDG"),
            uinv: r("UINV"),
            uinv_a: r("UINV_A"),
            mass: r("MASS"),
            thet: r("THET"),
            dstr: r("DSTR"),
            ctau: r("CTAU"),
            delt: r("DELT"),
            tstr: r("TSTR"),
            tau: r("TAU"),
            dis: r("DIS"),
            ctq: r("CTQ"),
            uslp: vec![0.0; r("CTQ").len()],
            vti: r("VTI"),
            ipan: idx("IPAN"),
            isys: idx("ISYS"),
        }
    };
    BoundaryLayer {
        sides: [side("1"), side("2")],
        iblte: [iblte[0] as usize - 1, iblte[1] as usize - 1],
        nbl: [nbl[0] as usize, nbl[1] as usize],
        itran: [itran[0].max(1) as usize - 1, itran[1].max(1) as usize - 1],
        xssitr: [d.reals("XSSITR")[0], d.reals("XSSITR")[1]],
        tforce: [false, false],
        ist: d.int("IST") as usize - 1,
        sst: d.real("SST"),
        sst_go: 0.0,
        sst_gp: 0.0,
        nsys: d.int("NSYS") as usize,
        wgap: d.reals("WGAP").to_vec(),
        xoctr: [d.reals("XOCTR")[0], d.reals("XOCTR")[1]],
        yoctr: [0.0, 0.0],
    }
}

/// Kernel constants exactly as SETBL sets them, for lift coefficient `cl`.
fn kernel(run: &str, cl: f64, wgap: &[f64]) -> Kernel {
    let flow = FlowConditions {
        reynolds: common::reynolds(run),
        ..FlowConditions::default()
    };
    let comp = Compressibility::at_cl(cl, &flow, 1.0);
    let p = KernelParams::new(
        comp.minf,
        comp.reinf,
        comp.tklam,
        comp.tkl_msq,
        flow.gamma,
        1.0,
        wgap[0],
    );
    Kernel::new(p, BlParams::default(), AmplificationModel::Envelope)
}

/// Which arrays a reference dump can verify.
#[derive(Clone, Copy, PartialEq)]
enum Arrays {
    /// Primary and derived arrays (the `mrchue` dump).
    All,
    /// Primary arrays only: SETBL recomputes TAU, DIS, CTQ, DELT, TSTR (xbl.f:278) and
    /// XSSITR (xbl.f:445) in its own station loop before the `setbl` dump, so those belong
    /// to the SETBL test.
    Primary,
}

fn compare_bl(what: &str, bl: &BoundaryLayer, d: &Dump, arrays: Arrays, ulps: u64) {
    let mut failures = Vec::new();
    for (is, k) in [(0, "1"), (1, "2")] {
        let s = &bl.sides[is];
        let n = d.reals(&format!("THET{k}")).len();
        let tag = |name: &str| format!("{what} {name}{k}");
        let mut list = vec![
            ("THET", &s.thet),
            ("DSTR", &s.dstr),
            ("CTAU", &s.ctau),
            ("UEDG", &s.uedg),
            ("MASS", &s.mass),
        ];
        if arrays == Arrays::All {
            list.extend([
                ("TAU", &s.tau),
                ("DIS", &s.dis),
                ("CTQ", &s.ctq),
                ("DELT", &s.delt),
                ("TSTR", &s.tstr),
            ]);
        }
        for (name, got) in list {
            let want = d.reals(&format!("{name}{k}"));
            let r = qfoil_golden::compare_slices(&tag(name), &got[1..n], &want[1..n], ulps, 0.0);
            if let Err(m) = r {
                failures.push(m.to_string());
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    let itran: Vec<i64> = bl.itran.iter().map(|&v| v as i64 + 1).collect();
    assert_eq!(itran, d.ints("ITRAN"), "{what} ITRAN");
    if arrays == Arrays::All {
        // SETBL also rewrites XSSITR (xbl.f:445)
        assert_golden!(ulps = ulps; "XSSITR" => bl.xssitr.to_vec(), d.reals("XSSITR"));
    }
}

fn check_run(run: &str, set: &DumpSet, ulps: u64) {
    let (_, pan) = common::setup(run);
    let init = set.first("viscal_init");
    let mut bl = bl_from_dump(init);
    // SETBL uses the CL current at that point, which is the inviscid CL
    let mut k = kernel(run, set.first("specal").real("CL"), init.reals("WGAP"));
    let inp = MarchInputs {
        pan: &pan,
        acrit: [9.0, 9.0],
        xstrip: [1.0, 1.0],
    };

    bl.mrchue(&mut k, &inp);
    compare_bl("MRCHUE", &bl, set.first("mrchue"), Arrays::All, ulps);

    bl.mrchdu(&mut k, &inp);
    compare_bl("MRCHDU", &bl, set.first("setbl"), Arrays::Primary, ulps);
}

#[test]
fn naca0012_march() {
    let run = "naca0012_re1e6_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
fn naca4412_march() {
    let run = "naca4412_re1e6_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
fn e387_march() {
    let run = "e387_re1e5_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
#[ignore = "needs full dumps: tools/gen_golden.sh"]
fn full_dump_marches() {
    for run in [
        "du91w2250_re1e6_a10.0",
        "naca0012_re1e6_a0.0",
        "naca0012_re1e6_a15.0",
        "naca4412_re1e6_a15.0",
    ] {
        check_run(run, &full_dump(run), 0);
    }
}

// ---- analytic station Jacobians (VS1, VS2 from BLSYS) against finite differences ----

use qfoil_core::bl::station::Regime;

/// Residuals (−VSREZ rows 1–3) of the interval between stations with primary variables
/// `v1`, `v2` = [N or √Cτ, θ, δ*, Ue, ξ], plus the analytic VS1/VS2 at that point.
fn interval(
    k0: &Kernel,
    regime: Regime,
    v1: [f64; 5],
    v2: [f64; 5],
) -> ([f64; 3], [[f64; 5]; 4], [[f64; 5]; 4]) {
    let mut k = k0.clone();
    k.simi = false;
    k.tran = false;
    k.turb = regime != Regime::Laminar;
    k.wake = regime == Regime::Wake;
    let set = |k: &mut Kernel, v: [f64; 5]| {
        let (ampl, ctau) = if regime == Regime::Laminar {
            (v[0], 0.0)
        } else {
            (0.0, v[0])
        };
        k.blprv(v[4], ampl, ctau, v[1], v[2], 0.0, v[3]);
        k.blkin();
    };
    set(&mut k, v1);
    k.blvar(regime);
    k.s1 = k.s2;
    set(&mut k, v2);
    k.blsys();
    let r = [-k.sys.vsrez[0], -k.sys.vsrez[1], -k.sys.vsrez[2]];
    (r, k.sys.vs1, k.sys.vs2)
}

fn check_jacobian(
    label: &str,
    k: &Kernel,
    regime: Regime,
    v1: [f64; 5],
    v2: [f64; 5],
    skip: &[(usize, usize, usize)],
) {
    let (_, vs1, vs2) = interval(k, regime, v1, v2);
    let mut failures = Vec::new();
    for (st, vs) in [(1usize, vs1), (2, vs2)] {
        for l in 0..5 {
            let base = if st == 1 { v1 } else { v2 };
            let h = 1e-7 * base[l].abs().max(1e-6);
            let eval = |dv: f64| {
                let (mut a, mut b) = (v1, v2);
                if st == 1 {
                    a[l] += dv
                } else {
                    b[l] += dv
                }
                interval(k, regime, a, b).0
            };
            let (rp, rm) = (eval(h), eval(-h));
            for row in 0..3 {
                if skip.contains(&(st, row, l)) {
                    continue;
                }
                let fd = (rp[row] - rm[row]) / (2.0 * h);
                let an = vs[row][l];
                let scale = vs[row].iter().fold(0.0_f64, |m, v| m.max(v.abs()));
                if (an - fd).abs() > 1e-5 * an.abs().max(fd.abs()).max(1e-3 * scale) {
                    failures.push(format!("{label}: d(row {row})/d(station {st} var {l}): analytic {an:e} vs fd {fd:e}"));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// Primary variables of a marched station.
fn vars(bl: &BoundaryLayer, is: usize, ibl: usize) -> [f64; 5] {
    let s = &bl.sides[is];
    [
        s.ctau[ibl],
        s.thet[ibl],
        s.dstr[ibl],
        s.uedg[ibl],
        s.xssi[ibl],
    ]
}

#[test]
fn station_jacobians_match_finite_differences() {
    let run = "naca0012_re1e6_a5.0";
    let set = fixture(run);
    let (_, pan) = common::setup(run);
    let init = set.first("viscal_init");
    let mut bl = bl_from_dump(init);
    let mut k = kernel(run, set.first("specal").real("CL"), init.reals("WGAP"));
    let inp = MarchInputs {
        pan: &pan,
        acrit: [9.0, 9.0],
        xstrip: [1.0, 1.0],
    };
    bl.mrchue(&mut k, &inp);
    k.p.amcrit = 9.0;

    // laminar (upper surface, ahead of transition)
    let i = bl.itran[0] / 2;
    check_jacobian(
        "laminar",
        &k,
        Regime::Laminar,
        vars(&bl, 0, i),
        vars(&bl, 0, i + 1),
        &[],
    );
    // wake (lower side, a few stations behind the TE)
    let w = bl.iblte[1] + 6;
    check_jacobian(
        "wake",
        &k,
        Regime::Wake,
        vars(&bl, 1, w),
        vars(&bl, 1, w + 1),
        &[],
    );
    // turbulent (upper surface, between transition and TE)
    let t = (bl.itran[0] + bl.iblte[0]) / 2;
    // PORTING_PLAN S15 (upstream XFOIL): the shear-lag row omits the Rθ dependence of
    // HKC = Hk - 1 - GCCON/Rθ (UQ_T1.., UQ_U1.. are formed in BLDIF but never used), so its
    // θ and Ue entries are approximate. Checked exactly with GCCON = 0 below.
    let approx = [(1, 0, 1), (1, 0, 3), (2, 0, 1), (2, 0, 3)];
    check_jacobian(
        "turbulent",
        &k,
        Regime::Turbulent,
        vars(&bl, 0, t),
        vars(&bl, 0, t + 1),
        &approx,
    );
    let mut k0 = k.clone();
    k0.bl.gccon = 0.0;
    check_jacobian(
        "turbulent, GCCON = 0",
        &k0,
        Regime::Turbulent,
        vars(&bl, 0, t),
        vars(&bl, 0, t + 1),
        &[],
    );
}
