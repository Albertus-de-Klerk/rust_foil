//! Step 4 parity: BL closures and amplification rates against the Fortran reference table
//! (`tools/fdrivers/closures.f`), plus finite-difference checks of every analytic partial.

use qfoil_core::bl::closure::{Correlation, cfl, cft, dil, dilw, dit, hct, hkin, hsl, hst};
use qfoil_core::bl::transition::{
    AmpStation, Amplification, AmplificationModel, IntervalRate, axset, dampl, dampl2,
};
use qfoil_golden::{assert_derivatives, assert_golden, closure_table};

const CFFAC: f64 = 1.0;

fn corr(c: Correlation) -> Vec<f64> {
    vec![c.v, c.hk, c.rt, c.msq]
}
fn amp(a: Amplification) -> Vec<f64> {
    vec![a.ax, a.hk, a.th, a.rt]
}

/// Every row of the table must be reproduced bit for bit.
#[test]
fn closures_bit_identical_to_reference() {
    let table = closure_table();
    let mut checked = 0;
    for (name, rows) in &table {
        for r in rows {
            let got: Vec<f64> = match name.as_str() {
                "HKIN" => {
                    let (a, b, c) = hkin(r[0], r[1]);
                    vec![a, b, c]
                }
                "HCT" => {
                    let (a, b, c) = hct(r[0], r[1]);
                    vec![a, b, c]
                }
                "DIL" => corr(dil(r[0], r[1]))[..3].to_vec(),
                "DILW" => corr(dilw(r[0], r[1]))[..3].to_vec(),
                "HSL" => corr(hsl(r[0], r[1], r[2])),
                "HST" => corr(hst(r[0], r[1], r[2])),
                "CFL" => corr(cfl(r[0], r[1], r[2])),
                "CFT" => corr(cft(r[0], r[1], r[2], CFFAC)),
                "DAMPL" => amp(dampl(r[0], r[1], r[2])),
                "DAMPL2" => amp(dampl2(r[0], r[1], r[2])),
                "DIT" => {
                    let (a, b, c, d, e) = dit(r[0], r[1], r[2], r[3]);
                    vec![a, b, c, d, e]
                }
                "AXSET" => {
                    let model = if r[0] == 0.0 {
                        AmplificationModel::Envelope
                    } else {
                        AmplificationModel::Modified
                    };
                    let s1 = AmpStation {
                        hk: r[1],
                        th: r[2],
                        rt: r[3],
                        ampl: r[4],
                    };
                    let s2 = AmpStation {
                        hk: r[5],
                        th: r[6],
                        rt: r[7],
                        ampl: r[8],
                    };
                    let IntervalRate { ax, d1, d2 } = axset(s1, s2, 9.0, model);
                    [vec![ax], d1.to_vec(), d2.to_vec()].concat()
                }
                other => panic!("unexpected routine {other}"),
            };
            let n_in = r.len() - got.len();
            let label = format!("{name}{:?}", &r[..n_in]);
            let want = r[n_in..].to_vec();
            let result = std::panic::catch_unwind(|| assert_golden!(ulps = 0; "out" => got, want));
            assert!(result.is_ok(), "{label}");
            checked += 1;
        }
    }
    assert_eq!(checked, 6331, "all table rows checked");
}

// ---- analytic derivatives vs finite differences (points inside each branch) ----

const HK_POINTS: [f64; 9] = [1.3, 2.2, 2.9, 3.7, 4.2, 4.8, 5.8, 7.5, 12.0];
const RT_POINTS: [f64; 5] = [120.0, 300.0, 900.0, 2.0e4, 3.0e5];

#[test]
fn laminar_and_turbulent_correlation_derivatives() {
    let v = |c: &Correlation| c.v;
    for hk in HK_POINTS {
        for rt in RT_POINTS {
            for msq in [0.0, 0.2] {
                let x = [hk, rt, msq];
                assert_derivatives!(|x: [f64; 3]| hst(x[0], x[1], x[2]), at x, value = v, rel = 1e-5,
                    { 0 => |c: &Correlation| c.hk, 1 => |c: &Correlation| c.rt, 2 => |c: &Correlation| c.msq });
                assert_derivatives!(|x: [f64; 3]| cft(x[0], x[1], x[2], CFFAC), at x, value = v, rel = 1e-5,
                    { 0 => |c: &Correlation| c.hk, 1 => |c: &Correlation| c.rt, 2 => |c: &Correlation| c.msq });
                assert_derivatives!(|x: [f64; 3]| hsl(x[0], x[1], x[2]), at x, value = v, rel = 1e-5,
                    { 0 => |c: &Correlation| c.hk });
                assert_derivatives!(|x: [f64; 3]| cfl(x[0], x[1], x[2]), at x, value = v, rel = 1e-5,
                    { 0 => |c: &Correlation| c.hk, 1 => |c: &Correlation| c.rt });
            }
            let x = [hk, rt];
            assert_derivatives!(|x: [f64; 2]| dil(x[0], x[1]), at x, value = v, rel = 1e-5,
                { 0 => |c: &Correlation| c.hk, 1 => |c: &Correlation| c.rt });
            // ∂/∂HK of DILW is wrong in the reference (S14); see dilw_hk_derivative_sign_error
            assert_derivatives!(|x: [f64; 2]| dilw(x[0], x[1]), at x, value = v, rel = 1e-5,
                { 1 => |c: &Correlation| c.rt });
        }
        for msq in [0.0, 0.1, 0.3] {
            let x = [hk, msq];
            assert_derivatives!(|x: [f64; 2]| hkin(x[0], x[1]), at x, value = |r: &(f64, f64, f64)| r.0, rel = 1e-6,
                { 0 => |r: &(f64, f64, f64)| r.1, 1 => |r: &(f64, f64, f64)| r.2 });
            assert_derivatives!(|x: [f64; 2]| hct(x[0], x[1]), at x, value = |r: &(f64, f64, f64)| r.0, rel = 1e-6,
                { 0 => |r: &(f64, f64, f64)| r.1, 1 => |r: &(f64, f64, f64)| r.2 });
        }
    }
    let x = [1.6, 0.6, 0.002, 0.12];
    type D5 = (f64, f64, f64, f64, f64);
    assert_derivatives!(|x: [f64; 4]| dit(x[0], x[1], x[2], x[3]), at x, value = |r: &D5| r.0, rel = 1e-6,
        { 0 => |r: &D5| r.1, 1 => |r: &D5| r.2, 2 => |r: &D5| r.3, 3 => |r: &D5| r.4 });
}

#[test]
fn amplification_rate_derivatives() {
    let v = |a: &Amplification| a.ax;
    // Rθ well above critical (ramp saturated) and inside the ramp
    for (hk, rt) in [
        (2.3, 2000.0),
        (2.6, 800.0),
        (3.2, 400.0),
        (3.8, 300.0),
        (5.0, 500.0),
        (8.0, 1.0e4),
    ] {
        let x = [hk, 1.0e-3, rt];
        for f in [dampl, dampl2] {
            assert_derivatives!(|x: [f64; 3]| f(x[0], x[1], x[2]), at x, value = v, rel = 1e-5, scale = 1.0,
                { 0 => |a: &Amplification| a.hk, 1 => |a: &Amplification| a.th, 2 => |a: &Amplification| a.rt });
        }
    }
}

#[test]
fn axset_derivatives() {
    for model in [AmplificationModel::Envelope, AmplificationModel::Modified] {
        for (hk, rt, a) in [(2.5, 1500.0, 3.0), (3.9, 600.0, 8.9), (4.4, 900.0, 8.95)] {
            let x = [hk, 1.0e-3, rt, a, hk + 0.1, 1.2e-3, rt * 1.1, a + 0.05];
            let f = |x: [f64; 8]| {
                let s1 = AmpStation {
                    hk: x[0],
                    th: x[1],
                    rt: x[2],
                    ampl: x[3],
                };
                let s2 = AmpStation {
                    hk: x[4],
                    th: x[5],
                    rt: x[6],
                    ampl: x[7],
                };
                axset(s1, s2, 9.0, model)
            };
            assert_derivatives!(f, at x, value = |r: &IntervalRate| r.ax, rel = 1e-5, scale = 1.0, {
                0 => |r: &IntervalRate| r.d1[0], 1 => |r: &IntervalRate| r.d1[1],
                2 => |r: &IntervalRate| r.d1[2], 3 => |r: &IntervalRate| r.d1[3],
                4 => |r: &IntervalRate| r.d2[0], 5 => |r: &IntervalRate| r.d2[1],
                6 => |r: &IntervalRate| r.d2[2], 7 => |r: &IntervalRate| r.d2[3],
            });
        }
    }
}

/// PORTING_PLAN S14: XFOIL 6.99 `DILW` (xblsys.f:2365) has
/// `RCD_HK = -1.10*(1-1/HK)*2/HK**3 - RCD/HK`, but the true derivative of
/// `1.1*(1-1/HK)**2/HK` is `+2.2*(1-1/HK)/HK**3 - RCD/HK`. The port keeps the reference
/// Jacobian. This test pins that behaviour and records the size of the error; it fails if
/// the derivative is ever corrected without a decision.
#[test]
fn dilw_hk_derivative_sign_error() {
    let (hk, rt) = (1.3, 1000.0);
    let d = dilw(hk, rt);
    let fd = qfoil_golden::fd::central_difference(|x: [f64; 2]| dilw(x[0], x[1]).v, [hk, rt], 0);
    // reference Jacobian, reconstructed with the wrong sign:
    let hs = hsl(hk, rt, 0.0);
    let rcd = 1.10 * (1.0 - 1.0 / hk).powi(2) / hk;
    let rcd_hk_ref = -1.10 * (1.0 - 1.0 / hk) * 2.0 / (hk * hk * hk) - rcd / hk;
    let di_hk_ref = 2.0 * rcd_hk_ref / (hs.v * rt) - (d.v / hs.v) * hs.hk;
    assert_eq!(d.hk, di_hk_ref, "port keeps the reference derivative");
    assert!(
        (d.hk - fd).abs() > 0.5 * fd.abs(),
        "reference derivative differs from the true one"
    );
}
