//! Step 3 parity: influence matrices, linear solve, inviscid solution, wake, DIJ.

mod common;

use qfoil_core::inviscid::GgcalcSystem;
use qfoil_core::settings::FlowConditions;
use qfoil_core::{InviscidSolution, Settings, Wake};
use qfoil_golden::{DumpSet, assert_golden, fixture, full_dump};

fn settings(run: &str) -> Settings {
    Settings {
        flow: FlowConditions {
            reynolds: common::reynolds(run),
            ..FlowConditions::default()
        },
        ..Settings::default()
    }
}

fn check_run(run: &str, set: &DumpSet, ulps: u64) {
    let (_, pan) = common::setup(run);
    let n = pan.len();

    // GGCALC assembly
    let sys = GgcalcSystem::assemble(&pan);
    let d = set.first("aij_raw");
    assert_golden!(ulps = ulps;
        "AIJ" => sys.aij.as_slice(), d.matrix("AIJ").2;
        "BIJ" => sys.bij.as_slice(), d.matrix("BIJ").2;
        "GAMU0" => sys.gamu[0], &d.matrix("GAMU").2[..n + 1];
        "GAMU1" => sys.gamu[1], &d.matrix("GAMU").2[n + 1..]);

    // GGCALC factor + solve
    let mut sol = InviscidSolution::new(&pan).unwrap();
    let d = set.first("ggcalc");
    let piv: Vec<i64> = sol.aij.pivots.iter().map(|&p| p as i64 + 1).collect();
    assert_eq!(piv, d.ints("AIJPIV"), "AIJPIV");
    assert_golden!(ulps = ulps;
        "AIJLU" => sol.aij.lu.as_slice(), d.matrix("AIJLU").2;
        "GAMU0" => sol.gamu[0], &d.matrix("GAMU").2[..n + 1];
        "GAMU1" => sol.gamu[1], &d.matrix("GAMU").2[n + 1..]);

    // SPECAL
    let d = set.first("specal");
    let alfa = d.real("ALFA");
    let op = sol.specal(&pan, alfa, &settings(run));
    assert_golden!(ulps = ulps;
        "GAM" => sol.gam, d.reals("GAM");
        "GAM_A" => sol.gam_a, d.reals("GAM_A");
        "QINV" => sol.qinv[..n], d.reals("QINV");
        "QINV_A" => sol.qinv_a[..n], d.reals("QINV_A");
        "CPI" => sol.cpi, d.reals("CPI");
        "PSIO" => sol.psio, d.real("PSIO");
        "GAMTE" => sol.gamte, d.real("GAMTE");
        "SIGTE" => sol.sigte, d.real("SIGTE");
        "CL" => op.forces.cl, d.real("CL");
        "CM" => op.forces.cm, d.real("CM");
        "CDP" => op.forces.cdp, d.real("CDP");
        "CL_ALF" => op.forces.cl_alf, d.real("CL_ALF");
        "CL_MSQ" => op.forces.cl_msq, d.real("CL_MSQ");
        "MINF" => op.comp.minf, d.real("MINF");
        "TKLAM" => op.comp.tklam, d.real("TKLAM"));

    // XYWAKE
    let wake = Wake::trace(&pan, &sol.strengths(1.0), 1.0);
    let d = set.first("xywake");
    assert_eq!(wake.len() as i64, d.int("NW"), "NW");
    // APANELW starts at APANEL(N), the airfoil TE panel, followed by the wake panels
    let apanelw: Vec<f64> = std::iter::once(pan.apanel[n - 1])
        .chain(wake.apanel[..wake.len() - 1].iter().copied())
        .collect();
    assert_golden!(ulps = ulps;
        "XW" => wake.x, d.reals("XW");
        "YW" => wake.y, d.reals("YW");
        "SW" => wake.s, d.reals("SW");
        "NXW" => wake.nx, d.reals("NXW");
        "NYW" => wake.ny, d.reals("NYW");
        "APANELW" => apanelw, d.reals("APANELW"));

    // QWCALC + QDCALC
    sol.set_wake_speeds(&pan, &wake);
    sol.qdcalc(&pan, &wake);
    let d = set.first("qdcalc");
    let m = n + wake.len();
    let q = d.matrix("QINVU").2;
    assert_golden!(ulps = ulps;
        "QINVU0" => sol.qinvu[0], &q[..m];
        "QINVU1" => sol.qinvu[1], &q[m..];
        "DIJ" => sol.dij.as_ref().unwrap().as_slice(), d.matrix("DIJ").2);
}

#[test]
fn naca0012() {
    let run = "naca0012_re1e6_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
fn naca4412() {
    let run = "naca4412_re1e6_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
fn e387() {
    let run = "e387_re1e5_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
#[ignore = "needs full dumps: tools/gen_golden.sh"]
fn full_dump_runs() {
    for run in [
        "du91w2250_re1e6_a10.0",
        "naca0012_re1e6_a0.0",
        "naca0012_re1e6_a15.0",
        "naca4412_re1e6_a15.0",
    ] {
        check_run(run, &full_dump(run), 0);
    }
}
