//! Phase 3 parity: every point of the 915-point golden matrix, via the public API.
//!
//! Run with `cargo test --release -p qfoil-core --test golden_polars -- --ignored`.

mod common;

use std::collections::BTreeMap;

use qfoil_core::polar::PreparedAirfoil;
use qfoil_core::settings::FlowConditions;
use qfoil_core::{PointStatus, Settings, ViscousSettings};
use qfoil_golden::{GoldenPolar, ulp_distance};

/// Per-case outcome: (points, bit-identical, flag mismatches, worst |ΔCL|, worst rel ΔCD).
fn run_case(case: &str, re: &str) -> (usize, usize, Vec<f64>, f64, f64) {
    let golden = GoldenPolar::read(case, re).unwrap();
    let (_, _) = common::setup(&format!("{case}_re{re}"));
    let airfoil = match case.strip_prefix("naca") {
        Some(c) => {
            qfoil_core::Airfoil::naca(qfoil_core::NacaDesignation::new(c.parse().unwrap()).unwrap())
                .unwrap()
        }
        None => qfoil_core::Airfoil::from_dat(
            &std::fs::read_to_string(qfoil_golden::airfoil_path(case)).unwrap(),
        )
        .unwrap(),
    };
    let settings = Settings {
        flow: FlowConditions {
            reynolds: golden.re,
            ..FlowConditions::default()
        },
        ..Settings::default()
    };
    let vs = ViscousSettings::default();
    let prep = PreparedAirfoil::new(&airfoil, &settings).unwrap();

    let results: BTreeMap<i32, _> = std::thread::scope(|scope| {
        let handles: Vec<_> = golden
            .points
            .iter()
            .map(|(&key, g)| {
                let (prep, settings, vs) = (&prep, &settings, &vs);
                scope.spawn(move || (key, prep.point(g.alpha, settings, Some(vs))))
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    let (mut exact, mut flag_mismatch, mut dcl, mut dcd) = (0, Vec::new(), 0.0_f64, 0.0_f64);
    for (key, g) in &golden.points {
        let p = &results[key];
        let conv = matches!(p.status, PointStatus::Converged { .. });
        if conv != g.converged {
            flag_mismatch.push(g.alpha);
        }
        let same = [
            (p.cl, g.cl),
            (p.cd, g.cd),
            (p.cdp, g.cdp),
            (p.cm, g.cm),
            (p.cdf, g.cdf),
            (p.xtr[0], g.xtr_top),
            (p.xtr[1], g.xtr_bot),
        ]
        .iter()
        .all(|&(a, b)| ulp_distance(a, b) == 0 || (a.is_nan() && b.is_nan()));
        if same {
            exact += 1;
        } else if std::env::var_os("GOLDEN_LIST").is_some() {
            println!(
                "DIFF {case} {re} {:+.1} conv {conv}/{} dCL {:.2e}",
                g.alpha,
                g.converged,
                p.cl - g.cl
            );
        }
        if conv && g.converged {
            dcl = dcl.max((p.cl - g.cl).abs());
            dcd = dcd.max(((p.cd - g.cd) / g.cd).abs());
        }
    }
    (golden.points.len(), exact, flag_mismatch, dcl, dcd)
}

#[test]
#[ignore = "915 viscous points: run with --release -- --ignored"]
fn all_golden_polars() {
    let mut report = Vec::new();
    let mut total = (0, 0, 0);
    for case in ["naca0012", "naca0020", "naca4412", "du91w2250", "e387"] {
        for re in ["1e5", "1e6", "5e6"] {
            let (n, exact, flags, dcl, dcd) = run_case(case, re);
            total.0 += n;
            total.1 += exact;
            total.2 += flags.len();
            report.push(format!("{case:10} Re {re}: {exact}/{n} bit-identical, flag mismatches {flags:?}, max|dCL| {dcl:.1e}, max rel dCD {dcd:.1e}"));
        }
    }
    println!("{}", report.join("\n"));
    println!(
        "TOTAL {}/{} bit-identical, {} flag mismatches",
        total.1, total.0, total.2
    );
    assert_eq!(
        total.1,
        total.0,
        "not all points bit-identical:\n{}",
        report.join("\n")
    );
}
