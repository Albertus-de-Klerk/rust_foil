//! Steps 7–8 parity: the complete VISCAL Newton loop (SETBL, BLSOLV, UPDATE, QVFUE, GAMQV,
//! STMOVE, CLCALC, CDCALC incl. QFoil GWAKE) at every dumped iteration, and the final point.

mod common;

use std::collections::BTreeMap;

use qfoil_core::settings::FlowConditions;
use qfoil_core::viscous::{Stage, ViscousSettings, viscal};
use qfoil_core::{InviscidSolution, Settings};
use qfoil_golden::{Dump, DumpSet, compare_slices, fixture, full_dump};

/// Dumps keyed by (stage, iteration). SETBL records carry ITER = 0, so their iteration is
/// their position; the trimmed fixtures keep SETBL 1, 2 and last (= final iteration).
fn index(set: &DumpSet) -> BTreeMap<(&'static str, usize), &Dump> {
    let mut m = BTreeMap::new();
    let last_iter = set.last("iter").int("ITER") as usize;
    let setbl: Vec<&Dump> = set.all("setbl").collect();
    let full = setbl.len() == last_iter;
    for (pos, d) in setbl.iter().enumerate() {
        let it = if full || pos < 2 { pos + 1 } else { last_iter };
        m.insert(("setbl", it), *d);
    }
    for tag in ["blsolv", "update", "iter"] {
        for d in set.all(tag) {
            let t: &'static str = match tag {
                "blsolv" => "blsolv",
                "update" => "update",
                _ => "iter",
            };
            m.insert((t, d.int("ITER") as usize), d);
        }
    }
    m
}

fn bl_mismatches(
    bl: &qfoil_core::bl::BoundaryLayer,
    d: &Dump,
    names: &[&str],
    ulps: u64,
) -> Vec<String> {
    let mut out = Vec::new();
    for (is, k) in [(0, "1"), (1, "2")] {
        let s = &bl.sides[is];
        for &name in names {
            let got: &[f64] = match name {
                "THET" => &s.thet,
                "DSTR" => &s.dstr,
                "CTAU" => &s.ctau,
                "UEDG" => &s.uedg,
                "MASS" => &s.mass,
                "XSSI" => &s.xssi,
                "TAU" => &s.tau,
                _ => unreachable!(),
            };
            let want = d.reals(&format!("{name}{k}"));
            if let Err(m) = compare_slices(
                &format!("{} {name}{k}", d.tag),
                &got[..want.len()],
                want,
                ulps,
                0.0,
            ) {
                out.push(m.to_string());
            }
        }
    }
    out
}

fn check_run(run: &str, set: &DumpSet, ulps: u64) {
    let (_, pan) = common::setup(run);
    let settings = Settings {
        flow: FlowConditions {
            reynolds: common::reynolds(run),
            ..FlowConditions::default()
        },
        ..Settings::default()
    };
    let vs = ViscousSettings::default();
    let mut sol = InviscidSolution::new(&pan).unwrap();
    let alfa = set.first("specal").real("ALFA");
    let mut op = sol.specal(&pan, alfa, &settings);

    let dumps = index(set);
    let mut failures: Vec<String> = Vec::new();
    let mut checked = 0;
    let res = viscal(&pan, &mut sol, &mut op, &settings, &vs, |v| {
        let key = match v.stage {
            Stage::Assembled => "setbl",
            Stage::Solved => "blsolv",
            Stage::Updated => "update",
            Stage::Iterated => "iter",
        };
        let Some(d) = dumps.get(&(key, v.iteration)) else {
            return;
        };
        checked += 1;
        let names: &[&str] = match v.stage {
            Stage::Assembled | Stage::Updated => &["THET", "DSTR", "CTAU", "UEDG", "MASS"],
            Stage::Solved => &[],
            Stage::Iterated => &["THET", "DSTR", "CTAU", "UEDG", "MASS", "XSSI", "TAU"],
        };
        failures.extend(bl_mismatches(v.bl, d, names, ulps));
        let scalars: Vec<(&str, f64)> = match v.stage {
            Stage::Updated => vec![
                ("RLX", v.update.unwrap().rlx),
                ("RMSBL", v.update.unwrap().rmsbl),
                ("CL", v.op.forces.cl),
            ],
            Stage::Iterated => vec![
                ("CL", v.op.forces.cl),
                ("CM", v.op.forces.cm),
                ("CD", v.op.cd),
                ("CDF", v.op.cdf),
                ("CDP", v.op.forces.cdp),
                ("SST", v.bl.sst),
            ],
            _ => vec![],
        };
        for (name, got) in scalars {
            if let Err(m) = compare_slices(
                &format!("{key} it{} {name}", v.iteration),
                &[got],
                &[d.real(name)],
                ulps,
                0.0,
            ) {
                failures.push(m.to_string());
            }
        }
        if v.stage == Stage::Iterated {
            let ist = v.bl.ist as i64 + 1;
            if ist != d.int("IST") {
                failures.push(format!(
                    "iter {} IST {ist} vs {}",
                    v.iteration,
                    d.int("IST")
                ));
            }
        }
        if v.stage == Stage::Solved {
            let flat: Vec<f64> = v
                .sys
                .vdel
                .iter()
                .flat_map(|b| b.iter().flatten().copied())
                .collect();
            if let Err(m) = compare_slices(
                &format!("blsolv it{} VDEL", v.iteration),
                &flat,
                d.matrix("VDEL").2,
                ulps,
                0.0,
            ) {
                failures.push(m.to_string());
            }
        }
    });
    assert!(failures.is_empty(), "{run}:\n{}", failures.join("\n"));
    assert!(
        checked >= 4,
        "{run}: only {checked} dumps matched an iteration"
    );

    // final point
    let f = set.last("final");
    assert_eq!(
        res.converged,
        f.int("CONVERGED") == 1,
        "{run}: converged flag"
    );
    let niter = if res.converged {
        f.int("ITER")
    } else {
        vs.max_iterations as i64
    };
    assert_eq!(res.iterations as i64, niter, "{run}: iterations");
    let e = set.last("viscal_end");
    let mut fails = Vec::new();
    for (name, got, want) in [
        ("CL", op.forces.cl, f.real("CL")),
        ("CM", op.forces.cm, f.real("CM")),
        ("CD", op.cd, f.real("CD")),
        ("CDF", op.cdf, f.real("CDF")),
        ("CDP", op.forces.cdp, f.real("CDP")),
        ("XOCTR1", res.bl.xoctr[0], f.reals("XOCTR")[0]),
        ("XOCTR2", res.bl.xoctr[1], f.reals("XOCTR")[1]),
    ] {
        if let Err(m) = compare_slices(name, &[got], &[want], ulps, 0.0) {
            fails.push(m.to_string());
        }
    }
    for (name, got, want) in [
        ("QVIS", &sol.qvis, e.reals("QVIS")),
        ("CPV", &res.cpv, e.reals("CPV")),
        ("CPI", &sol.cpi, e.reals("CPI")),
    ] {
        if let Err(m) = compare_slices(name, got, want, ulps, 0.0) {
            fails.push(m.to_string());
        }
    }
    assert!(fails.is_empty(), "{run} final:\n{}", fails.join("\n"));
}

#[test]
fn naca0012_viscal() {
    let run = "naca0012_re1e6_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
fn naca4412_viscal() {
    let run = "naca4412_re1e6_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
fn e387_viscal() {
    let run = "e387_re1e5_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
#[ignore = "needs full dumps: tools/gen_golden.sh"]
fn full_dump_viscal() {
    for run in [
        "du91w2250_re1e6_a10.0",
        "naca0012_re1e6_a0.0",
        "naca0012_re1e6_a15.0",
        "naca4412_re1e6_a15.0",
    ] {
        check_run(run, &full_dump(run), 0);
    }
}

/// Debug helper: `GOLDEN_RUN=<dump dir name> cargo test --release --test golden_viscal debug_run -- --ignored`
#[test]
#[ignore = "debug helper"]
fn debug_run() {
    let run = std::env::var("GOLDEN_RUN").expect("set GOLDEN_RUN");
    check_run(&run, &full_dump(&run), 0);
}
