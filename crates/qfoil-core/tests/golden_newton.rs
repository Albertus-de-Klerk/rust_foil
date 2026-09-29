//! Step 6 parity: BL set-up for a fresh point (STFIND, IBLPAN, XICALC, IBLSYS, UICALC,
//! QFoil D9 initialisation) and the first Newton iteration (SETBL, BLSOLV, UPDATE).

mod common;

use qfoil_core::bl::BoundaryLayer;
use qfoil_core::bl::march::MarchInputs;
use qfoil_core::bl::station::{Kernel, KernelParams};
use qfoil_core::bl::transition::AmplificationModel;
use qfoil_core::newton::{NewtonSystem, SetblInputs, UpdateInputs};
use qfoil_core::settings::{BlParams, FlowConditions};
use qfoil_core::{InviscidSolution, Paneling, Settings, Wake};
use qfoil_golden::{Dump, DumpSet, assert_golden, fixture, full_dump};

struct Case {
    pan: Paneling,
    wake: Wake,
    sol: InviscidSolution,
    flow: FlowConditions,
    cl: f64,
    alfa: f64,
}

/// Inviscid solution, wake and DIJ as VISCAL has them before its Newton loop.
fn prepare(run: &str, set: &DumpSet) -> Case {
    let (_, pan) = common::setup(run);
    let flow = FlowConditions {
        reynolds: common::reynolds(run),
        ..FlowConditions::default()
    };
    let settings = Settings {
        flow: flow.clone(),
        ..Settings::default()
    };
    let mut sol = InviscidSolution::new(&pan).unwrap();
    let alfa = set.first("specal").real("ALFA");
    let op = sol.specal(&pan, alfa, &settings);
    let wake = Wake::trace(&pan, &sol.strengths(1.0), 1.0);
    sol.set_wake_speeds(&pan, &wake);
    sol.qdcalc(&pan, &wake);
    Case {
        pan,
        wake,
        sol,
        flow,
        cl: op.forces.cl,
        alfa,
    }
}

fn side_arrays(bl: &BoundaryLayer, d: &Dump, names: &[&str], ulps: u64) {
    let mut failures = Vec::new();
    for (is, k) in [(0, "1"), (1, "2")] {
        let s = &bl.sides[is];
        for &name in names {
            let got: &[f64] = match name {
                "XSSI" => &s.xssi,
                "UINV" => &s.uinv,
                "UINV_A" => &s.uinv_a,
                "UEDG" => &s.uedg,
                "CTAU" => &s.ctau,
                "THET" => &s.thet,
                "DSTR" => &s.dstr,
                "MASS" => &s.mass,
                "TAU" => &s.tau,
                "DIS" => &s.dis,
                "CTQ" => &s.ctq,
                "DELT" => &s.delt,
                "TSTR" => &s.tstr,
                "VTI" => &s.vti,
                other => panic!("{other}"),
            };
            let want = d.reals(&format!("{name}{k}"));
            if let Err(m) = qfoil_golden::compare_slices(
                &format!("{name}{k}"),
                &got[..want.len()],
                want,
                ulps,
                0.0,
            ) {
                failures.push(m.to_string());
            }
        }
    }
    assert!(failures.is_empty(), "{}:\n{}", d.tag, failures.join("\n"));
}

fn flat(blocks: &[[[f64; 3]; 2]]) -> Vec<f64> {
    blocks
        .iter()
        .flat_map(|b| b.iter().flatten().copied())
        .collect()
}

fn check_run(run: &str, set: &DumpSet, ulps: u64) {
    let c = prepare(run, set);

    // ---- BL set-up (VISCAL before the Newton loop)
    let mut bl = BoundaryLayer::new(&c.pan, &c.wake, &c.sol);
    let d = set.first("viscal_init");
    assert_eq!(bl.ist as i64 + 1, d.int("IST"), "IST");
    assert_eq!(bl.nbl.map(|v| v as i64).to_vec(), d.ints("NBL"), "NBL");
    assert_eq!(
        bl.iblte.map(|v| v as i64 + 1).to_vec(),
        d.ints("IBLTE"),
        "IBLTE"
    );
    assert_eq!(bl.nsys as i64, d.int("NSYS"), "NSYS");
    for (is, k) in [(0, "1"), (1, "2")] {
        let n = bl.nbl[is];
        let ipan: Vec<i64> = bl.sides[is].ipan[1..n]
            .iter()
            .map(|&v| v as i64 + 1)
            .collect();
        assert_eq!(ipan, d.ints(&format!("IPAN{k}"))[1..], "IPAN{k}");
        let isys: Vec<i64> = bl.sides[is].isys[1..n]
            .iter()
            .map(|&v| v as i64 + 1)
            .collect();
        assert_eq!(isys, d.ints(&format!("ISYS{k}"))[1..], "ISYS{k}");
    }
    assert_golden!(ulps = ulps; "SST" => bl.sst, d.real("SST"); "WGAP" => bl.wgap, d.reals("WGAP"));
    side_arrays(
        &bl,
        d,
        &["XSSI", "UINV", "UINV_A", "UEDG", "CTAU", "VTI"],
        ulps,
    );

    // ---- iteration 1: SETBL (MRCHUE + MRCHDU + assembly)
    let dij = c.sol.dij.as_ref().unwrap();
    let mut k = Kernel::new(
        KernelParams::default(),
        BlParams::default(),
        AmplificationModel::Envelope,
    );
    let mut lblini = false;
    let inp = SetblInputs {
        march: MarchInputs {
            pan: &c.pan,
            acrit: [9.0, 9.0],
            xstrip: [1.0, 1.0],
        },
        dij,
        flow: &c.flow,
        clmr: c.cl,
        lalfa: true,
    };
    let mut sys: NewtonSystem = bl.setbl(&mut k, &mut lblini, &inp);
    let d = set.first("setbl");
    assert_golden!(ulps = ulps;
        "VA" => flat(&sys.va), d.matrix("VA").2;
        "VB" => flat(&sys.vb), d.matrix("VB").2;
        "VDEL" => flat(&sys.vdel), d.matrix("VDEL").2;
        "XSSITR" => bl.xssitr.to_vec(), d.reals("XSSITR"));
    side_arrays(
        &bl,
        d,
        &[
            "THET", "DSTR", "CTAU", "UEDG", "MASS", "TAU", "DIS", "CTQ", "DELT", "TSTR",
        ],
        ulps,
    );
    if let Some(vm) = set
        .dumps
        .values()
        .find(|x| x.tag == "setbl" && x.records.contains_key("VM"))
    {
        assert_golden!(ulps = ulps; "VM" => sys.vm, vm.matrix("VM").2);
    }

    // ---- BLSOLV
    let s = &c.pan.nodes.s;
    let ivte1 = bl.sides[0].isys[bl.iblte[0]];
    let ivz = bl.sides[1].isys[bl.iblte[1] + 1];
    sys.blsolv(0.01, s[s.len() - 1] - s[0], ivte1, ivz);
    let d = set.first("blsolv");
    assert_golden!(ulps = ulps; "VDEL" => flat(&sys.vdel), d.matrix("VDEL").2);

    // ---- UPDATE
    let (mut cl, mut alfa) = (c.cl, c.alfa);
    let upd = UpdateInputs {
        x: &c.pan.nodes.x,
        y: &c.pan.nodes.y,
        dij,
        minf: 0.0,
        minf_cl: 0.0,
        gamm1: c.flow.gamma - 1.0,
        lalfa: true,
        clspec: 0.0,
    };
    let r = bl.update(&sys, &upd, &mut cl, &mut alfa);
    let d = set.first("update");
    side_arrays(&bl, d, &["THET", "DSTR", "CTAU", "UEDG", "MASS"], ulps);
    assert_golden!(ulps = ulps;
        "CL" => cl, d.real("CL");
        "RLX" => r.rlx, d.real("RLX");
        "RMSBL" => r.rmsbl, d.real("RMSBL");
        "RMXBL" => r.rmxbl, d.real("RMXBL"));
}

#[test]
fn naca0012_iteration1() {
    let run = "naca0012_re1e6_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
fn naca4412_iteration1() {
    let run = "naca4412_re1e6_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
fn e387_iteration1() {
    let run = "e387_re1e5_a5.0";
    check_run(run, &fixture(run), 0);
}

#[test]
#[ignore = "needs full dumps: tools/gen_golden.sh"]
fn full_dump_iteration1() {
    for run in [
        "du91w2250_re1e6_a10.0",
        "naca0012_re1e6_a0.0",
        "naca0012_re1e6_a15.0",
        "naca4412_re1e6_a15.0",
        "naca0012_re1e6_a5.0_vm",
    ] {
        let data = run.trim_end_matches("_vm");
        let mut set = full_dump(data);
        if run.ends_with("_vm") {
            // the VM run only holds `setbl` records with VM; merge them in
            let vm = full_dump(run);
            let first_vm = vm.dumps.values().next().unwrap().clone();
            let key = set.dumps.values().find(|d| d.tag == "setbl").unwrap().seq;
            set.dumps
                .get_mut(&key)
                .unwrap()
                .records
                .extend(first_vm.records);
        }
        check_run(data, &set, 0);
    }
}
