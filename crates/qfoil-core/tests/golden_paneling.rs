//! Steps 1–2 parity: buffer geometry and panel nodes against the reference dumps
//! (`pangen` for NACA sections, `abcopy` for loaded files).

use qfoil_core::{Airfoil, Geometry, NacaDesignation, Paneling, PanelingMode};
use qfoil_golden::{Dump, DumpSet, airfoil_path, assert_golden, fixture, full_dump};

/// Bit-level parity is required for geometry.
const ULPS: u64 = 0;

fn check(dump: &Dump, geom: &Geometry, pan: &Paneling) {
    let n = pan.len();
    assert_eq!(dump.int("N") as usize, n, "panel node count");
    assert_eq!(
        dump.int("NB") as usize,
        geom.buffer.len(),
        "buffer point count"
    );
    assert_eq!(dump.int("SHARP") == 1, pan.trailing_edge.sharp, "SHARP");
    let b = &geom.buffer;
    let p = &pan.nodes;
    assert_golden!(ulps = ULPS;
        "XB" => b.x, dump.reals("XB");
        "YB" => b.y, dump.reals("YB");
        "SB" => b.s, dump.reals("SB");
        "X" => p.x, dump.reals("X");
        "Y" => p.y, dump.reals("Y");
        "S" => p.s, dump.reals("S");
        "XP" => p.xp, dump.reals("XP");
        "YP" => p.yp, dump.reals("YP");
        "NX" => pan.nx, dump.reals("NX");
        "NY" => pan.ny, dump.reals("NY");
        "APANEL" => pan.apanel, dump.reals("APANEL");
        "SLE" => pan.sle, dump.real("SLE");
        "XLE" => pan.le.0, dump.real("XLE");
        "YLE" => pan.le.1, dump.real("YLE");
        "XTE" => pan.te.0, dump.real("XTE");
        "YTE" => pan.te.1, dump.real("YTE");
        "CHORD" => pan.chord, dump.real("CHORD");
        "ANTE" => pan.trailing_edge.ante, dump.real("ANTE");
        "ASTE" => pan.trailing_edge.aste, dump.real("ASTE");
        "DSTE" => pan.trailing_edge.dste, dump.real("DSTE");
    );
    if dump.records.contains_key("XBP") {
        assert_golden!(ulps = ULPS;
            "XBP" => b.xp, dump.reals("XBP");
            "YBP" => b.yp, dump.reals("YBP"));
    }
}

fn naca_case(set: &DumpSet, code: u32) {
    let airfoil = Airfoil::naca(NacaDesignation::new(code).unwrap()).unwrap();
    let geom = Geometry::new(&airfoil).unwrap();
    let pan = Paneling::new(&geom, &PanelingMode::Auto).unwrap();
    check(set.first("pangen"), &geom, &pan);
}

fn file_case(set: &DumpSet, name: &str) {
    let text = std::fs::read_to_string(airfoil_path(name)).unwrap();
    let airfoil = Airfoil::from_dat(&text).unwrap();
    let geom = Geometry::new(&airfoil).unwrap();
    let pan = Paneling::new(&geom, &PanelingMode::Auto).unwrap();
    check(set.first("abcopy"), &geom, &pan);
}

#[test]
fn naca0012_pangen() {
    naca_case(&fixture("naca0012_re1e6_a5.0"), 12);
}

#[test]
fn naca4412_pangen() {
    naca_case(&fixture("naca4412_re1e6_a5.0"), 4412);
}

#[test]
fn e387_abcopy() {
    file_case(&fixture("e387_re1e5_a5.0"), "e387");
}

#[test]
#[ignore = "needs full dumps: tools/gen_golden.sh"]
fn du91w2250_abcopy() {
    file_case(&full_dump("du91w2250_re1e6_a10.0"), "du91w2250");
}
