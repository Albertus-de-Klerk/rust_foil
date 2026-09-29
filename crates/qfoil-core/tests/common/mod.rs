//! Shared setup for golden parity tests: rebuild each reference run's airfoil.

#![allow(dead_code)]

use qfoil_core::{Airfoil, Geometry, NacaDesignation, Paneling, PanelingMode};
use qfoil_golden::airfoil_path;

/// Geometry and panels for a golden run named `<case>_re<Re>_a<alpha>`.
pub fn setup(run: &str) -> (Geometry, Paneling) {
    let case = run.split('_').next().unwrap();
    let airfoil = match case.strip_prefix("naca") {
        Some(code) => Airfoil::naca(NacaDesignation::new(code.parse().unwrap()).unwrap()).unwrap(),
        None => Airfoil::from_dat(&std::fs::read_to_string(airfoil_path(case)).unwrap()).unwrap(),
    };
    let geom = Geometry::new(&airfoil).unwrap();
    let pan = Paneling::new(&geom, &PanelingMode::Auto).unwrap();
    (geom, pan)
}

/// Reynolds number encoded in a run name.
pub fn reynolds(run: &str) -> f64 {
    run.split('_')
        .find_map(|p| p.strip_prefix("re"))
        .unwrap()
        .parse()
        .unwrap()
}
