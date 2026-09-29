//! Reader for `tests/golden/polars/<case>_re<Re>.tsv`.

use std::collections::BTreeMap;
use std::path::Path;

use crate::{GoldenError, golden_root};

/// One operating point of a golden polar (full precision, from the `final` dump).
#[derive(Debug, Clone, PartialEq)]
pub struct GoldenPoint {
    /// Angle of attack in degrees.
    pub alpha: f64,
    /// QFoil's `LVCONV` flag.
    pub converged: bool,
    /// Newton iterations used (the iteration limit if unconverged). `None` if the run died.
    pub niter: Option<u32>,
    /// Lift coefficient.
    pub cl: f64,
    /// Total drag coefficient.
    pub cd: f64,
    /// Pressure drag (surface-integrated, CLCALC's CDP).
    pub cdp: f64,
    /// Moment coefficient about (0.25, 0).
    pub cm: f64,
    /// Friction drag.
    pub cdf: f64,
    /// Transition x/c on the top side.
    pub xtr_top: f64,
    /// Transition x/c on the bottom side.
    pub xtr_bot: f64,
    /// Final RMS Newton change.
    pub rmsbl: f64,
}

/// A golden polar keyed by α in centidegrees, so iteration is in α order.
#[derive(Debug, Clone)]
pub struct GoldenPolar {
    /// Case id (`naca0012`, ...).
    pub case: String,
    /// Reynolds number.
    pub re: f64,
    /// Points by α × 100 (exact for the 0.5° grid).
    pub points: BTreeMap<i32, GoldenPoint>,
}

impl GoldenPolar {
    /// Reads `polars/<case>_re<re_label>.tsv`, e.g. `("naca0012", "1e6")`.
    pub fn read(case: &str, re_label: &str) -> Result<Self, GoldenError> {
        let path = golden_root()
            .join("polars")
            .join(format!("{case}_re{re_label}.tsv"));
        let text = std::fs::read_to_string(&path).map_err(|source| GoldenError::Io {
            path: path.clone(),
            source,
        })?;
        let re = re_label
            .parse()
            .map_err(|_| fmt_err(&path, 0, "bad Re label"))?;
        let mut points = BTreeMap::new();
        for (ln, line) in text.lines().enumerate().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 11 {
                return Err(fmt_err(&path, ln + 1, "expected 11 columns"));
            }
            let num = |k: usize| f[k].parse::<f64>().unwrap_or(f64::NAN);
            let alpha = num(0);
            let point = GoldenPoint {
                alpha,
                converged: f[1] == "1",
                niter: f[2].parse().ok(),
                cl: num(3),
                cd: num(4),
                cdp: num(5),
                cm: num(6),
                cdf: num(7),
                xtr_top: num(8),
                xtr_bot: num(9),
                rmsbl: num(10),
            };
            points.insert(alpha_key(alpha), point);
        }
        Ok(Self {
            case: case.to_owned(),
            re,
            points,
        })
    }

    /// Point at `alpha` degrees.
    pub fn at(&self, alpha: f64) -> Option<&GoldenPoint> {
        self.points.get(&alpha_key(alpha))
    }
}

/// α in centidegrees, the key of [`GoldenPolar::points`].
pub fn alpha_key(alpha: f64) -> i32 {
    (alpha * 100.0).round() as i32
}

fn fmt_err(path: &Path, line: usize, message: &str) -> GoldenError {
    GoldenError::Format {
        path: path.to_owned(),
        line,
        message: message.to_owned(),
    }
}
