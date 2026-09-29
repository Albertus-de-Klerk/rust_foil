//! NACA 4- and 5-digit sections. Port of XFOIL `NACA`, `NACA4`, `NACA5` (naca.f).

use std::fmt;
use std::str::FromStr;

use super::{Airfoil, AirfoilSource};
use crate::error::GeometryError;
use crate::fortran::powi;
use crate::limits::IQX;

/// Points per side: `NSIDE = IQX/3` (xfoil.f `NACA`).
///
/// This depends on QFoil's raised `IQX = 1400` (PORTING_PLAN D13): 466 points per side and
/// 931 buffer points, compared with 123 per side in stock XFOIL 6.99.
pub const NSIDE: usize = IQX / 3;

/// TE point bunching exponent (`AN`).
const AN: f64 = 1.5;

/// A NACA designation, validated like XFOIL's `NACA` command.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NacaDesignation(u32);

impl NacaDesignation {
    /// Creates a designation. Four digits (`0..=9999`) or five digits with mean-line code
    /// 210–250 (`21000..=25099`).
    pub fn new(code: u32) -> Result<Self, GeometryError> {
        let ok = code <= 9999
            || ((21000..=25099).contains(&code)
                && (210..=250).contains(&(code / 100))
                && (code / 100) % 10 == 0);
        if ok {
            Ok(Self(code))
        } else {
            Err(GeometryError::NacaDesignation(code))
        }
    }

    /// The numeric designation.
    pub fn code(self) -> u32 {
        self.0
    }
}

impl FromStr for NacaDesignation {
    type Err = GeometryError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let digits = s
            .trim()
            .trim_start_matches(['N', 'n', 'A', 'a', 'C', 'c', ' ']);
        digits
            .parse()
            .map_err(|_| GeometryError::NacaDesignation(0))
            .and_then(Self::new)
    }
}

impl fmt::Display for NacaDesignation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0 <= 9999 {
            write!(f, "NACA {:04}", self.0)
        } else {
            write!(f, "NACA {:05}", self.0)
        }
    }
}

/// Thickness distribution shared by NACA4 and NACA5.
fn thickness(xx: f64, t: f64) -> f64 {
    (0.29690 * xx.sqrt() - 0.12600 * xx - 0.35160 * powi(xx, 2) + 0.28430 * powi(xx, 3)
        - 0.10150 * powi(xx, 4))
        * t
        / 0.20
}

/// Chordwise stations with TE bunching.
fn stations() -> Vec<f64> {
    let anp = AN + 1.0;
    (0..NSIDE)
        .map(|i| {
            if i == NSIDE - 1 {
                1.0
            } else {
                let frac = i as f64 / (NSIDE - 1) as f64;
                1.0 - anp * frac * (1.0 - frac).powf(AN) - (1.0 - frac).powf(anp)
            }
        })
        .collect()
}

/// Assembles the buffer: upper surface from TE to LE, then lower surface LE to TE.
fn assemble(name: String, xx: &[f64], yt: &[f64], yc: &[f64]) -> Airfoil {
    let upper = (0..NSIDE).rev().map(|i| (xx[i], yc[i] + yt[i]));
    let lower = (1..NSIDE).map(|i| (xx[i], yc[i] - yt[i]));
    let (x, y) = upper.chain(lower).unzip();
    Airfoil {
        name,
        x,
        y,
        source: AirfoilSource::Naca,
    }
}

pub(super) fn generate(des: NacaDesignation) -> Result<Airfoil, GeometryError> {
    let code = des.code();
    let xx = stations();
    let digits = |k: u32| (code / 10u32.pow(k)) % 10;
    let t = f64::from(digits(1) * 10 + digits(0)) / 100.0;
    let yt: Vec<f64> = xx.iter().map(|&x| thickness(x, t)).collect();

    let yc: Vec<f64> = if code <= 9999 {
        // NACA4: M = N4/100, P = N3/10
        let m = f64::from(digits(3)) / 100.0;
        let p = f64::from(digits(2)) / 10.0;
        xx.iter()
            .map(|&x| {
                if x < p {
                    m / powi(p, 2) * (2.0 * p * x - powi(x, 2))
                } else {
                    m / powi(1.0 - p, 2) * ((1.0 - 2.0 * p) + 2.0 * p * x - powi(x, 2))
                }
            })
            .collect()
    } else {
        // NACA5: tabulated mean-line constants for the first three digits
        let (m, c) = match code / 100 {
            210 => (0.0580, 361.4),
            220 => (0.1260, 51.64),
            230 => (0.2025, 15.957),
            240 => (0.2900, 6.643),
            250 => (0.3910, 3.230),
            _ => return Err(GeometryError::NacaDesignation(code)),
        };
        xx.iter()
            .map(|&x| {
                if x < m {
                    (c / 6.0) * (powi(x, 3) - 3.0 * m * powi(x, 2) + m * m * (3.0 - m) * x)
                } else {
                    (c / 6.0) * powi(m, 3) * (1.0 - x)
                }
            })
            .collect()
    };
    Ok(assemble(des.to_string(), &xx, &yt, &yc))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designation_validation() {
        assert!(NacaDesignation::new(12).is_ok());
        assert!(NacaDesignation::new(23012).is_ok());
        assert!(NacaDesignation::new(26012).is_err());
        assert!(NacaDesignation::new(21512).is_err());
        assert_eq!("NACA 4412".parse::<NacaDesignation>().unwrap().code(), 4412);
        assert_eq!(NacaDesignation::new(12).unwrap().to_string(), "NACA 0012");
    }

    #[test]
    fn buffer_shape() {
        let a = generate(NacaDesignation::new(12).unwrap()).unwrap();
        assert_eq!(a.len(), 2 * NSIDE - 1);
        assert_eq!((a.x[0], a.x[NSIDE - 1]), (1.0, 0.0));
        // Closed TE thickness of the NACA formula: 2 * 0.00126 * 0.12/0.2 * ... ~ 0.00252
        assert!((a.y[0] - a.y[a.len() - 1] - 0.00252).abs() < 1e-5);
    }
}
