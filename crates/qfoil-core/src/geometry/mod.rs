//! Airfoil input and the buffer geometry (XFOIL `LOAD`, `NACA`, `AREAD`, `LEFIND`).

mod dat;
mod naca;

pub use naca::NacaDesignation;

use crate::error::{GeometryError, ParseError};
use crate::spline::{Curve, d2val, deval, seval};

/// Airfoil coordinates as supplied by the user (file order, not yet processed).
#[derive(Debug, Clone, PartialEq)]
pub struct Airfoil {
    /// Airfoil name (first line of a labelled file, or `NACA xxxx`).
    pub name: String,
    /// x coordinates.
    pub x: Vec<f64>,
    /// y coordinates.
    pub y: Vec<f64>,
    /// How the coordinates were produced; this selects QFoil's panelling path.
    pub source: AirfoilSource,
}

/// Origin of an [`Airfoil`], which decides how QFoil panels it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AirfoilSource {
    /// Read from a coordinate file (`LOAD`): nodes are used as panels (ABCOPY).
    File,
    /// Built-in NACA generator (`NACA`): re-panelled by PANGEN.
    Naca,
}

impl Airfoil {
    /// Parses a coordinate file (plain, labelled or single-element MSES). Port of XFOIL `AREAD`.
    pub fn from_dat(text: &str) -> Result<Self, ParseError> {
        dat::parse(text)
    }

    /// Generates a NACA 4- or 5-digit section. Port of XFOIL `NACA` / `NACA4` / `NACA5`.
    pub fn naca(designation: NacaDesignation) -> Result<Self, GeometryError> {
        naca::generate(designation)
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Whether there are no points.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}

/// QFoil buffer-array limit (`IBX = 4*IQX`).
pub const MAX_BUFFER_POINTS: usize = 4 * crate::limits::IQX;

/// The buffer airfoil: input points in counter-clockwise order, splined in arc length.
///
/// Corresponds to XFOIL's `XB, YB, SB, XBP, YBP` (COMMON `/CR14/`).
#[derive(Debug, Clone, PartialEq)]
pub struct Geometry {
    /// Airfoil name.
    pub name: String,
    /// Splined buffer points.
    pub buffer: Curve,
    /// Whether the input was clockwise and has been reversed (`LCLOCK`).
    pub reversed: bool,
    /// Panelling path implied by the input source.
    pub source: AirfoilSource,
}

impl Geometry {
    /// Builds the buffer airfoil. Port of the numerical part of XFOIL `LOAD` (orientation
    /// check, SCALC, SEGSPL). NACA input skips the orientation check, as `NACA` does.
    ///
    /// `GEOPAR` (thickness, camber, inertias) is not ported: it only produces printed
    /// information and does not affect the solution.
    pub fn new(airfoil: &Airfoil) -> Result<Self, GeometryError> {
        let nb = airfoil.len();
        if nb < 4 {
            return Err(GeometryError::TooFewPoints(nb));
        }
        if nb > MAX_BUFFER_POINTS {
            return Err(GeometryError::TooManyPoints {
                what: "buffer airfoil",
                got: nb,
                max: MAX_BUFFER_POINTS,
            });
        }
        let (mut x, mut y) = (airfoil.x.clone(), airfoil.y.clone());
        let mut reversed = false;
        if airfoil.source == AirfoilSource::File {
            // Signed area, positive for counter-clockwise ordering.
            let mut area = 0.0;
            for i in 0..nb {
                let ip = if i == nb - 1 { 0 } else { i + 1 };
                area += 0.5 * (y[i] + y[ip]) * (x[i] - x[ip]);
            }
            if area < 0.0 {
                x.reverse();
                y.reverse();
                reversed = true;
            }
        }
        Ok(Self {
            name: airfoil.name.clone(),
            buffer: Curve::new(x, y)?,
            reversed,
            source: airfoil.source,
        })
    }
}

/// Finds the leading-edge arc length: the point where the surface tangent is normal to
/// the line to the trailing-edge midpoint. Port of XFOIL `LEFIND`.
///
/// Returns `s[i]` of a sharp (doubled-point) LE directly. If the Newton iteration fails,
/// the first-guess node is returned, as in the Fortran.
pub fn lefind(x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64]) -> f64 {
    let n = x.len();
    let dseps = (s[n - 1] - s[0]) * 1.0e-5;
    let xte = 0.5 * (x[0] + x[n - 1]);
    let yte = 0.5 * (y[0] + y[n - 1]);

    // Fortran: DO I=3, N-2 ... GO TO 11 ; on fall-through the index is N-1.
    // 0-based: i in 2..=n-3, fall-through i = n-2.
    let mut i = n - 2;
    for k in 2..=n.saturating_sub(3) {
        let dxte = x[k] - xte;
        let dyte = y[k] - yte;
        let dx = x[k + 1] - x[k];
        let dy = y[k + 1] - y[k];
        if dxte * dx + dyte * dy < 0.0 {
            i = k;
            break;
        }
    }
    let mut sle = s[i];
    if s[i] == s[i - 1] {
        return sle;
    }
    for _ in 0..50 {
        let xle = seval(sle, x, xp, s);
        let yle = seval(sle, y, yp, s);
        let dxds = deval(sle, x, xp, s);
        let dyds = deval(sle, y, yp, s);
        let dxdd = d2val(sle, x, xp, s);
        let dydd = d2val(sle, y, yp, s);
        let xchord = xle - xte;
        let ychord = yle - yte;
        let res = xchord * dxds + ychord * dyds;
        let ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;
        let lim = 0.02 * (xchord + ychord).abs();
        let dsle = (-res / ress).max(-lim).min(lim);
        sle += dsle;
        if dsle.abs() < dseps {
            return sle;
        }
    }
    s[i]
}

impl Curve {
    /// Leading-edge arc length of this curve ([`lefind`]).
    pub fn leading_edge(&self) -> f64 {
        lefind(&self.x, &self.xp, &self.y, &self.yp, &self.s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clockwise_file_is_reversed() {
        let a = Airfoil {
            name: "tri".into(),
            x: vec![1.0, 0.5, 0.0, 0.5, 1.0],
            y: vec![0.0, -0.1, 0.0, 0.1, 0.001],
            source: AirfoilSource::File,
        };
        let g = Geometry::new(&a).unwrap();
        assert!(g.reversed);
        assert_eq!(g.buffer.y[1], 0.1);
    }
}
