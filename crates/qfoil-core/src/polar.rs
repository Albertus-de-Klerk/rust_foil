//! Polar sweeps: the public entry point [`analyse_polar`].
//!
//! Every angle of attack is solved independently from a cold start, which is how QBlade
//! runs QFoil (one process per α). QFoil 0.9 cannot sweep within one session: see
//! PORTING_PLAN S10. The α-independent work (panelling, the vortex influence
//! factorisation GGCALC) is done once per airfoil.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt;

use crate::error::GeometryError;
use crate::geometry::{Airfoil, Geometry};
use crate::inviscid::InviscidSolution;
use crate::linalg::LinalgError;
use crate::paneling::Paneling;
use crate::settings::Settings;
use crate::viscous::{ViscousSettings, viscal};

/// Errors that prevent a polar from being computed at all.
///
/// Non-convergence of individual points is not an error; see [`PointStatus`].
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum AnalysisError {
    /// Invalid airfoil or panelling.
    #[error(transparent)]
    Geometry(#[from] GeometryError),
    /// The panel system exceeds QFoil's size limits.
    #[error(transparent)]
    Linalg(#[from] LinalgError),
    /// Invalid settings.
    #[error("invalid settings: {0}")]
    Settings(String),
}

/// An angle of attack in degrees, totally ordered (for use as a [`BTreeMap`] key).
#[derive(Debug, Clone, Copy)]
pub struct Alpha(pub f64);

impl PartialEq for Alpha {
    fn eq(&self, other: &Self) -> bool {
        self.0.total_cmp(&other.0) == Ordering::Equal
    }
}
impl Eq for Alpha {}
impl PartialOrd for Alpha {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Alpha {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}
impl fmt::Display for Alpha {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}°", self.0)
    }
}

/// Angles of attack to analyse.
#[derive(Debug, Clone, PartialEq)]
pub enum AlphaSchedule {
    /// `start, start+step, ...` up to `end` (QFoil ASEQ point count:
    /// `INT((end-start)/step + 0.5) + 1`).
    Range {
        /// First α (deg).
        start: f64,
        /// Last α (deg).
        end: f64,
        /// Increment (deg); its sign is taken from `end - start`.
        step: f64,
    },
    /// Explicit list (deg).
    List(Vec<f64>),
}

impl AlphaSchedule {
    /// The α values in degrees.
    pub fn values(&self) -> Vec<f64> {
        match self {
            Self::List(v) => v.clone(),
            &Self::Range { start, end, step } => {
                let step = if end < start { -step.abs() } else { step.abs() };
                let n = if step == 0.0 {
                    1
                } else {
                    ((end - start) / step + 0.5) as usize + 1
                };
                (0..n).map(|i| start + step * i as f64).collect()
            }
        }
    }
}

/// Everything [`analyse_polar`] needs besides the airfoil.
#[derive(Debug, Clone, PartialEq)]
pub struct PolarSettings {
    /// Angles of attack.
    pub alpha: AlphaSchedule,
    /// Freestream, panelling and geometric settings.
    pub settings: Settings,
    /// Viscous parameters; `None` for an inviscid polar.
    pub viscous: Option<ViscousSettings>,
}

/// Outcome of one operating point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PointStatus {
    /// Inviscid point (always "converged").
    Inviscid,
    /// The Newton iteration reached `RMSBL < 1e-4`.
    Converged {
        /// Iterations used.
        iterations: usize,
    },
    /// The iteration limit was reached. QFoil drops such points from its polar file.
    NotConverged {
        /// Iterations used.
        iterations: usize,
        /// Final RMS Newton change.
        rms: f64,
    },
}

impl PointStatus {
    /// Whether the point belongs in a polar (inviscid or converged).
    pub fn is_converged(&self) -> bool {
        !matches!(self, Self::NotConverged { .. })
    }
}

/// Results of one operating point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolarPoint {
    /// Angle of attack (deg).
    pub alpha: f64,
    /// Lift coefficient.
    pub cl: f64,
    /// Total drag coefficient (wake momentum deficit, QFoil GWAKE correction).
    pub cd: f64,
    /// Pressure drag from surface pressure integration (the polar file's `CDp`).
    pub cdp: f64,
    /// Friction drag.
    pub cdf: f64,
    /// Moment coefficient about the reference point.
    pub cm: f64,
    /// Transition x/c on the top and bottom surface (1 if laminar to the TE).
    pub xtr: [f64; 2],
    /// Fractional panel index of transition, top and bottom (the polar file's `Itr`).
    pub itr: [f64; 2],
    /// Convergence.
    pub status: PointStatus,
}

/// A polar: operating points in α order.
#[derive(Debug, Clone, PartialEq)]
pub struct Polar {
    /// Airfoil name.
    pub name: String,
    /// Settings used.
    pub settings: PolarSettings,
    /// Points keyed by α.
    pub points: BTreeMap<Alpha, PolarPoint>,
}

impl Polar {
    /// Converged (or inviscid) points in α order.
    pub fn converged(&self) -> impl Iterator<Item = &PolarPoint> {
        self.points.values().filter(|p| p.status.is_converged())
    }
}

/// π/180 as XFOIL computes it (`DTOR = PI/180.0`).
const DTOR: f64 = crate::fortran::PI / 180.0;

/// Airfoil preparation shared by all points of a polar.
#[derive(Debug, Clone)]
pub struct PreparedAirfoil {
    /// Panels.
    pub paneling: Paneling,
    /// Inviscid solution with the factored influence matrix (α-independent part).
    pub base: InviscidSolution,
}

impl PreparedAirfoil {
    /// Buffer geometry, panels and GGCALC.
    pub fn new(airfoil: &Airfoil, settings: &Settings) -> Result<Self, AnalysisError> {
        let geom = Geometry::new(airfoil)?;
        let paneling = Paneling::new(&geom, &settings.paneling)?;
        let base = InviscidSolution::new(&paneling)?;
        Ok(Self { paneling, base })
    }

    /// Analyses one angle of attack (deg) from a cold start, as QFoil's `ALFA` does in a
    /// fresh session.
    pub fn point(
        &self,
        alpha_deg: f64,
        settings: &Settings,
        viscous: Option<&ViscousSettings>,
    ) -> PolarPoint {
        let pan = &self.paneling;
        let mut sol = self.base.clone();
        let alfa = DTOR * alpha_deg;
        let mut op = sol.specal(pan, alfa, settings);
        let Some(vs) = viscous else {
            return PolarPoint {
                alpha: alpha_deg,
                cl: op.forces.cl,
                cd: 0.0,
                cdp: op.forces.cdp,
                cdf: 0.0,
                cm: op.forces.cm,
                xtr: [1.0, 1.0],
                itr: [0.0, 0.0],
                status: PointStatus::Inviscid,
            };
        };
        let res = viscal(pan, &mut sol, &mut op, settings, vs, |_| {});
        let status = if res.converged {
            PointStatus::Converged {
                iterations: res.iterations,
            }
        } else {
            PointStatus::NotConverged {
                iterations: res.iterations,
                rms: res.last.map_or(f64::NAN, |r| r.rmsbl),
            }
        };
        PolarPoint {
            alpha: alpha_deg,
            cl: op.forces.cl,
            cd: op.cd,
            cdp: op.forces.cdp,
            cdf: op.cdf,
            cm: op.forces.cm,
            xtr: res.bl.xoctr,
            itr: res.bl.tindex,
            status,
        }
    }
}

/// Computes a polar for `airfoil`.
///
/// Each α is an independent cold-start solve (QBlade semantics), so results do not
/// depend on the order or grouping of the angles.
pub fn analyse_polar(airfoil: &Airfoil, settings: &PolarSettings) -> Result<Polar, AnalysisError> {
    if let Some(vs) = &settings.viscous
        && settings.settings.flow.reynolds <= 0.0
    {
        return Err(AnalysisError::Settings(format!(
            "viscous analysis needs a positive Reynolds number (got {}, Ncrit {:?})",
            settings.settings.flow.reynolds, vs.ncrit
        )));
    }
    let prepared = PreparedAirfoil::new(airfoil, &settings.settings)?;
    let points = settings
        .alpha
        .values()
        .into_iter()
        .map(|a| {
            (
                Alpha(a),
                prepared.point(a, &settings.settings, settings.viscous.as_ref()),
            )
        })
        .collect();
    Ok(Polar {
        name: airfoil.name.clone(),
        settings: settings.clone(),
        points,
    })
}
