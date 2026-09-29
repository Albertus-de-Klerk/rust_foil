//! State of the current operating point (XFOIL `/CR09/` quantities that change during a
//! solve). This is the `OperatingPoint` struct approved in PORTING_PLAN §4.

use crate::inviscid::forces::{Compressibility, Forces};

/// Angle of attack, integrated coefficients and freestream state of one operating point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OperatingPoint {
    /// Angle of attack (rad).
    pub alfa: f64,
    /// Freestream speed (always 1 in XFOIL).
    pub qinf: f64,
    /// CL, CM, CDP and their derivatives.
    pub forces: Forces,
    /// Mach/Re and Karman–Tsien terms at the current CL.
    pub comp: Compressibility,
    /// Total drag (viscous; zero for an inviscid point).
    pub cd: f64,
    /// Friction drag (viscous).
    pub cdf: f64,
}

impl OperatingPoint {
    /// An inviscid point from SPECAL.
    pub fn inviscid(alfa: f64, qinf: f64, forces: Forces, comp: Compressibility) -> Self {
        Self {
            alfa,
            qinf,
            forces,
            comp,
            cd: 0.0,
            cdf: 0.0,
        }
    }
}
