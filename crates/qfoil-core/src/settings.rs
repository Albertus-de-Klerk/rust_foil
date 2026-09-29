//! User settings (XFOIL `INIT` / `BLPINI` defaults, OPER/VPAR parameters).

use crate::paneling::PanelingMode;

/// How the freestream Reynolds number varies with CL (`RETYP`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReynoldsType {
    /// Fixed Re (type 1).
    #[default]
    Fixed,
    /// Re ∝ 1/√CL, fixed lift (type 2).
    InverseSqrtCl,
    /// Re ∝ 1/CL, fixed lift and dynamic pressure (type 3).
    InverseCl,
}

/// How the freestream Mach number varies with CL (`MATYP`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MachType {
    /// Fixed Mach (type 1).
    #[default]
    Fixed,
    /// M ∝ 1/√CL, fixed lift (type 2).
    InverseSqrtCl,
    /// Fixed Mach, type 3 (identical to type 1 in MRCL).
    Fixed3,
}

/// Freestream conditions (`MINF1`, `REINF1`, `MATYP`, `RETYP`, `GAMMA`).
#[derive(Debug, Clone, PartialEq)]
pub struct FlowConditions {
    /// Mach number at CL = 1 (`MINF1`; the actual Mach for fixed-Mach polars).
    pub mach: f64,
    /// Reynolds number at CL = 1 (`REINF1`; the actual Re for fixed-Re polars).
    pub reynolds: f64,
    /// Mach dependence on CL.
    pub mach_type: MachType,
    /// Re dependence on CL.
    pub reynolds_type: ReynoldsType,
    /// Ratio of specific heats (`GAMMA`).
    pub gamma: f64,
}

impl Default for FlowConditions {
    fn default() -> Self {
        Self {
            mach: 0.0,
            reynolds: 0.0,
            mach_type: MachType::Fixed,
            reynolds_type: ReynoldsType::Fixed,
            gamma: 1.4,
        }
    }
}

/// All analysis settings.
#[derive(Debug, Clone, PartialEq)]
pub struct Settings {
    /// Freestream conditions.
    pub flow: FlowConditions,
    /// Panelling of the input airfoil.
    pub paneling: PanelingMode,
    /// Wake length / chord (`WAKLEN`).
    pub wake_length: f64,
    /// Moment reference point (`XCMREF`, `YCMREF`).
    pub moment_ref: (f64, f64),
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            flow: FlowConditions::default(),
            paneling: PanelingMode::Auto,
            wake_length: 1.0,
            moment_ref: (0.25, 0.0),
        }
    }
}
