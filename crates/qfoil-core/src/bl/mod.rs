//! Integral boundary-layer model: closures, transition, station equations and marching.
//!
//! Index convention: BL station arrays are 0-based with Rust index `ibl` = Fortran
//! `IBL - 1`; index 0 is XFOIL's dummy stagnation station `IBL = 1`. Station indices stored
//! in [`BoundaryLayer`] (`iblte`, `itran`) use the same shift, so Fortran comparisons such
//! as `IBL.LT.ITRAN(IS)` translate unchanged. Side 0 is the upper surface (Fortran
//! `IS = 1`), side 1 the lower surface plus wake.

pub mod closure;
pub mod coupling;
pub mod equations;
pub mod march;
pub mod station;
pub mod transition;

/// BL arrays of one side (COMMON `/CR15/`, `/CI05/` per-side parts).
#[allow(missing_docs)]
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Side {
    /// BL arc length ξ from the stagnation point.
    pub xssi: Vec<f64>,
    /// Edge speed (incompressible).
    pub uedg: Vec<f64>,
    /// Inviscid edge speed.
    pub uinv: Vec<f64>,
    /// dUINV/dα.
    pub uinv_a: Vec<f64>,
    /// Mass defect δ*·Ue.
    pub mass: Vec<f64>,
    /// Momentum thickness θ.
    pub thet: Vec<f64>,
    /// Displacement thickness δ*.
    pub dstr: Vec<f64>,
    /// Amplification exponent N (laminar stations) or √Cτ (turbulent stations).
    pub ctau: Vec<f64>,
    /// BL thickness δ.
    pub delt: Vec<f64>,
    /// Kinetic-energy thickness θ*.
    pub tstr: Vec<f64>,
    /// Wall shear stress.
    pub tau: Vec<f64>,
    /// Dissipation.
    pub dis: Vec<f64>,
    /// Equilibrium √Cτ.
    pub ctq: Vec<f64>,
    /// Slip-velocity factor 1.6/(1+Us), output only (`USLP`).
    pub uslp: Vec<f64>,
    /// ±1: sign relating BL speed and panel tangential speed.
    pub vti: Vec<f64>,
    /// Panel node of each station (0-based, airfoil then wake nodes).
    pub ipan: Vec<usize>,
    /// Newton-system row of each station (0-based).
    pub isys: Vec<usize>,
}

impl Side {
    /// All-zero arrays for `n` stations.
    pub fn zeros(n: usize) -> Self {
        let z = || vec![0.0; n];
        Self {
            xssi: z(),
            uedg: z(),
            uinv: z(),
            uinv_a: z(),
            mass: z(),
            thet: z(),
            dstr: z(),
            ctau: z(),
            delt: z(),
            tstr: z(),
            tau: z(),
            dis: z(),
            ctq: z(),
            uslp: z(),
            vti: z(),
            ipan: vec![0; n],
            isys: vec![0; n],
        }
    }
}

/// The station arrays UPDATE mirrors from the lower wake to the upper wake.
#[derive(Debug, Clone, Copy)]
pub(crate) struct StationValues {
    ctau: f64,
    thet: f64,
    dstr: f64,
    uedg: f64,
    tau: f64,
    dis: f64,
    ctq: f64,
    delt: f64,
    tstr: f64,
}

impl Side {
    pub(crate) fn clone_station(&self, i: usize) -> StationValues {
        StationValues {
            ctau: self.ctau[i],
            thet: self.thet[i],
            dstr: self.dstr[i],
            uedg: self.uedg[i],
            tau: self.tau[i],
            dis: self.dis[i],
            ctq: self.ctq[i],
            delt: self.delt[i],
            tstr: self.tstr[i],
        }
    }

    pub(crate) fn set_station(&mut self, i: usize, v: &StationValues) {
        self.ctau[i] = v.ctau;
        self.thet[i] = v.thet;
        self.dstr[i] = v.dstr;
        self.uedg[i] = v.uedg;
        self.tau[i] = v.tau;
        self.dis[i] = v.dis;
        self.ctq[i] = v.ctq;
        self.delt[i] = v.delt;
        self.tstr[i] = v.tstr;
    }
}

/// Boundary-layer state of both sides and the wake.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryLayer {
    /// Upper side (0) and lower side plus wake (1).
    pub sides: [Side; 2],
    /// TE station index per side (Fortran `IBLTE - 1`).
    pub iblte: [usize; 2],
    /// Number of stations per side, including the dummy station 0 (Fortran `NBL`).
    pub nbl: [usize; 2],
    /// Transition station index per side (Fortran `ITRAN - 1`).
    pub itran: [usize; 2],
    /// Transition ξ per side (`XSSITR`).
    pub xssitr: [f64; 2],
    /// Whether transition was forced (`TFORCE`).
    pub tforce: [bool; 2],
    /// Stagnation panel: stagnation point between nodes `ist` and `ist+1` (Fortran `IST-1`).
    pub ist: usize,
    /// Stagnation point arc length (`SST`).
    pub sst: f64,
    /// dSST/dGAM(IST) (`SST_GO`).
    pub sst_go: f64,
    /// dSST/dGAM(IST+1) (`SST_GP`).
    pub sst_gp: f64,
    /// Number of Newton-system rows (`NSYS`).
    pub nsys: usize,
    /// "Dead air" thickness in the wake behind a blunt TE (`WGAP`, NW).
    pub wgap: Vec<f64>,
    /// Transition x/c per side (`XOCTR`, set by SETBL).
    pub xoctr: [f64; 2],
    /// Transition y/c per side (`YOCTR`).
    pub yoctr: [f64; 2],
    /// Fractional panel index of transition per side (`TINDEX`, the polar file's `Itr`).
    pub tindex: [f64; 2],
}
