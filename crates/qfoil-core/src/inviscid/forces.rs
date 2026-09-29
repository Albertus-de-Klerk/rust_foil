//! Compressibility state and force integration. Port of XFOIL `MRCL`, `COMSET`, `CPCALC`,
//! `CLCALC`.

use crate::fortran::{pow, powi};
use crate::settings::{FlowConditions, MachType, ReynoldsType};

/// Freestream state for the current CL (`MINF`, `REINF`, Karman–Tsien terms, sonic values).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Compressibility {
    /// Mach number (`MINF`).
    pub minf: f64,
    /// dM/dCL (`MINF_CL`).
    pub minf_cl: f64,
    /// Reynolds number (`REINF`).
    pub reinf: f64,
    /// dRe/dCL (`REINF_CL`).
    pub reinf_cl: f64,
    /// Karman–Tsien parameter (`TKLAM`).
    pub tklam: f64,
    /// dTKLAM/d(M²) (`TKL_MSQ`).
    pub tkl_msq: f64,
    /// Sonic Cp (`CPSTAR`).
    pub cpstar: f64,
    /// Sonic speed (`QSTAR`).
    pub qstar: f64,
}

/// Sets Mach and Reynolds numbers for lift coefficient `cls`. Port of XFOIL `MRCL`.
///
/// Returns `(M, dM/dCL, Re, dRe/dCL)`. The Fortran's warnings about limiting M to 0.99 and
/// Re to 100·Re₁ are not printed; the limits are applied.
pub fn mrcl(cls: f64, flow: &FlowConditions) -> (f64, f64, f64, f64) {
    let cla = cls.max(0.000001);
    let (mut minf, mut m_cls) = match flow.mach_type {
        MachType::Fixed | MachType::Fixed3 => (flow.mach, 0.0),
        MachType::InverseSqrtCl => {
            let m = flow.mach / cla.sqrt();
            (m, -0.5 * m / cla)
        }
    };
    let (mut reinf, mut r_cls) = match flow.reynolds_type {
        ReynoldsType::Fixed => (flow.reynolds, 0.0),
        ReynoldsType::InverseSqrtCl => {
            let r = flow.reynolds / cla.sqrt();
            (r, -0.5 * r / cla)
        }
        ReynoldsType::InverseCl => {
            let r = flow.reynolds / cla;
            (r, -r / cla)
        }
    };
    if minf >= 0.99 {
        minf = 0.99;
        m_cls = 0.0;
    }
    let rrat = if flow.reynolds > 0.0 {
        reinf / flow.reynolds
    } else {
        1.0
    };
    if rrat > 100.0 {
        reinf = flow.reynolds * 100.0;
        r_cls = 0.0;
    }
    (minf, m_cls, reinf, r_cls)
}

impl Compressibility {
    /// Mach/Re for `cl` (MRCL) and the dependent terms (COMSET).
    pub fn at_cl(cl: f64, flow: &FlowConditions, qinf: f64) -> Self {
        let (minf, minf_cl, reinf, reinf_cl) = mrcl(cl, flow);
        let mut c = Self {
            minf,
            minf_cl,
            reinf,
            reinf_cl,
            ..Self::default()
        };
        c.comset(flow.gamma, qinf);
        c
    }

    /// Karman–Tsien parameter and sonic conditions. Port of XFOIL `COMSET`.
    pub fn comset(&mut self, gamma: f64, qinf: f64) {
        let gamm1 = gamma - 1.0;
        let minf = self.minf;
        let beta = (1.0 - powi(minf, 2)).sqrt();
        let beta_msq = -0.5 / beta;
        self.tklam = powi(minf, 2) / powi(1.0 + beta, 2);
        self.tkl_msq = 1.0 / powi(1.0 + beta, 2) - 2.0 * self.tklam / (1.0 + beta) * beta_msq;
        if minf == 0.0 {
            self.cpstar = -999.0;
            self.qstar = 999.0;
        } else {
            self.cpstar = 2.0 / (gamma * powi(minf, 2))
                * (pow(
                    (1.0 + 0.5 * gamm1 * powi(minf, 2)) / (1.0 + 0.5 * gamm1),
                    gamma / gamm1,
                ) - 1.0);
            self.qstar =
                qinf / minf * ((1.0 + 0.5 * gamm1 * powi(minf, 2)) / (1.0 + 0.5 * gamm1)).sqrt();
        }
    }
}

/// Karman–Tsien Cp from surface speed. Port of XFOIL `CPCALC`.
///
/// Returns `false` if the compressibility correction became invalid (denominator ≤ 0),
/// where the Fortran prints a warning.
pub fn cpcalc(q: &[f64], qinf: f64, minf: f64, cp: &mut [f64]) -> bool {
    let beta = (1.0 - powi(minf, 2)).sqrt();
    let bfac = 0.5 * powi(minf, 2) / (1.0 + beta);
    let mut ok = true;
    for (c, &qi) in cp.iter_mut().zip(q) {
        let cpinc = 1.0 - powi(qi / qinf, 2);
        let den = beta + bfac * cpinc;
        *c = cpinc / den;
        if den <= 0.0 {
            ok = false;
        }
    }
    ok
}

/// Integrated forces from [`clcalc`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Forces {
    /// Lift coefficient.
    pub cl: f64,
    /// Moment coefficient about the reference point.
    pub cm: f64,
    /// Pressure drag from surface pressure integration.
    pub cdp: f64,
    /// dCL/dα.
    pub cl_alf: f64,
    /// dCL/d(M²).
    pub cl_msq: f64,
}

/// Integrates surface pressures for CL, CM, CDP. Port of XFOIL `CLCALC`.
pub fn clcalc(
    x: &[f64],
    y: &[f64],
    gam: &[f64],
    gam_a: &[f64],
    alfa: f64,
    minf: f64,
    qinf: f64,
    (xref, yref): (f64, f64),
) -> Forces {
    let n = x.len();
    let sa = alfa.sin();
    let ca = alfa.cos();
    let beta = (1.0 - powi(minf, 2)).sqrt();
    let beta_msq = -0.5 / beta;
    let bfac = 0.5 * powi(minf, 2) / (1.0 + beta);
    let bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

    let mut f = Forces::default();

    let i = 0;
    let cginc = 1.0 - powi(gam[i] / qinf, 2);
    let mut cpg1 = cginc / (beta + bfac * cginc);
    let mut cpg1_msq = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
    let cpi_gam = -2.0 * gam[i] / powi(qinf, 2);
    let cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
    let mut cpg1_alf = cpc_cpi * cpi_gam * gam_a[i];

    for i in 0..n {
        let ip = if i == n - 1 { 0 } else { i + 1 };
        let cginc = 1.0 - powi(gam[ip] / qinf, 2);
        let cpg2 = cginc / (beta + bfac * cginc);
        let cpg2_msq = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
        let cpi_gam = -2.0 * gam[ip] / powi(qinf, 2);
        let cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
        let cpg2_alf = cpc_cpi * cpi_gam * gam_a[ip];

        let dx = (x[ip] - x[i]) * ca + (y[ip] - y[i]) * sa;
        let dy = (y[ip] - y[i]) * ca - (x[ip] - x[i]) * sa;
        let dg = cpg2 - cpg1;
        let ax = (0.5 * (x[ip] + x[i]) - xref) * ca + (0.5 * (y[ip] + y[i]) - yref) * sa;
        let ay = (0.5 * (y[ip] + y[i]) - yref) * ca - (0.5 * (x[ip] + x[i]) - xref) * sa;
        let ag = 0.5 * (cpg2 + cpg1);
        let dx_alf = -(x[ip] - x[i]) * sa + (y[ip] - y[i]) * ca;
        let ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
        let ag_msq = 0.5 * (cpg2_msq + cpg1_msq);

        f.cl += dx * ag;
        f.cdp -= dy * ag;
        f.cm = f.cm - dx * (ag * ax + dg * dx / 12.0) - dy * (ag * ay + dg * dy / 12.0);
        f.cl_alf = f.cl_alf + dx * ag_alf + ag * dx_alf;
        f.cl_msq += dx * ag_msq;

        cpg1 = cpg2;
        cpg1_alf = cpg2_alf;
        cpg1_msq = cpg2_msq;
    }
    f
}
