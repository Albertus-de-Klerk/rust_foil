//! BL station variables and the local-system kernel (XBL.INC state). Port of XFOIL
//! `BLPRV`, `BLKIN`, `BLVAR`, `BLMID`.
//!
//! XFOIL keeps two stations in COMMON `/V_VAR1/` and `/V_VAR2/` (73 reals each,
//! EQUIVALENCEd to `COM1`/`COM2`). Here they are [`Station`] values `s1`, `s2` in a
//! [`Kernel`]; `COM1 = COM2` becomes `s1 = s2`. Field names drop the station digit:
//! `HK2_T2` is `s2.hk_t`.

use super::closure::{cfl, cft, dil, dilw, hct, hkin, hsl, hst};
use super::transition::AmplificationModel;
use crate::fortran::powi;
use crate::settings::BlParams;

/// Primary and secondary variables of one BL station with their sensitivities
/// (`X, U, T, D, S, AMPL, ...`; see XBL.INC).
#[allow(missing_docs)]
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Station {
    pub x: f64,
    pub u: f64,
    pub t: f64,
    pub d: f64,
    pub s: f64,
    pub ampl: f64,
    pub u_uei: f64,
    pub u_ms: f64,
    pub dw: f64,
    pub h: f64,
    pub h_t: f64,
    pub h_d: f64,
    pub m: f64,
    pub m_u: f64,
    pub m_ms: f64,
    pub r: f64,
    pub r_u: f64,
    pub r_ms: f64,
    pub v: f64,
    pub v_u: f64,
    pub v_ms: f64,
    pub v_re: f64,
    pub hk: f64,
    pub hk_u: f64,
    pub hk_t: f64,
    pub hk_d: f64,
    pub hk_ms: f64,
    pub hs: f64,
    pub hs_u: f64,
    pub hs_t: f64,
    pub hs_d: f64,
    pub hs_ms: f64,
    pub hs_re: f64,
    pub hc: f64,
    pub hc_u: f64,
    pub hc_t: f64,
    pub hc_d: f64,
    pub hc_ms: f64,
    pub rt: f64,
    pub rt_u: f64,
    pub rt_t: f64,
    pub rt_ms: f64,
    pub rt_re: f64,
    pub cf: f64,
    pub cf_u: f64,
    pub cf_t: f64,
    pub cf_d: f64,
    pub cf_ms: f64,
    pub cf_re: f64,
    pub di: f64,
    pub di_u: f64,
    pub di_t: f64,
    pub di_d: f64,
    pub di_s: f64,
    pub di_ms: f64,
    pub di_re: f64,
    pub us: f64,
    pub us_u: f64,
    pub us_t: f64,
    pub us_d: f64,
    pub us_ms: f64,
    pub us_re: f64,
    pub cq: f64,
    pub cq_u: f64,
    pub cq_t: f64,
    pub cq_d: f64,
    pub cq_ms: f64,
    pub cq_re: f64,
    pub de: f64,
    pub de_u: f64,
    pub de_t: f64,
    pub de_d: f64,
    pub de_ms: f64,
}

/// Flow regime of BLVAR/BLMID (`ITYP` 1–3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Regime {
    /// Laminar (ITYP = 1).
    Laminar,
    /// Turbulent (ITYP = 2).
    Turbulent,
    /// Turbulent wake (ITYP = 3).
    Wake,
}

/// Freestream-dependent constants of the BL equations (COMMON `/V_VAR/`), set in SETBL.
#[allow(missing_docs)]
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct KernelParams {
    pub dwte: f64,
    pub qinfbl: f64,
    pub tkbl: f64,
    pub tkbl_ms: f64,
    pub rstbl: f64,
    pub rstbl_ms: f64,
    pub hstinv: f64,
    pub hstinv_ms: f64,
    pub reybl: f64,
    pub reybl_ms: f64,
    pub reybl_re: f64,
    pub gambl: f64,
    pub gm1bl: f64,
    /// Sutherland viscosity ratio. Set only by XFOIL's plotting code, so 0 in analysis.
    pub hvrat: f64,
    pub bule: f64,
    pub xiforc: f64,
    pub amcrit: f64,
}

impl KernelParams {
    /// Port of the freestream part of XFOIL `SETBL` (after MRCL/COMSET for the current CL).
    pub fn new(
        minf: f64,
        reinf: f64,
        tklam: f64,
        tkl_msq: f64,
        gamma: f64,
        qinf: f64,
        dwte: f64,
    ) -> Self {
        let gm1bl = gamma - 1.0;
        let qinfbl = qinf;
        let rstbl = (1.0 + 0.5 * gm1bl * powi(minf, 2)).powf(1.0 / gm1bl);
        let rstbl_ms = 0.5 * rstbl / (1.0 + 0.5 * gm1bl * powi(minf, 2));
        let hstinv = gm1bl * powi(minf / qinfbl, 2) / (1.0 + 0.5 * gm1bl * powi(minf, 2));
        let hstinv_ms = gm1bl * powi(1.0 / qinfbl, 2) / (1.0 + 0.5 * gm1bl * powi(minf, 2))
            - 0.5 * gm1bl * hstinv / (1.0 + 0.5 * gm1bl * powi(minf, 2));
        let hvrat = 0.0;
        let herat = 1.0 - 0.5 * powi(qinfbl, 2) * hstinv;
        let herat_ms = -0.5 * powi(qinfbl, 2) * hstinv_ms;
        let reybl = reinf * powi(herat, 3).sqrt() * (1.0 + hvrat) / (herat + hvrat);
        let reybl_re = powi(herat, 3).sqrt() * (1.0 + hvrat) / (herat + hvrat);
        let reybl_ms = reybl * (1.5 / herat - 1.0 / (herat + hvrat)) * herat_ms;
        Self {
            dwte,
            qinfbl,
            tkbl: tklam,
            tkbl_ms: tkl_msq,
            rstbl,
            rstbl_ms,
            hstinv,
            hstinv_ms,
            reybl,
            reybl_ms,
            reybl_re,
            gambl: gamma,
            gm1bl,
            hvrat,
            bule: 0.0,
            xiforc: 0.0,
            amcrit: 0.0,
        }
    }
}

/// Midpoint skin friction and its sensitivities (`/V_VARA/` CFM part).
#[allow(missing_docs)]
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Midpoint {
    pub cfm: f64,
    pub cfm_ms: f64,
    pub cfm_re: f64,
    pub cfm_u1: f64,
    pub cfm_t1: f64,
    pub cfm_d1: f64,
    pub cfm_u2: f64,
    pub cfm_t2: f64,
    pub cfm_d2: f64,
}

/// Transition location and its sensitivities (`/V_VARA/` XT part).
#[allow(missing_docs)]
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct TransitionPoint {
    pub xt: f64,
    pub xt_a1: f64,
    pub xt_ms: f64,
    pub xt_re: f64,
    pub xt_xf: f64,
    pub xt_x1: f64,
    pub xt_t1: f64,
    pub xt_d1: f64,
    pub xt_u1: f64,
    pub xt_x2: f64,
    pub xt_t2: f64,
    pub xt_d2: f64,
    pub xt_u2: f64,
}

/// The local 4×5 Newton system of one interval (`/V_SYS/`). Rows are the equations
/// (amplification or shear lag, momentum, shape, and the 4th slot used by the march);
/// columns are the variables (N or Cτ, θ, δ*, Ue, ξ).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct LocalSystem {
    /// ∂R/∂(station 1 variables), `VS1(K,L)` = `vs1[K-1][L-1]`.
    pub vs1: [[f64; 5]; 4],
    /// ∂R/∂(station 2 variables).
    pub vs2: [[f64; 5]; 4],
    /// Residuals (`VSREZ`).
    pub vsrez: [f64; 4],
    /// ∂R/∂Re (`VSR`).
    pub vsr: [f64; 4],
    /// ∂R/∂M² (`VSM`).
    pub vsm: [f64; 4],
    /// ∂R/∂ξ_forced (`VSX`).
    pub vsx: [f64; 4],
}

impl LocalSystem {
    /// Zeroes all coefficients (the reset at the top of BLDIF/TESYS).
    pub fn clear(&mut self) {
        *self = Self::default();
    }
}

/// Station-2 rows of a [`LocalSystem`] viewed as an `(i, j)`-indexed matrix for GAUSS.
pub(crate) struct Rows<'a>(pub &'a mut [[f64; 5]; 4]);

impl std::ops::Index<(usize, usize)> for Rows<'_> {
    type Output = f64;
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        &self.0[i][j]
    }
}
impl std::ops::IndexMut<(usize, usize)> for Rows<'_> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut f64 {
        &mut self.0[i][j]
    }
}

/// The BL equation kernel: the two stations plus all state XBL.INC keeps in COMMON.
#[derive(Debug, Clone, PartialEq)]
pub struct Kernel {
    /// Upstream station ("1").
    pub s1: Station,
    /// Current station ("2").
    pub s2: Station,
    /// Freestream constants.
    pub p: KernelParams,
    /// Model constants (`/BLPAR/`).
    pub bl: BlParams,
    /// e^N amplification model (`IDAMPV`).
    pub model: AmplificationModel,
    /// Similarity (stagnation) station.
    pub simi: bool,
    /// Transition occurs in the current interval.
    pub tran: bool,
    /// Current interval is turbulent.
    pub turb: bool,
    /// Current interval is in the wake.
    pub wake: bool,
    /// Transition is forced (trip).
    pub trforc: bool,
    /// Transition is free.
    pub trfree: bool,
    /// Midpoint skin friction.
    pub mid: Midpoint,
    /// Transition point.
    pub xt: TransitionPoint,
    /// Local Newton system.
    pub sys: LocalSystem,
}

impl Kernel {
    /// A kernel with zeroed stations.
    pub fn new(p: KernelParams, bl: BlParams, model: AmplificationModel) -> Self {
        Self {
            s1: Station::default(),
            s2: Station::default(),
            p,
            bl,
            model,
            simi: false,
            tran: false,
            turb: false,
            wake: false,
            trforc: false,
            trfree: false,
            mid: Midpoint::default(),
            xt: TransitionPoint::default(),
            sys: LocalSystem::default(),
        }
    }

    /// Sets the primary station-2 variables. Port of XFOIL `BLPRV`.
    pub fn blprv(
        &mut self,
        xsi: f64,
        ami: f64,
        cti: f64,
        thi: f64,
        dsi: f64,
        dswaki: f64,
        uei: f64,
    ) {
        let p = &self.p;
        let s = &mut self.s2;
        s.x = xsi;
        s.ampl = ami;
        s.s = cti;
        s.t = thi;
        s.d = dsi - dswaki;
        s.dw = dswaki;
        let den = 1.0 - p.tkbl * powi(uei / p.qinfbl, 2);
        s.u = uei * (1.0 - p.tkbl) / den;
        s.u_uei = (1.0 + p.tkbl * (2.0 * s.u * uei / powi(p.qinfbl, 2) - 1.0)) / den;
        s.u_ms = (s.u * powi(uei / p.qinfbl, 2) - uei) * p.tkbl_ms / den;
    }

    /// Turbulence-independent secondary variables of station 2. Port of XFOIL `BLKIN`.
    pub fn blkin(&mut self) {
        let p = &self.p;
        let s = &mut self.s2;
        let u = s.u;

        // edge Mach number squared
        s.m = u * u * p.hstinv / (p.gm1bl * (1.0 - 0.5 * u * u * p.hstinv));
        let tr = 1.0 + 0.5 * p.gm1bl * s.m;
        s.m_u = 2.0 * s.m * tr / u;
        s.m_ms = u * u * tr / (p.gm1bl * (1.0 - 0.5 * u * u * p.hstinv)) * p.hstinv_ms;

        // edge static density (isentropic)
        s.r = p.rstbl * tr.powf(-1.0 / p.gm1bl);
        s.r_u = -s.r / tr * 0.5 * s.m_u;
        s.r_ms = -s.r / tr * 0.5 * s.m_ms + p.rstbl_ms * tr.powf(-1.0 / p.gm1bl);

        // shape parameter
        s.h = s.d / s.t;
        s.h_d = 1.0 / s.t;
        s.h_t = -s.h / s.t;

        // static/stagnation enthalpy ratio
        let herat = 1.0 - 0.5 * u * u * p.hstinv;
        let he_u = -u * p.hstinv;
        let he_ms = -0.5 * u * u * p.hstinv_ms;

        // molecular viscosity
        s.v = powi(herat, 3).sqrt() * (1.0 + p.hvrat) / (herat + p.hvrat) / p.reybl;
        let v_he = s.v * (1.5 / herat - 1.0 / (herat + p.hvrat));
        s.v_u = v_he * he_u;
        s.v_ms = -s.v / p.reybl * p.reybl_ms + v_he * he_ms;
        s.v_re = -s.v / p.reybl * p.reybl_re;

        // kinematic shape parameter
        let (hk, hk_h, hk_m) = hkin(s.h, s.m);
        s.hk = hk;
        s.hk_u = hk_m * s.m_u;
        s.hk_t = hk_h * s.h_t;
        s.hk_d = hk_h * s.h_d;
        s.hk_ms = hk_m * s.m_ms;

        // momentum-thickness Reynolds number
        s.rt = s.r * u * s.t / s.v;
        s.rt_u = s.rt * (1.0 / u + s.r_u / s.r - s.v_u / s.v);
        s.rt_t = s.rt / s.t;
        s.rt_ms = s.rt * (s.r_ms / s.r - s.v_ms / s.v);
        s.rt_re = s.rt * (-s.v_re / s.v);
    }

    /// All secondary variables of station 2 and their sensitivities. Port of XFOIL `BLVAR`.
    pub fn blvar(&mut self, ityp: Regime) {
        let bl = &self.bl;
        let s = &mut self.s2;

        // HK is clamped, its sensitivities are not (as in the Fortran)
        s.hk = if ityp == Regime::Wake {
            s.hk.max(1.00005)
        } else {
            s.hk.max(1.05000)
        };

        // density thickness shape parameter H**
        let (hc, hc_hk, hc_m) = hct(s.hk, s.m);
        s.hc = hc;
        s.hc_u = hc_hk * s.hk_u + hc_m * s.m_u;
        s.hc_t = hc_hk * s.hk_t;
        s.hc_d = hc_hk * s.hk_d;
        s.hc_ms = hc_hk * s.hk_ms + hc_m * s.m_ms;

        // KE thickness shape parameter from H - H* correlations
        let c = if ityp == Regime::Laminar {
            hsl(s.hk, s.rt, s.m)
        } else {
            hst(s.hk, s.rt, s.m)
        };
        s.hs = c.v;
        s.hs_u = c.hk * s.hk_u + c.rt * s.rt_u + c.msq * s.m_u;
        s.hs_t = c.hk * s.hk_t + c.rt * s.rt_t;
        s.hs_d = c.hk * s.hk_d;
        s.hs_ms = c.hk * s.hk_ms + c.rt * s.rt_ms + c.msq * s.m_ms;
        s.hs_re = c.rt * s.rt_re;

        // normalised slip velocity Us
        s.us = 0.5 * s.hs * (1.0 - (s.hk - 1.0) / (bl.gbcon * s.h));
        let us_hs = 0.5 * (1.0 - (s.hk - 1.0) / (bl.gbcon * s.h));
        let us_hk = 0.5 * s.hs * (-1.0 / (bl.gbcon * s.h));
        let us_h = 0.5 * s.hs * (s.hk - 1.0) / (bl.gbcon * powi(s.h, 2));
        s.us_u = us_hs * s.hs_u + us_hk * s.hk_u;
        s.us_t = us_hs * s.hs_t + us_hk * s.hk_t + us_h * s.h_t;
        s.us_d = us_hs * s.hs_d + us_hk * s.hk_d + us_h * s.h_d;
        s.us_ms = us_hs * s.hs_ms + us_hk * s.hk_ms;
        s.us_re = us_hs * s.hs_re;
        let clamp_us = |s: &mut Station, v: f64| {
            s.us = v;
            s.us_u = 0.0;
            s.us_t = 0.0;
            s.us_d = 0.0;
            s.us_ms = 0.0;
            s.us_re = 0.0;
        };
        if ityp != Regime::Wake && s.us > 0.95 {
            clamp_us(s, 0.98);
        }
        if ityp == Regime::Wake && s.us > 0.99995 {
            clamp_us(s, 0.99995);
        }

        // equilibrium wake-layer shear coefficient (Ctau)eq^1/2  (12 Oct 94)
        let (mut hkc, mut hkc_hk, mut hkc_rt) = (s.hk - 1.0, 1.0, 0.0);
        if ityp == Regime::Turbulent {
            let gcc = bl.gccon;
            hkc = s.hk - 1.0 - gcc / s.rt;
            hkc_hk = 1.0;
            hkc_rt = gcc / powi(s.rt, 2);
            if hkc < 0.01 {
                hkc = 0.01;
                hkc_hk = 0.0;
                hkc_rt = 0.0;
            }
        }
        let ctcon = bl.ctcon;
        let hkb = s.hk - 1.0;
        let usb = 1.0 - s.us;
        let (hs, h, hk) = (s.hs, s.h, s.hk);
        s.cq = (ctcon * hs * hkb * powi(hkc, 2) / (usb * h * powi(hk, 2))).sqrt();
        let cq = s.cq;
        let cq_hs = ctcon * hkb * powi(hkc, 2) / (usb * h * powi(hk, 2)) * 0.5 / cq;
        let cq_us = ctcon * hs * hkb * powi(hkc, 2) / (usb * h * powi(hk, 2)) / usb * 0.5 / cq;
        let cq_hk = ctcon * hs * powi(hkc, 2) / (usb * h * powi(hk, 2)) * 0.5 / cq
            - ctcon * hs * hkb * powi(hkc, 2) / (usb * h * powi(hk, 3)) * 2.0 * 0.5 / cq
            + ctcon * hs * hkb * hkc / (usb * h * powi(hk, 2)) * 2.0 * 0.5 / cq * hkc_hk;
        let cq_rt = ctcon * hs * hkb * hkc / (usb * h * powi(hk, 2)) * 2.0 * 0.5 / cq * hkc_rt;
        let cq_h = -ctcon * hs * hkb * powi(hkc, 2) / (usb * h * powi(hk, 2)) / h * 0.5 / cq;
        s.cq_u = cq_hs * s.hs_u + cq_us * s.us_u + cq_hk * s.hk_u;
        s.cq_t = cq_hs * s.hs_t + cq_us * s.us_t + cq_hk * s.hk_t;
        s.cq_d = cq_hs * s.hs_d + cq_us * s.us_d + cq_hk * s.hk_d;
        s.cq_ms = cq_hs * s.hs_ms + cq_us * s.us_ms + cq_hk * s.hk_ms;
        s.cq_re = cq_hs * s.hs_re + cq_us * s.us_re;
        s.cq_u += cq_rt * s.rt_u;
        s.cq_t = s.cq_t + cq_h * s.h_t + cq_rt * s.rt_t;
        s.cq_d += cq_h * s.h_d;
        s.cq_ms += cq_rt * s.rt_ms;
        s.cq_re += cq_rt * s.rt_re;

        // skin friction
        let cf = match ityp {
            Regime::Wake => Default::default(),
            Regime::Laminar => cfl(s.hk, s.rt, s.m),
            Regime::Turbulent => {
                let t = cft(s.hk, s.rt, s.m, bl.cffac);
                let l = cfl(s.hk, s.rt, s.m);
                // laminar Cf above turbulent only for unreasonably small Rθ
                if l.v > t.v { l } else { t }
            }
        };
        s.cf = cf.v;
        s.cf_u = cf.hk * s.hk_u + cf.rt * s.rt_u + cf.msq * s.m_u;
        s.cf_t = cf.hk * s.hk_t + cf.rt * s.rt_t;
        s.cf_d = cf.hk * s.hk_d;
        s.cf_ms = cf.hk * s.hk_ms + cf.rt * s.rt_ms + cf.msq * s.m_ms;
        s.cf_re = cf.rt * s.rt_re;

        // dissipation function 2 CD / H*
        match ityp {
            Regime::Laminar => {
                let d = dil(s.hk, s.rt);
                s.di = d.v;
                s.di_u = d.hk * s.hk_u + d.rt * s.rt_u;
                s.di_t = d.hk * s.hk_t + d.rt * s.rt_t;
                s.di_d = d.hk * s.hk_d;
                s.di_s = 0.0;
                s.di_ms = d.hk * s.hk_ms + d.rt * s.rt_ms;
                s.di_re = d.rt * s.rt_re;
            }
            Regime::Turbulent => {
                // turbulent wall contribution
                let ct = cft(s.hk, s.rt, s.m, bl.cffac);
                let cf2t_u = ct.hk * s.hk_u + ct.rt * s.rt_u + ct.msq * s.m_u;
                let cf2t_t = ct.hk * s.hk_t + ct.rt * s.rt_t;
                let cf2t_d = ct.hk * s.hk_d;
                let cf2t_ms = ct.hk * s.hk_ms + ct.rt * s.rt_ms + ct.msq * s.m_ms;
                let cf2t_re = ct.rt * s.rt_re;

                s.di = (0.5 * ct.v * s.us) * 2.0 / s.hs;
                let di_hs = -(0.5 * ct.v * s.us) * 2.0 / powi(s.hs, 2);
                let di_us = (0.5 * ct.v) * 2.0 / s.hs;
                let di_cf2t = (0.5 * s.us) * 2.0 / s.hs;
                s.di_s = 0.0;
                s.di_u = di_hs * s.hs_u + di_us * s.us_u + di_cf2t * cf2t_u;
                s.di_t = di_hs * s.hs_t + di_us * s.us_t + di_cf2t * cf2t_t;
                s.di_d = di_hs * s.hs_d + di_us * s.us_d + di_cf2t * cf2t_d;
                s.di_ms = di_hs * s.hs_ms + di_us * s.us_ms + di_cf2t * cf2t_ms;
                s.di_re = di_hs * s.hs_re + di_us * s.us_re + di_cf2t * cf2t_re;

                // minimum Hk for the wake layer to still exist; correct wall dissipation
                let grt = s.rt.ln();
                let hmin = 1.0 + 2.1 / grt;
                let hm_rt = -(2.1 / powi(grt, 2)) / s.rt;
                let fl = (s.hk - 1.0) / (hmin - 1.0);
                let fl_hk = 1.0 / (hmin - 1.0);
                let fl_rt = (-fl / (hmin - 1.0)) * hm_rt;
                let tfl = fl.tanh();
                let dfac = 0.5 + 0.5 * tfl;
                let df_fl = 0.5 * (1.0 - powi(tfl, 2));
                let df_hk = df_fl * fl_hk;
                let df_rt = df_fl * fl_rt;
                s.di_s *= dfac;
                s.di_u = s.di_u * dfac + s.di * (df_hk * s.hk_u + df_rt * s.rt_u);
                s.di_t = s.di_t * dfac + s.di * (df_hk * s.hk_t + df_rt * s.rt_t);
                s.di_d = s.di_d * dfac + s.di * (df_hk * s.hk_d);
                s.di_ms = s.di_ms * dfac + s.di * (df_hk * s.hk_ms + df_rt * s.rt_ms);
                s.di_re = s.di_re * dfac + s.di * (df_rt * s.rt_re);
                s.di *= dfac;
            }
            Regime::Wake => {
                s.di = 0.0;
                s.di_s = 0.0;
                s.di_u = 0.0;
                s.di_t = 0.0;
                s.di_d = 0.0;
                s.di_ms = 0.0;
                s.di_re = 0.0;
            }
        }

        // turbulent outer-layer contribution
        if ityp != Regime::Laminar {
            let dd = powi(s.s, 2) * (0.995 - s.us) * 2.0 / s.hs;
            let dd_hs = -powi(s.s, 2) * (0.995 - s.us) * 2.0 / powi(s.hs, 2);
            let dd_us = -powi(s.s, 2) * 2.0 / s.hs;
            let dd_s = s.s * 2.0 * (0.995 - s.us) * 2.0 / s.hs;
            s.di += dd;
            s.di_s = dd_s;
            s.di_u = s.di_u + dd_hs * s.hs_u + dd_us * s.us_u;
            s.di_t = s.di_t + dd_hs * s.hs_t + dd_us * s.us_t;
            s.di_d = s.di_d + dd_hs * s.hs_d + dd_us * s.us_d;
            s.di_ms = s.di_ms + dd_hs * s.hs_ms + dd_us * s.us_ms;
            s.di_re = s.di_re + dd_hs * s.hs_re + dd_us * s.us_re;

            // laminar stress contribution to the outer layer
            let dd = 0.15 * powi(0.995 - s.us, 2) / s.rt * 2.0 / s.hs;
            let dd_us = -0.15 * (0.995 - s.us) * 2.0 / s.rt * 2.0 / s.hs;
            let dd_hs = -dd / s.hs;
            let dd_rt = -dd / s.rt;
            s.di += dd;
            s.di_u = s.di_u + dd_hs * s.hs_u + dd_us * s.us_u + dd_rt * s.rt_u;
            s.di_t = s.di_t + dd_hs * s.hs_t + dd_us * s.us_t + dd_rt * s.rt_t;
            s.di_d = s.di_d + dd_hs * s.hs_d + dd_us * s.us_d;
            s.di_ms = s.di_ms + dd_hs * s.hs_ms + dd_us * s.us_ms + dd_rt * s.rt_ms;
            s.di_re = s.di_re + dd_hs * s.hs_re + dd_us * s.us_re + dd_rt * s.rt_re;
        }

        // laminar CD above turbulent (only for unreasonably small Rθ): use laminar
        let laminar_floor = match ityp {
            Regime::Turbulent => Some(dil(s.hk, s.rt)),
            Regime::Wake => Some(dilw(s.hk, s.rt)),
            Regime::Laminar => None,
        };
        if let Some(l) = laminar_floor
            && l.v > s.di
        {
            s.di = l.v;
            s.di_s = 0.0;
            s.di_u = l.hk * s.hk_u + l.rt * s.rt_u;
            s.di_t = l.hk * s.hk_t + l.rt * s.rt_t;
            s.di_d = l.hk * s.hk_d;
            s.di_ms = l.hk * s.hk_ms + l.rt * s.rt_ms;
            s.di_re = l.rt * s.rt_re;
        }

        if ityp == Regime::Wake {
            // double dissipation for the two wake halves
            s.di *= 2.0;
            s.di_s *= 2.0;
            s.di_u *= 2.0;
            s.di_t *= 2.0;
            s.di_d *= 2.0;
            s.di_ms *= 2.0;
            s.di_re *= 2.0;
        }

        // BL thickness Delta from simplified Green's correlation
        s.de = (3.15 + 1.72 / (s.hk - 1.0)) * s.t + s.d;
        let de_hk = (-1.72 / powi(s.hk - 1.0, 2)) * s.t;
        s.de_u = de_hk * s.hk_u;
        s.de_t = de_hk * s.hk_t + (3.15 + 1.72 / (s.hk - 1.0));
        s.de_d = de_hk * s.hk_d + 1.0;
        s.de_ms = de_hk * s.hk_ms;
        const HDMAX: f64 = 12.0;
        if s.de > HDMAX * s.t {
            s.de = HDMAX * s.t;
            s.de_u = 0.0;
            s.de_t = HDMAX;
            s.de_d = 0.0;
            s.de_ms = 0.0;
        }
    }

    /// Midpoint skin friction `CFM`. Port of XFOIL `BLMID`.
    ///
    /// At the similarity station the station-1 kinematic variables are first set from
    /// station 2, as in the Fortran.
    pub fn blmid(&mut self, ityp: Regime) {
        if self.simi {
            let (s1, s2) = (&mut self.s1, &self.s2);
            s1.hk = s2.hk;
            s1.hk_t = s2.hk_t;
            s1.hk_d = s2.hk_d;
            s1.hk_u = s2.hk_u;
            s1.hk_ms = s2.hk_ms;
            s1.rt = s2.rt;
            s1.rt_t = s2.rt_t;
            s1.rt_u = s2.rt_u;
            s1.rt_ms = s2.rt_ms;
            s1.rt_re = s2.rt_re;
            s1.m = s2.m;
            s1.m_u = s2.m_u;
            s1.m_ms = s2.m_ms;
        }
        let (s1, s2) = (&self.s1, &self.s2);
        let hka = 0.5 * (s1.hk + s2.hk);
        let rta = 0.5 * (s1.rt + s2.rt);
        let ma = 0.5 * (s1.m + s2.m);
        let c = match ityp {
            Regime::Wake => Default::default(),
            Regime::Laminar => cfl(hka, rta, ma),
            Regime::Turbulent => {
                let t = cft(hka, rta, ma, self.bl.cffac);
                let l = cfl(hka, rta, ma);
                if l.v > t.v { l } else { t }
            }
        };
        let (cfm_hka, cfm_rta, cfm_ma) = (c.hk, c.rt, c.msq);
        self.mid = Midpoint {
            cfm: c.v,
            cfm_u1: 0.5 * (cfm_hka * s1.hk_u + cfm_ma * s1.m_u + cfm_rta * s1.rt_u),
            cfm_t1: 0.5 * (cfm_hka * s1.hk_t + cfm_rta * s1.rt_t),
            cfm_d1: 0.5 * (cfm_hka * s1.hk_d),
            cfm_u2: 0.5 * (cfm_hka * s2.hk_u + cfm_ma * s2.m_u + cfm_rta * s2.rt_u),
            cfm_t2: 0.5 * (cfm_hka * s2.hk_t + cfm_rta * s2.rt_t),
            cfm_d2: 0.5 * (cfm_hka * s2.hk_d),
            cfm_ms: 0.5
                * (cfm_hka * s1.hk_ms
                    + cfm_ma * s1.m_ms
                    + cfm_rta * s1.rt_ms
                    + cfm_hka * s2.hk_ms
                    + cfm_ma * s2.m_ms
                    + cfm_rta * s2.rt_ms),
            cfm_re: 0.5 * (cfm_rta * s1.rt_re + cfm_rta * s2.rt_re),
        };
    }
}
