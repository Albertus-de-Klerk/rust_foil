//! Boundary-layer closure relations with analytic derivatives. Port of the correlation
//! routines of XFOIL `xblsys.f`: `HKIN`, `HCT`, `HSL`, `HST`, `CFL`, `CFT`, `DIL`, `DILW`,
//! `DIT`.
//!
//! Each returns a value with its partial derivatives; field names follow the Fortran
//! (`hs_hk` is ∂HS/∂HK). Operation order follows the Fortran for bit-identical results.

use crate::fortran::{pow, powi};

/// A correlation `f(HK, RT, M²)` with its partials.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Correlation {
    /// Value.
    pub v: f64,
    /// ∂/∂HK.
    pub hk: f64,
    /// ∂/∂Rθ.
    pub rt: f64,
    /// ∂/∂M².
    pub msq: f64,
}

/// Kinematic shape parameter `HK(H, M²)` (Whitfield). Port of XFOIL `HKIN`.
///
/// Returns `(HK, ∂HK/∂H, ∂HK/∂M²)`.
pub fn hkin(h: f64, msq: f64) -> (f64, f64, f64) {
    let hk = (h - 0.29 * msq) / (1.0 + 0.113 * msq);
    let hk_h = 1.0 / (1.0 + 0.113 * msq);
    let hk_msq = (-0.29 - 0.113 * hk) / (1.0 + 0.113 * msq);
    (hk, hk_h, hk_msq)
}

/// Density shape parameter `H**` (Whitfield). Port of XFOIL `HCT`.
///
/// Returns `(HC, ∂HC/∂HK, ∂HC/∂M²)`.
pub fn hct(hk: f64, msq: f64) -> (f64, f64, f64) {
    let hc = msq * (0.064 / (hk - 0.8) + 0.251);
    let hc_hk = msq * (-0.064 / powi(hk - 0.8, 2));
    let hc_msq = 0.064 / (hk - 0.8) + 0.251;
    (hc, hc_hk, hc_msq)
}

/// Laminar energy shape factor `H*`. Port of XFOIL `HSL`.
pub fn hsl(hk: f64, _rt: f64, _msq: f64) -> Correlation {
    let (hs, hs_hk) = if hk < 4.35 {
        let tmp = hk - 4.35;
        let hs = 0.0111 * powi(tmp, 2) / (hk + 1.0) - 0.0278 * powi(tmp, 3) / (hk + 1.0) + 1.528
            - 0.0002 * powi(tmp * hk, 2);
        let hs_hk = 0.0111 * (2.0 * tmp - powi(tmp, 2) / (hk + 1.0)) / (hk + 1.0)
            - 0.0278 * (3.0 * powi(tmp, 2) - powi(tmp, 3) / (hk + 1.0)) / (hk + 1.0)
            - 0.0002 * 2.0 * tmp * hk * (tmp + hk);
        (hs, hs_hk)
    } else {
        let hs2 = 0.015;
        let hs = hs2 * powi(hk - 4.35, 2) / hk + 1.528;
        let hs_hk = hs2 * 2.0 * (hk - 4.35) / hk - hs2 * powi(hk - 4.35, 2) / powi(hk, 2);
        (hs, hs_hk)
    };
    Correlation {
        v: hs,
        hk: hs_hk,
        rt: 0.0,
        msq: 0.0,
    }
}

/// Minimum turbulent `H*` (`HSMIN`). QFoil: 1.505 (XFOIL 6.99: 1.500), PORTING_PLAN D4.
pub const HSMIN: f64 = 1.505;
/// Separated-branch asymptote coefficient (`DHSINF`). QFoil: 0.04 (XFOIL 6.99: 0.015), D4.
pub const DHSINF: f64 = 0.04;

/// Turbulent energy shape factor `H*`. Port of XFOIL/QFoil `HST` (QFoil constants, D4).
pub fn hst(hk: f64, rt: f64, msq: f64) -> Correlation {
    // limited Rθ dependence for Rθ < 200
    let (ho, ho_rt) = if rt > 400.0 {
        (3.0 + 400.0 / rt, -400.0 / powi(rt, 2))
    } else {
        (4.0, 0.0)
    };
    let (rtz, rtz_rt) = if rt > 200.0 { (rt, 1.0) } else { (200.0, 0.0) };

    let (hs, hs_hk, hs_rt) = if hk < ho {
        // attached branch (arctan(y+) + Schlichting profiles, 29 Nov 91)
        let hr = (ho - hk) / (ho - 1.0);
        let hr_hk = -1.0 / (ho - 1.0);
        let hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt;
        let a = 2.0 - HSMIN - 4.0 / rtz;
        let hs = a * powi(hr, 2) * 1.5 / (hk + 0.5) + HSMIN + 4.0 / rtz;
        let hs_hk =
            -a * powi(hr, 2) * 1.5 / powi(hk + 0.5, 2) + a * hr * 2.0 * 1.5 / (hk + 0.5) * hr_hk;
        let hs_rt = a * hr * 2.0 * 1.5 / (hk + 0.5) * hr_rt
            + (powi(hr, 2) * 1.5 / (hk + 0.5) - 1.0) * 4.0 / powi(rtz, 2) * rtz_rt;
        (hs, hs_hk, hs_rt)
    } else {
        // separated branch
        let grt = rtz.ln();
        let hdif = hk - ho;
        let rtmp = hk - ho + 4.0 / grt;
        let htmp = 0.007 * grt / powi(rtmp, 2) + DHSINF / hk;
        let htmp_hk = -0.014 * grt / powi(rtmp, 3) - DHSINF / powi(hk, 2);
        let htmp_rt = -0.014 * grt / powi(rtmp, 3) * (-ho_rt - 4.0 / powi(grt, 2) / rtz * rtz_rt)
            + 0.007 / powi(rtmp, 2) / rtz * rtz_rt;
        let hs = powi(hdif, 2) * htmp + HSMIN + 4.0 / rtz;
        let hs_hk = hdif * 2.0 * htmp + powi(hdif, 2) * htmp_hk;
        let hs_rt =
            powi(hdif, 2) * htmp_rt - 4.0 / powi(rtz, 2) * rtz_rt + hdif * 2.0 * htmp * (-ho_rt);
        (hs, hs_hk, hs_rt)
    };

    // Whitfield's minor additional compressibility correction
    let fm = 1.0 + 0.014 * msq;
    let hs = (hs + 0.028 * msq) / fm;
    Correlation {
        v: hs,
        hk: hs_hk / fm,
        rt: hs_rt / fm,
        msq: 0.028 / fm - 0.014 * hs / fm,
    }
}

/// Laminar skin friction (Falkner–Skan). Port of XFOIL `CFL`.
pub fn cfl(hk: f64, rt: f64, _msq: f64) -> Correlation {
    let (cf, cf_hk) = if hk < 5.5 {
        let tmp = powi(5.5 - hk, 3) / (hk + 1.0);
        let cf = (0.0727 * tmp - 0.07) / rt;
        let cf_hk = (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / rt;
        (cf, cf_hk)
    } else {
        let tmp = 1.0 - 1.0 / (hk - 4.5);
        let cf = (0.015 * powi(tmp, 2) - 0.07) / rt;
        let cf_hk = (0.015 * tmp * 2.0 / powi(hk - 4.5, 2)) / rt;
        (cf, cf_hk)
    };
    Correlation {
        v: cf,
        hk: cf_hk,
        rt: -cf / rt,
        msq: 0.0,
    }
}

/// XFOIL's truncated ln(10) in CFT. Not `LN_10`: the literal must match for bit parity.
#[allow(clippy::approx_constant)]
const LN10_CFT: f64 = 2.3026;

/// Turbulent skin friction (Coles). Port of XFOIL `CFT`; `cffac` is `/BLPAR/ CFFAC`.
pub fn cft(hk: f64, rt: f64, msq: f64, cffac: f64) -> Correlation {
    const GAM: f64 = 1.4;
    let gm1 = GAM - 1.0;
    let fc = (1.0 + 0.5 * gm1 * msq).sqrt();
    let grt = (rt / fc).ln().max(3.0);
    let gex = -1.74 - 0.31 * hk;
    let arg = (-1.33 * hk).max(-20.0);
    let thk = (4.0 - hk / 0.875).tanh();

    let cfo = cffac * 0.3 * arg.exp() * pow(grt / LN10_CFT, gex);
    let cf = (cfo + 1.1e-4 * (thk - 1.0)) / fc;
    let cf_hk =
        (-1.33 * cfo - 0.31 * (grt / LN10_CFT).ln() * cfo - 1.1e-4 * (1.0 - powi(thk, 2)) / 0.875)
            / fc;
    let cf_rt = gex * cfo / (fc * grt) / rt;
    let cf_msq =
        gex * cfo / (fc * grt) * (-0.25 * gm1 / powi(fc, 2)) - 0.25 * gm1 * cf / powi(fc, 2);
    Correlation {
        v: cf,
        hk: cf_hk,
        rt: cf_rt,
        msq: cf_msq,
    }
}

/// Laminar dissipation `2 CD/H*` (Falkner–Skan). Port of XFOIL `DIL`.
pub fn dil(hk: f64, rt: f64) -> Correlation {
    let (di, di_hk) = if hk < 4.0 {
        let di = (0.00205 * pow(4.0 - hk, 5.5) + 0.207) / rt;
        let di_hk = (-0.00205 * 5.5 * pow(4.0 - hk, 4.5)) / rt;
        (di, di_hk)
    } else {
        let hkb = hk - 4.0;
        let den = 1.0 + 0.02 * powi(hkb, 2);
        let di = (-0.0016 * powi(hkb, 2) / den + 0.207) / rt;
        let di_hk = (-0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * powi(hkb, 2) / powi(den, 2))) / rt;
        (di, di_hk)
    };
    Correlation {
        v: di,
        hk: di_hk,
        rt: -di / rt,
        msq: 0.0,
    }
}

/// Laminar wake dissipation `2 CD/H*`. Port of XFOIL `DILW`.
pub fn dilw(hk: f64, rt: f64) -> Correlation {
    let hs = hsl(hk, rt, 0.0);
    let rcd = 1.10 * powi(1.0 - 1.0 / hk, 2) / hk;
    let rcd_hk = -1.10 * (1.0 - 1.0 / hk) * 2.0 / powi(hk, 3) - rcd / hk;
    let di = 2.0 * rcd / (hs.v * rt);
    let di_hk = 2.0 * rcd_hk / (hs.v * rt) - (di / hs.v) * hs.hk;
    let di_rt = -di / rt - (di / hs.v) * hs.rt;
    Correlation {
        v: di,
        hk: di_hk,
        rt: di_rt,
        msq: 0.0,
    }
}

/// Turbulent dissipation `2 CD/H*` and its partials. Port of XFOIL `DIT`.
///
/// Returns `(DI, ∂/∂HS, ∂/∂US, ∂/∂CF, ∂/∂ST)`. Not called by BLVAR in XFOIL 6.99, which
/// inlines its own turbulent dissipation; ported for completeness.
pub fn dit(hs: f64, us: f64, cf: f64, st: f64) -> (f64, f64, f64, f64, f64) {
    let di = (0.5 * cf * us + st * st * (1.0 - us)) * 2.0 / hs;
    let di_hs = -(0.5 * cf * us + st * st * (1.0 - us)) * 2.0 / powi(hs, 2);
    let di_us = (0.5 * cf - st * st) * 2.0 / hs;
    let di_cf = (0.5 * us) * 2.0 / hs;
    let di_st = (2.0 * st * (1.0 - us)) * 2.0 / hs;
    (di, di_hs, di_us, di_cf, di_st)
}
