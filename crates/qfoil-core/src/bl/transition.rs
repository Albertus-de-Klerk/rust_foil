//! e^N amplification rates. Port of XFOIL `DAMPL`, `DAMPL2` and `AXSET` (xblsys.f).

use crate::fortran::powi;

/// Spatial amplification rate `AX = dN/dx` with partials.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Amplification {
    /// dN/dx.
    pub ax: f64,
    /// ∂AX/∂HK.
    pub hk: f64,
    /// ∂AX/∂θ.
    pub th: f64,
    /// ∂AX/∂Rθ.
    pub rt: f64,
}

/// Amplification model (`IDAMP`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AmplificationModel {
    /// Envelope e^N method (`DAMPL`, IDAMP = 0, the default).
    #[default]
    Envelope,
    /// Modified envelope with the Orr–Sommerfeld maximum for separated profiles
    /// (`DAMPL2`, IDAMP = 1).
    Modified,
}

/// XFOIL's truncated ln(10) in DAMPL/DAMPL2. Not `LN_10`: the literal must match for bit parity.
#[allow(clippy::approx_constant)]
const LN10_DAMPL: f64 = 2.3025851;

/// Half-width of the log10(Rθ/Rθcrit) ramp (`DGR`).
const DGR: f64 = 0.08;

/// Shared Falkner–Skan part of DAMPL/DAMPL2. `extra_af` adds DAMPL2's `0.1*EXP(-20*HMI)` term.
fn envelope(hk: f64, th: f64, rt: f64, extra_af: bool) -> (Amplification, f64, f64) {
    let hmi = 1.0 / (hk - 1.0);
    let hmi_hk = -powi(hmi, 2);

    // log10(critical Rθ) – H correlation for Falkner–Skan profiles
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;
    let bb = (14.0 * hmi - 9.24).tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;
    let grcrit = aa + 0.7 * (bb + 1.0);
    let grc_hk = aa_hk + 0.7 * bb_hk;

    let gr = rt.log10();
    let gr_rt = 1.0 / (LN10_DAMPL * rt);

    if gr < grcrit - DGR {
        // no amplification below the critical Rθ
        return (Amplification::default(), gr, gr_rt);
    }

    // steep cubic ramp turning AX on smoothly around Rθcrit
    let rnorm = (gr - (grcrit - DGR)) / (2.0 * DGR);
    let rn_hk = -grc_hk / (2.0 * DGR);
    let rn_rt = gr_rt / (2.0 * DGR);
    let (rfac, rfac_hk, rfac_rt) = if rnorm >= 1.0 {
        (1.0, 0.0, 0.0)
    } else {
        let rfac = 3.0 * powi(rnorm, 2) - 2.0 * powi(rnorm, 3);
        let rfac_rn = 6.0 * rnorm - 6.0 * powi(rnorm, 2);
        (rfac, rfac_rn * rn_hk, rfac_rn * rn_rt)
    };

    // envelope slope dN/dRθ for Falkner–Skan profiles
    let arg = 3.87 * hmi - 2.52;
    let arg_hk = 3.87 * hmi_hk;
    let ex = (-powi(arg, 2)).exp();
    let ex_hk = ex * (-2.0 * arg * arg_hk);
    let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    let dadr_hk = 0.028 - 0.0345 * ex_hk;

    // m(H) correlation (1 March 91): conversion d/dRθ -> d/dx
    let (af, af_hmi) = if extra_af {
        let brg = -20.0 * hmi;
        (
            -0.05 + 2.7 * hmi - 5.5 * powi(hmi, 2) + 3.0 * powi(hmi, 3) + 0.1 * brg.exp(),
            2.7 - 11.0 * hmi + 9.0 * powi(hmi, 2) - 2.0 * brg.exp(),
        )
    } else {
        (
            -0.05 + 2.7 * hmi - 5.5 * powi(hmi, 2) + 3.0 * powi(hmi, 3),
            2.7 - 11.0 * hmi + 9.0 * powi(hmi, 2),
        )
    };
    let af_hk = af_hmi * hmi_hk;

    let ax = (af * dadr / th) * rfac;
    let a = Amplification {
        ax,
        hk: (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk,
        th: -ax / th,
        rt: (af * dadr / th) * rfac_rt,
    };
    (a, gr, gr_rt)
}

/// Envelope amplification rate. Port of XFOIL `DAMPL`.
pub fn dampl(hk: f64, th: f64, rt: f64) -> Amplification {
    envelope(hk, th, rt, false).0
}

/// Modified envelope amplification rate. Port of XFOIL `DAMPL2`.
pub fn dampl2(hk: f64, th: f64, rt: f64) -> Amplification {
    const HK1: f64 = 3.5;
    const HK2: f64 = 4.0;
    let (ax1, gr, gr_rt) = envelope(hk, th, rt, true);
    if hk < HK1 {
        return ax1;
    }

    // non-envelope maximum-amplification correction for separated profiles
    let hnorm = (hk - HK1) / (HK2 - HK1);
    let hn_hk = 1.0 / (HK2 - HK1);
    let (hfac, hf_hk) = if hnorm >= 1.0 {
        (1.0, 0.0)
    } else {
        (
            3.0 * powi(hnorm, 2) - 2.0 * powi(hnorm, 3),
            (6.0 * hnorm - 6.0 * powi(hnorm, 2)) * hn_hk,
        )
    };

    let gr0 = 0.30 + 0.35 * (-0.15 * (hk - 5.0)).exp();
    let gr0_hk = -0.35 * (-0.15 * (hk - 5.0)).exp() * 0.15;
    let tnr = (1.2 * (gr - gr0)).tanh();
    let tnr_rt = (1.0 - powi(tnr, 2)) * 1.2 * gr_rt;
    let tnr_hk = -(1.0 - powi(tnr, 2)) * 1.2 * gr0_hk;

    let mut ax2 = Amplification {
        ax: (0.086 * tnr - 0.25 / (hk - 1.0).powf(1.5)) / th,
        hk: (0.086 * tnr_hk + 1.5 * 0.25 / (hk - 1.0).powf(2.5)) / th,
        rt: (0.086 * tnr_rt) / th,
        th: 0.0,
    };
    ax2.th = -ax2.ax / th;
    if ax2.ax < 0.0 {
        ax2 = Amplification::default();
    }

    // blend the two rates
    Amplification {
        ax: hfac * ax2.ax + (1.0 - hfac) * ax1.ax,
        hk: hfac * ax2.hk + (1.0 - hfac) * ax1.hk + hf_hk * (ax2.ax - ax1.ax),
        rt: hfac * ax2.rt + (1.0 - hfac) * ax1.rt,
        th: hfac * ax2.th + (1.0 - hfac) * ax1.th,
    }
}

/// Amplification rate for `model`.
pub fn rate(model: AmplificationModel, hk: f64, th: f64, rt: f64) -> Amplification {
    match model {
        AmplificationModel::Envelope => dampl(hk, th, rt),
        AmplificationModel::Modified => dampl2(hk, th, rt),
    }
}

/// Inputs of one station for [`axset`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AmpStation {
    /// Kinematic shape parameter.
    pub hk: f64,
    /// Momentum thickness.
    pub th: f64,
    /// Momentum-thickness Reynolds number.
    pub rt: f64,
    /// Amplification exponent N.
    pub ampl: f64,
}

/// Averaged amplification rate over an interval, with partials w.r.t. both stations.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct IntervalRate {
    /// Averaged dN/dx.
    pub ax: f64,
    /// Partials w.r.t. station 1: (HK, θ, Rθ, N).
    pub d1: [f64; 4],
    /// Partials w.r.t. station 2: (HK, θ, Rθ, N).
    pub d2: [f64; 4],
}

/// RMS-averaged amplification over interval 1–2 plus a small term that keeps dN/dx > 0
/// near N = Ncrit. Port of XFOIL `AXSET` (2nd-order version).
pub fn axset(
    s1: AmpStation,
    s2: AmpStation,
    acrit: f64,
    model: AmplificationModel,
) -> IntervalRate {
    let a1 = rate(model, s1.hk, s1.th, s1.rt);
    let a2 = rate(model, s2.hk, s2.th, s2.rt);

    // rms average
    let axsq = 0.5 * (powi(a1.ax, 2) + powi(a2.ax, 2));
    let (axa, axa_ax1, axa_ax2) = if axsq <= 0.0 {
        (0.0, 0.0, 0.0)
    } else {
        let axa = axsq.sqrt();
        (axa, 0.5 * a1.ax / axa, 0.5 * a2.ax / axa)
    };

    // small additional term to ensure dN/dx > 0 near N = Ncrit
    let arg = (20.0 * (acrit - 0.5 * (s1.ampl + s2.ampl))).min(20.0);
    let (exn, exn_a1, exn_a2) = if arg <= 0.0 {
        (1.0, 0.0, 0.0)
    } else {
        let exn = (-arg).exp();
        (exn, 20.0 * 0.5 * exn, 20.0 * 0.5 * exn)
    };
    let tt = s1.th + s2.th;
    let dax = exn * 0.002 / tt;
    let dax_a1 = exn_a1 * 0.002 / tt;
    let dax_a2 = exn_a2 * 0.002 / tt;
    let dax_t1 = -dax / tt;
    let dax_t2 = -dax / tt;

    IntervalRate {
        ax: axa + dax,
        d1: [
            axa_ax1 * a1.hk,
            axa_ax1 * a1.th + dax_t1,
            axa_ax1 * a1.rt,
            dax_a1,
        ],
        d2: [
            axa_ax2 * a2.hk,
            axa_ax2 * a2.th + dax_t2,
            axa_ax2 * a2.rt,
            dax_a2,
        ],
    }
}
