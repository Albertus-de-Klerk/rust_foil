//! Discretised BL equations of one interval. Port of XFOIL/QFoil `BLDIF`, `TRCHEK2`,
//! `TRDIF`, `BLSYS`, `TESYS`.

use super::station::{Kernel, LocalSystem, Regime};
use super::transition::{AmpStation, axset};
use crate::fortran::powi;

/// Interval type of BLDIF (`ITYP` 0–3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interval {
    /// Similarity (stagnation) station (ITYP = 0).
    Similarity,
    /// Laminar interval (ITYP = 1).
    Laminar,
    /// Turbulent interval (ITYP = 2).
    Turbulent,
    /// Wake interval (ITYP = 3).
    Wake,
}

impl Kernel {
    /// Newton-system coefficients and residuals of the current interval.
    /// Port of XFOIL/QFoil `BLDIF` (QFoil changes D1–D3 in the shear-lag equation).
    pub fn bldif(&mut self, ityp: Interval) {
        let (s1, s2) = (self.s1, self.s2);
        let bl = &self.bl;
        let mid = self.mid;

        let (xlog, ulog, tlog, hlog, ddlog) = if ityp == Interval::Similarity {
            // prescribed similarity logarithmic differences
            (1.0, self.p.bule, 0.5 * (1.0 - self.p.bule), 0.0, 0.0)
        } else {
            (
                (s2.x / s1.x).ln(),
                (s2.u / s1.u).ln(),
                (s2.t / s1.t).ln(),
                (s2.hs / s1.hs).ln(),
                1.0,
            )
        };

        self.sys.clear();
        let sys = &mut self.sys;

        // local upwinding, based on the local change in log(Hk-1)
        let hupwt = 1.0;
        let (hdcon, hd_hk1) = if ityp == Interval::Wake {
            (hupwt / powi(s2.hk, 2), 0.0)
        } else {
            (5.0 * hupwt / powi(s2.hk, 2), 0.0)
        };
        let hd_hk2 = -hdcon * 2.0 / s2.hk;
        let arg = ((s2.hk - 1.0) / (s1.hk - 1.0)).abs();
        let hl = arg.ln();
        let hl_hk1 = -1.0 / (s1.hk - 1.0);
        let hl_hk2 = 1.0 / (s2.hk - 1.0);
        // UPW = 0.5 trapezoidal ... 1.0 backward Euler
        let hlsq = powi(hl, 2).min(15.0);
        let ehh = (-hlsq * hdcon).exp();
        let upw = 1.0 - 0.5 * ehh;
        let upw_hl = ehh * hl * hdcon;
        let upw_hd = 0.5 * ehh * hlsq;
        let upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1;
        let upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;
        let upw_u1 = upw_hk1 * s1.hk_u;
        let upw_t1 = upw_hk1 * s1.hk_t;
        let upw_d1 = upw_hk1 * s1.hk_d;
        let upw_u2 = upw_hk2 * s2.hk_u;
        let upw_t2 = upw_hk2 * s2.hk_t;
        let upw_d2 = upw_hk2 * s2.hk_d;
        let upw_ms = upw_hk1 * s1.hk_ms + upw_hk2 * s2.hk_ms;

        match ityp {
            Interval::Similarity => {
                // LE point: zero amplification
                sys.vs2[0][0] = 1.0;
                sys.vsr[0] = 0.0;
                sys.vsrez[0] = -s2.ampl;
            }
            Interval::Laminar => {
                // amplification equation, averaged rate over X1..X2
                let a = axset(
                    AmpStation {
                        hk: s1.hk,
                        th: s1.t,
                        rt: s1.rt,
                        ampl: s1.ampl,
                    },
                    AmpStation {
                        hk: s2.hk,
                        th: s2.t,
                        rt: s2.rt,
                        ampl: s2.ampl,
                    },
                    self.p.amcrit,
                    self.model,
                );
                let (ax, [ax_hk1, ax_t1, ax_rt1, ax_a1], [ax_hk2, ax_t2, ax_rt2, ax_a2]) =
                    (a.ax, a.d1, a.d2);
                let rezc = s2.ampl - s1.ampl - ax * (s2.x - s1.x);
                let z_ax = -(s2.x - s1.x);
                sys.vs1[0][0] = z_ax * ax_a1 - 1.0;
                sys.vs1[0][1] = z_ax * (ax_hk1 * s1.hk_t + ax_t1 + ax_rt1 * s1.rt_t);
                sys.vs1[0][2] = z_ax * (ax_hk1 * s1.hk_d);
                sys.vs1[0][3] = z_ax * (ax_hk1 * s1.hk_u + ax_rt1 * s1.rt_u);
                sys.vs1[0][4] = ax;
                sys.vs2[0][0] = z_ax * ax_a2 + 1.0;
                sys.vs2[0][1] = z_ax * (ax_hk2 * s2.hk_t + ax_t2 + ax_rt2 * s2.rt_t);
                sys.vs2[0][2] = z_ax * (ax_hk2 * s2.hk_d);
                sys.vs2[0][3] = z_ax * (ax_hk2 * s2.hk_u + ax_rt2 * s2.rt_u);
                sys.vs2[0][4] = -ax;
                sys.vsm[0] = z_ax
                    * (ax_hk1 * s1.hk_ms
                        + ax_rt1 * s1.rt_ms
                        + ax_hk2 * s2.hk_ms
                        + ax_rt2 * s2.rt_ms);
                sys.vsr[0] = z_ax * (ax_rt1 * s1.rt_re + ax_rt2 * s2.rt_re);
                sys.vsx[0] = 0.0;
                sys.vsrez[0] = -rezc;
            }
            Interval::Turbulent | Interval::Wake => {
                // shear-lag equation
                let sa = (1.0 - upw) * s1.s + upw * s2.s;
                let cqa = (1.0 - upw) * s1.cq + upw * s2.cq;
                let cfa = (1.0 - upw) * s1.cf + upw * s2.cf;
                let hka = (1.0 - upw) * s1.hk + upw * s2.hk;
                let usa = 0.5 * (s1.us + s2.us);
                let rta = 0.5 * (s1.rt + s2.rt);
                let dea = 0.5 * (s1.de + s2.de);
                let da = 0.5 * (s1.d + s2.d);

                // QFoil D3: wake dissipation length DLCON/4 (XFOIL: DLCON)
                let ald = if ityp == Interval::Wake {
                    bl.dlcon / 4.0
                } else {
                    1.0
                };

                // equilibrium 1/Ue dUe/dx (12 Oct 94)
                let (hkc, hkc_hka, hkc_rta) = if ityp == Interval::Turbulent {
                    let gcc = bl.gccon;
                    let hkc = hka - 1.0 - gcc / rta;
                    if hkc < 0.01 {
                        (0.01, 0.0, 0.0)
                    } else {
                        (hkc, 1.0, gcc / powi(rta, 2))
                    }
                } else {
                    (hka - 1.0, 1.0, 0.0)
                };
                let hr = hkc / (bl.gacon * ald * hka);
                let hr_hka = hkc_hka / (bl.gacon * ald * hka) - hr / hka;
                let hr_rta = hkc_rta / (bl.gacon * ald * hka);

                let uq = (0.5 * cfa - powi(hr, 2)) / (bl.gbcon * da);
                let uq_hka = -2.0 * hr * hr_hka / (bl.gbcon * da);
                let uq_rta = -2.0 * hr * hr_rta / (bl.gbcon * da);
                let uq_cfa = 0.5 / (bl.gbcon * da);
                let uq_da = -uq / da;
                let uq_upw = uq_cfa * (s2.cf - s1.cf) + uq_hka * (s2.hk - s1.hk);

                let mut uq_t1 =
                    (1.0 - upw) * (uq_cfa * s1.cf_t + uq_hka * s1.hk_t) + uq_upw * upw_t1;
                let mut uq_d1 =
                    (1.0 - upw) * (uq_cfa * s1.cf_d + uq_hka * s1.hk_d) + uq_upw * upw_d1;
                let mut uq_u1 =
                    (1.0 - upw) * (uq_cfa * s1.cf_u + uq_hka * s1.hk_u) + uq_upw * upw_u1;
                let mut uq_t2 = upw * (uq_cfa * s2.cf_t + uq_hka * s2.hk_t) + uq_upw * upw_t2;
                let mut uq_d2 = upw * (uq_cfa * s2.cf_d + uq_hka * s2.hk_d) + uq_upw * upw_d2;
                let mut uq_u2 = upw * (uq_cfa * s2.cf_u + uq_hka * s2.hk_u) + uq_upw * upw_u2;
                let mut uq_ms = (1.0 - upw) * (uq_cfa * s1.cf_ms + uq_hka * s1.hk_ms)
                    + uq_upw * upw_ms
                    + upw * (uq_cfa * s2.cf_ms + uq_hka * s2.hk_ms);
                let mut uq_re = (1.0 - upw) * uq_cfa * s1.cf_re + upw * uq_cfa * s2.cf_re;
                uq_t1 += 0.5 * uq_rta * s1.rt_t;
                uq_d1 += 0.5 * uq_da;
                uq_u1 += 0.5 * uq_rta * s1.rt_u;
                uq_t2 += 0.5 * uq_rta * s2.rt_t;
                uq_d2 += 0.5 * uq_da;
                uq_u2 += 0.5 * uq_rta * s2.rt_u;
                uq_ms = uq_ms + 0.5 * uq_rta * s1.rt_ms + 0.5 * uq_rta * s2.rt_ms;
                uq_re = uq_re + 0.5 * uq_rta * s1.rt_re + 0.5 * uq_rta * s2.rt_re;
                // UQ's T/D/U/MS/RE partials are formed but not used by the shear-lag
                // row, exactly as in XFOIL (Z_HKA, Z_CFA and Z_DA carry them).
                let _ = (uq_t1, uq_d1, uq_u1, uq_t2, uq_d2, uq_u2, uq_ms, uq_re);

                // QFoil D1: shear-lag coefficient Kc(Hk) = 4.65 - 0.95 tanh(0.5 (Hk - 3.5))
                // (XFOIL: SCC = SCCON*1.333/(1+USA)); SCC_HKA is the new Jacobian term (D2).
                let tanh_arg = 0.5 * (hka - 3.5);
                let tanh_val = tanh_arg.tanh();
                let xkcvar = 4.65 - 0.95 * tanh_val;
                let sech2_val = 1.0 - powi(tanh_val, 2);
                let dxk_dh = -0.95 * sech2_val * 0.5;
                let scc = xkcvar * 1.333 / (1.0 + usa);
                let scc_usa = -scc / (1.0 + usa);
                let scc_hka = dxk_dh * 1.333 / (1.0 + usa);

                let slog = (s2.s / s1.s).ln();
                let dxi = s2.x - s1.x;
                let duxcon = bl.duxcon;
                let rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog
                    + dea * 2.0 * (uq * dxi - ulog) * duxcon;

                let z_cfa = dea * 2.0 * uq_cfa * dxi * duxcon;
                // QFoil D2: + SCC_HKA*(CQA - SA*ALD)*DXI
                let z_hka = dea * 2.0 * uq_hka * dxi * duxcon + scc_hka * (cqa - sa * ald) * dxi;
                let z_da = dea * 2.0 * uq_da * dxi * duxcon;
                let z_sl = -dea * 2.0;
                let z_ul = -dea * 2.0 * duxcon;
                let z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq * duxcon;
                let z_usa = scc_usa * (cqa - sa * ald) * dxi;
                let z_cqa = scc * dxi;
                let z_sa = -scc * dxi * ald;
                let z_dea = 2.0 * ((uq * dxi - ulog) * duxcon - slog);

                let z_upw = z_cqa * (s2.cq - s1.cq)
                    + z_sa * (s2.s - s1.s)
                    + z_cfa * (s2.cf - s1.cf)
                    + z_hka * (s2.hk - s1.hk);
                let z_de1 = 0.5 * z_dea;
                let z_de2 = 0.5 * z_dea;
                let z_us1 = 0.5 * z_usa;
                let z_us2 = 0.5 * z_usa;
                let z_d1 = 0.5 * z_da;
                let z_d2 = 0.5 * z_da;
                let z_u1 = -z_ul / s1.u;
                let z_u2 = z_ul / s2.u;
                let z_x1 = -z_dxi;
                let z_x2 = z_dxi;
                let z_s1 = (1.0 - upw) * z_sa - z_sl / s1.s;
                let z_s2 = upw * z_sa + z_sl / s2.s;
                let z_cq1 = (1.0 - upw) * z_cqa;
                let z_cq2 = upw * z_cqa;
                let z_cf1 = (1.0 - upw) * z_cfa;
                let z_cf2 = upw * z_cfa;
                let z_hk1 = (1.0 - upw) * z_hka;
                let z_hk2 = upw * z_hka;

                sys.vs1[0][0] = z_s1;
                sys.vs1[0][1] = z_upw * upw_t1 + z_de1 * s1.de_t + z_us1 * s1.us_t;
                sys.vs1[0][2] = z_d1 + z_upw * upw_d1 + z_de1 * s1.de_d + z_us1 * s1.us_d;
                sys.vs1[0][3] = z_u1 + z_upw * upw_u1 + z_de1 * s1.de_u + z_us1 * s1.us_u;
                sys.vs1[0][4] = z_x1;
                sys.vs2[0][0] = z_s2;
                sys.vs2[0][1] = z_upw * upw_t2 + z_de2 * s2.de_t + z_us2 * s2.us_t;
                sys.vs2[0][2] = z_d2 + z_upw * upw_d2 + z_de2 * s2.de_d + z_us2 * s2.us_d;
                sys.vs2[0][3] = z_u2 + z_upw * upw_u2 + z_de2 * s2.de_u + z_us2 * s2.us_u;
                sys.vs2[0][4] = z_x2;
                sys.vsm[0] = z_upw * upw_ms
                    + z_de1 * s1.de_ms
                    + z_us1 * s1.us_ms
                    + z_de2 * s2.de_ms
                    + z_us2 * s2.us_ms;

                sys.vs1[0][1] = sys.vs1[0][1] + z_cq1 * s1.cq_t + z_cf1 * s1.cf_t + z_hk1 * s1.hk_t;
                sys.vs1[0][2] = sys.vs1[0][2] + z_cq1 * s1.cq_d + z_cf1 * s1.cf_d + z_hk1 * s1.hk_d;
                sys.vs1[0][3] = sys.vs1[0][3] + z_cq1 * s1.cq_u + z_cf1 * s1.cf_u + z_hk1 * s1.hk_u;
                sys.vs2[0][1] = sys.vs2[0][1] + z_cq2 * s2.cq_t + z_cf2 * s2.cf_t + z_hk2 * s2.hk_t;
                sys.vs2[0][2] = sys.vs2[0][2] + z_cq2 * s2.cq_d + z_cf2 * s2.cf_d + z_hk2 * s2.hk_d;
                sys.vs2[0][3] = sys.vs2[0][3] + z_cq2 * s2.cq_u + z_cf2 * s2.cf_u + z_hk2 * s2.hk_u;
                sys.vsm[0] = sys.vsm[0]
                    + z_cq1 * s1.cq_ms
                    + z_cf1 * s1.cf_ms
                    + z_hk1 * s1.hk_ms
                    + z_cq2 * s2.cq_ms
                    + z_cf2 * s2.cf_ms
                    + z_hk2 * s2.hk_ms;
                sys.vsr[0] =
                    z_cq1 * s1.cq_re + z_cf1 * s1.cf_re + z_cq2 * s2.cq_re + z_cf2 * s2.cf_re;
                sys.vsx[0] = 0.0;
                sys.vsrez[0] = -rezc;
            }
        }

        // ---- momentum equation
        let ha = 0.5 * (s1.h + s2.h);
        let ma = 0.5 * (s1.m + s2.m);
        let xa = 0.5 * (s1.x + s2.x);
        let ta = 0.5 * (s1.t + s2.t);
        let hwa = 0.5 * (s1.dw / s1.t + s2.dw / s2.t);

        // Cf term, using the central value CFM for better drag accuracy
        let cfx = 0.50 * mid.cfm * xa / ta + 0.25 * (s1.cf * s1.x / s1.t + s2.cf * s2.x / s2.t);
        let cfx_xa = 0.50 * mid.cfm / ta;
        let cfx_ta = -0.50 * mid.cfm * xa / powi(ta, 2);
        let cfx_x1 = 0.25 * s1.cf / s1.t + cfx_xa * 0.5;
        let cfx_x2 = 0.25 * s2.cf / s2.t + cfx_xa * 0.5;
        let cfx_t1 = -0.25 * s1.cf * s1.x / powi(s1.t, 2) + cfx_ta * 0.5;
        let cfx_t2 = -0.25 * s2.cf * s2.x / powi(s2.t, 2) + cfx_ta * 0.5;
        let cfx_cf1 = 0.25 * s1.x / s1.t;
        let cfx_cf2 = 0.25 * s2.x / s2.t;
        let cfx_cfm = 0.50 * xa / ta;

        let btmp = ha + 2.0 - ma + hwa;
        let rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;
        let z_cfx = -xlog * 0.5;
        let z_ha = ulog;
        let z_hwa = ulog;
        let z_ma = -ulog;
        let z_xl = -ddlog * 0.5 * cfx;
        let z_ul = ddlog * btmp;
        let z_tl = ddlog;
        let z_cfm = z_cfx * cfx_cfm;
        let z_cf1 = z_cfx * cfx_cf1;
        let z_cf2 = z_cfx * cfx_cf2;
        let z_t1 = -z_tl / s1.t + z_cfx * cfx_t1 + z_hwa * 0.5 * (-s1.dw / powi(s1.t, 2));
        let z_t2 = z_tl / s2.t + z_cfx * cfx_t2 + z_hwa * 0.5 * (-s2.dw / powi(s2.t, 2));
        let z_x1 = -z_xl / s1.x + z_cfx * cfx_x1;
        let z_x2 = z_xl / s2.x + z_cfx * cfx_x2;
        let z_u1 = -z_ul / s1.u;
        let z_u2 = z_ul / s2.u;

        sys.vs1[1][1] = 0.5 * z_ha * s1.h_t + z_cfm * mid.cfm_t1 + z_cf1 * s1.cf_t + z_t1;
        sys.vs1[1][2] = 0.5 * z_ha * s1.h_d + z_cfm * mid.cfm_d1 + z_cf1 * s1.cf_d;
        sys.vs1[1][3] = 0.5 * z_ma * s1.m_u + z_cfm * mid.cfm_u1 + z_cf1 * s1.cf_u + z_u1;
        sys.vs1[1][4] = z_x1;
        sys.vs2[1][1] = 0.5 * z_ha * s2.h_t + z_cfm * mid.cfm_t2 + z_cf2 * s2.cf_t + z_t2;
        sys.vs2[1][2] = 0.5 * z_ha * s2.h_d + z_cfm * mid.cfm_d2 + z_cf2 * s2.cf_d;
        sys.vs2[1][3] = 0.5 * z_ma * s2.m_u + z_cfm * mid.cfm_u2 + z_cf2 * s2.cf_u + z_u2;
        sys.vs2[1][4] = z_x2;
        sys.vsm[1] = 0.5 * z_ma * s1.m_ms
            + z_cfm * mid.cfm_ms
            + z_cf1 * s1.cf_ms
            + 0.5 * z_ma * s2.m_ms
            + z_cf2 * s2.cf_ms;
        sys.vsr[1] = z_cfm * mid.cfm_re + z_cf1 * s1.cf_re + z_cf2 * s2.cf_re;
        sys.vsx[1] = 0.0;
        sys.vsrez[1] = -rezt;

        // ---- shape parameter equation
        let xot1 = s1.x / s1.t;
        let xot2 = s2.x / s2.t;
        let ha = 0.5 * (s1.h + s2.h);
        let hsa = 0.5 * (s1.hs + s2.hs);
        let hca = 0.5 * (s1.hc + s2.hc);
        let hwa = 0.5 * (s1.dw / s1.t + s2.dw / s2.t);

        let dix = (1.0 - upw) * s1.di * xot1 + upw * s2.di * xot2;
        let cfx = (1.0 - upw) * s1.cf * xot1 + upw * s2.cf * xot2;
        let dix_upw = s2.di * xot2 - s1.di * xot1;
        let cfx_upw = s2.cf * xot2 - s1.cf * xot1;

        let btmp = 2.0 * hca / hsa + 1.0 - ha - hwa;
        let rezh = hlog + btmp * ulog + xlog * (0.5 * cfx - dix);
        let z_cfx = xlog * 0.5;
        let z_dix = -xlog;
        let z_hca = 2.0 * ulog / hsa;
        let z_ha = -ulog;
        let z_hwa = -ulog;
        let z_xl = ddlog * (0.5 * cfx - dix);
        let z_ul = ddlog * btmp;
        let z_hl = ddlog;
        let z_upw = z_cfx * cfx_upw + z_dix * dix_upw;
        let z_hs1 = -hca * ulog / powi(hsa, 2) - z_hl / s1.hs;
        let z_hs2 = -hca * ulog / powi(hsa, 2) + z_hl / s2.hs;
        let z_cf1 = (1.0 - upw) * z_cfx * xot1;
        let z_cf2 = upw * z_cfx * xot2;
        let z_di1 = (1.0 - upw) * z_dix * xot1;
        let z_di2 = upw * z_dix * xot2;
        let mut z_t1 = (1.0 - upw) * (z_cfx * s1.cf + z_dix * s1.di) * (-xot1 / s1.t);
        let mut z_t2 = upw * (z_cfx * s2.cf + z_dix * s2.di) * (-xot2 / s2.t);
        let z_x1 = (1.0 - upw) * (z_cfx * s1.cf + z_dix * s1.di) / s1.t - z_xl / s1.x;
        let z_x2 = upw * (z_cfx * s2.cf + z_dix * s2.di) / s2.t + z_xl / s2.x;
        let z_u1 = -z_ul / s1.u;
        let z_u2 = z_ul / s2.u;
        z_t1 += z_hwa * 0.5 * (-s1.dw / powi(s1.t, 2));
        z_t2 += z_hwa * 0.5 * (-s2.dw / powi(s2.t, 2));

        sys.vs1[2][0] = z_di1 * s1.di_s;
        sys.vs1[2][1] = z_hs1 * s1.hs_t + z_cf1 * s1.cf_t + z_di1 * s1.di_t + z_t1;
        sys.vs1[2][2] = z_hs1 * s1.hs_d + z_cf1 * s1.cf_d + z_di1 * s1.di_d;
        sys.vs1[2][3] = z_hs1 * s1.hs_u + z_cf1 * s1.cf_u + z_di1 * s1.di_u + z_u1;
        sys.vs1[2][4] = z_x1;
        sys.vs2[2][0] = z_di2 * s2.di_s;
        sys.vs2[2][1] = z_hs2 * s2.hs_t + z_cf2 * s2.cf_t + z_di2 * s2.di_t + z_t2;
        sys.vs2[2][2] = z_hs2 * s2.hs_d + z_cf2 * s2.cf_d + z_di2 * s2.di_d;
        sys.vs2[2][3] = z_hs2 * s2.hs_u + z_cf2 * s2.cf_u + z_di2 * s2.di_u + z_u2;
        sys.vs2[2][4] = z_x2;
        sys.vsm[2] = z_hs1 * s1.hs_ms
            + z_cf1 * s1.cf_ms
            + z_di1 * s1.di_ms
            + z_hs2 * s2.hs_ms
            + z_cf2 * s2.cf_ms
            + z_di2 * s2.di_ms;
        sys.vsr[2] = z_hs1 * s1.hs_re
            + z_cf1 * s1.cf_re
            + z_di1 * s1.di_re
            + z_hs2 * s2.hs_re
            + z_cf2 * s2.cf_re
            + z_di2 * s2.di_re;

        sys.vs1[2][1] = sys.vs1[2][1] + 0.5 * (z_hca * s1.hc_t + z_ha * s1.h_t) + z_upw * upw_t1;
        sys.vs1[2][2] = sys.vs1[2][2] + 0.5 * (z_hca * s1.hc_d + z_ha * s1.h_d) + z_upw * upw_d1;
        sys.vs1[2][3] = sys.vs1[2][3] + 0.5 * (z_hca * s1.hc_u) + z_upw * upw_u1;
        sys.vs2[2][1] = sys.vs2[2][1] + 0.5 * (z_hca * s2.hc_t + z_ha * s2.h_t) + z_upw * upw_t2;
        sys.vs2[2][2] = sys.vs2[2][2] + 0.5 * (z_hca * s2.hc_d + z_ha * s2.h_d) + z_upw * upw_d2;
        sys.vs2[2][3] = sys.vs2[2][3] + 0.5 * (z_hca * s2.hc_u) + z_upw * upw_u2;
        sys.vsm[2] =
            sys.vsm[2] + 0.5 * (z_hca * s1.hc_ms) + z_upw * upw_ms + 0.5 * (z_hca * s2.hc_ms);
        sys.vsx[2] = 0.0;
        sys.vsrez[2] = -rezh;
    }

    /// Solves the implicit amplification equation over X1..X2 and, if transition occurs,
    /// sets XT and its sensitivities. Port of XFOIL `TRCHEK2` (TRCHEK calls it).
    pub fn trchek(&mut self) {
        const DAEPS: f64 = 5.0e-5;
        let amcrit = self.p.amcrit;
        let xiforc = self.p.xiforc;
        let s1 = self.s1;
        let c2sav = self.s2;

        let amp = |s: &super::station::Station, a: f64| AmpStation {
            hk: s.hk,
            th: s.t,
            rt: s.rt,
            ampl: a,
        };
        let a0 = axset(
            amp(&s1, s1.ampl),
            amp(&c2sav, c2sav.ampl),
            amcrit,
            self.model,
        );
        let mut ampl2 = s1.ampl + a0.ax * (c2sav.x - s1.x);

        let (x1, x2) = (s1.x, c2sav.x);
        let (t1, t2, d1, d2, u1, u2) = (s1.t, c2sav.t, s1.d, c2sav.d, s1.u, c2sav.u);

        // quantities of the final iterate, needed after the loop
        let mut it = TrIter::default();

        for _ in 0..30 {
            // weighting factors defining the "T" point from 1 and 2
            let (amplt, amplt_a2, sfa, sfa_a1, sfa_a2) = if ampl2 <= amcrit {
                (ampl2, 1.0, 1.0, 0.0, 0.0)
            } else {
                let sfa = (amcrit - s1.ampl) / (ampl2 - s1.ampl);
                (
                    amcrit,
                    0.0,
                    sfa,
                    (sfa - 1.0) / (ampl2 - s1.ampl),
                    (-sfa) / (ampl2 - s1.ampl),
                )
            };
            let (sfx, sfx_x1, sfx_x2, sfx_xf) = if xiforc < x2 {
                let sfx = (xiforc - x1) / (x2 - x1);
                (
                    sfx,
                    (sfx - 1.0) / (x2 - x1),
                    (-sfx) / (x2 - x1),
                    1.0 / (x2 - x1),
                )
            } else {
                (1.0, 0.0, 0.0, 0.0)
            };
            // weighting factor from free or forced transition
            let (wf2, wf2_a1, wf2_a2, wf2_x1, wf2_x2, wf2_xf) = if sfa < sfx {
                (sfa, sfa_a1, sfa_a2, 0.0, 0.0, 0.0)
            } else {
                (sfx, 0.0, 0.0, sfx_x1, sfx_x2, sfx_xf)
            };
            let wf1 = 1.0 - wf2;
            let (wf1_a1, wf1_a2, wf1_x1, wf1_x2, wf1_xf) =
                (-wf2_a1, -wf2_a2, -wf2_x1, -wf2_x2, -wf2_xf);

            // interpolate BL variables to XT
            let xt = x1 * wf1 + x2 * wf2;
            let tt = t1 * wf1 + t2 * wf2;
            let dt = d1 * wf1 + d2 * wf2;
            let ut = u1 * wf1 + u2 * wf2;
            let xt_a2 = x1 * wf1_a2 + x2 * wf2_a2;
            let tt_a2 = t1 * wf1_a2 + t2 * wf2_a2;
            let dt_a2 = d1 * wf1_a2 + d2 * wf2_a2;
            let ut_a2 = u1 * wf1_a2 + u2 * wf2_a2;

            // laminar secondary "T" variables: temporarily clobber station 2 for BLKIN
            self.s2.x = xt;
            self.s2.t = tt;
            self.s2.d = dt;
            self.s2.u = ut;
            self.blkin();
            let st = self.s2;
            // restore station 2 except AMPL2
            self.s2 = c2sav;
            self.s2.ampl = ampl2;

            // amplification rate over X1..XT
            let a = axset(amp(&s1, s1.ampl), amp(&st, amplt), amcrit, self.model);
            it = TrIter {
                wf1,
                wf2,
                wf1_a1,
                wf2_a1,
                wf1_x1,
                wf2_x1,
                wf1_x2,
                wf2_x2,
                wf1_xf,
                wf2_xf,
                xt,
                xt_a2,
                tt_a2,
                dt_a2,
                ut_a2,
                amplt_a2,
                st,
                a,
            };
            self.xt.xt = xt;
            if a.ax <= 0.0 {
                break; // no amplification here
            }
            let [ax_hkt, ax_tt, ax_rtt, ax_at] = a.d2;
            let ax_a2 = (ax_hkt * st.hk_t + ax_tt + ax_rtt * st.rt_t) * tt_a2
                + (ax_hkt * st.hk_d) * dt_a2
                + (ax_hkt * st.hk_u + ax_rtt * st.rt_u) * ut_a2
                + ax_at * amplt_a2;

            // residual of the implicit AMPL2 definition
            let res = ampl2 - s1.ampl - a.ax * (x2 - x1);
            let res_a2 = 1.0 - ax_a2 * (x2 - x1);
            let da2 = -res / res_a2;
            let mut rlx = 1.0;
            let dxt = xt_a2 * da2;
            if rlx * (dxt / (x2 - x1)).abs() > 0.05 {
                rlx = 0.05 * ((x2 - x1) / dxt).abs();
            }
            if rlx * da2.abs() > 1.0 {
                rlx = 1.0 * (1.0 / da2).abs();
            }
            if da2.abs() < DAEPS {
                break;
            }
            if (ampl2 > amcrit && ampl2 + rlx * da2 < amcrit)
                || (ampl2 < amcrit && ampl2 + rlx * da2 > amcrit)
            {
                // do not step across AMCRIT
                ampl2 = amcrit;
            } else {
                ampl2 += rlx * da2;
            }
            self.s2.ampl = ampl2;
        }
        self.s2.ampl = ampl2;

        // free or forced transition?
        self.trfree = ampl2 >= amcrit;
        self.trforc = xiforc > x1 && xiforc <= x2;
        self.tran = self.trforc || self.trfree;
        if !self.tran {
            return;
        }
        if self.trfree && self.trforc {
            self.trforc = xiforc < it.xt;
            self.trfree = xiforc >= it.xt;
        }
        if self.trforc {
            // forced transition: XT prescribed
            self.xt = super::station::TransitionPoint {
                xt: xiforc,
                xt_xf: 1.0,
                ..Default::default()
            };
            return;
        }

        // free transition: sensitivities of XT
        let TrIter {
            wf1,
            wf2,
            wf1_a1,
            wf2_a1,
            wf1_x1,
            wf2_x1,
            wf1_x2,
            wf2_x2,
            wf1_xf,
            wf2_xf,
            xt,
            xt_a2,
            tt_a2,
            dt_a2,
            ut_a2,
            amplt_a2,
            st,
            a,
        } = it;
        let (tt_t1, dt_d1, ut_u1) = (wf1, wf1, wf1);
        let (tt_t2, dt_d2, ut_u2) = (wf2, wf2, wf2);
        let xt_a1 = x1 * wf1_a1 + x2 * wf2_a1;
        let tt_a1 = t1 * wf1_a1 + t2 * wf2_a1;
        let dt_a1 = d1 * wf1_a1 + d2 * wf2_a1;
        let ut_a1 = u1 * wf1_a1 + u2 * wf2_a1;
        let xt_x1 = x1 * wf1_x1 + x2 * wf2_x1 + wf1;
        let tt_x1 = t1 * wf1_x1 + t2 * wf2_x1;
        let dt_x1 = d1 * wf1_x1 + d2 * wf2_x1;
        let ut_x1 = u1 * wf1_x1 + u2 * wf2_x1;
        let xt_x2 = x1 * wf1_x2 + x2 * wf2_x2 + wf2;
        let tt_x2 = t1 * wf1_x2 + t2 * wf2_x2;
        let dt_x2 = d1 * wf1_x2 + d2 * wf2_x2;
        let ut_x2 = u1 * wf1_x2 + u2 * wf2_x2;
        let tt_xf = t1 * wf1_xf + t2 * wf2_xf;
        let dt_xf = d1 * wf1_xf + d2 * wf2_xf;
        let ut_xf = u1 * wf1_xf + u2 * wf2_xf;

        // AX = AX(HK1, T1, RT1, A1, HKT, TT, RTT, AT)
        let [ax_hk1, ax_t1, ax_rt1, ax_a1] = a.d1;
        let [ax_hkt, ax_tt, ax_rtt, ax_at] = a.d2;
        let gt = ax_hkt * st.hk_t + ax_tt + ax_rtt * st.rt_t;
        let gd = ax_hkt * st.hk_d;
        let gu = ax_hkt * st.hk_u + ax_rtt * st.rt_u;
        let ax_t1 = ax_hk1 * s1.hk_t + ax_t1 + ax_rt1 * s1.rt_t + gt * tt_t1;
        let ax_d1 = ax_hk1 * s1.hk_d + gd * dt_d1;
        let ax_u1 = ax_hk1 * s1.hk_u + ax_rt1 * s1.rt_u + gu * ut_u1;
        let ax_a1 = ax_a1 + gt * tt_a1 + gd * dt_a1 + gu * ut_a1;
        let ax_x1 = gt * tt_x1 + gd * dt_x1 + gu * ut_x1;
        let ax_t2 = gt * tt_t2;
        let ax_d2 = gd * dt_d2;
        let ax_u2 = gu * ut_u2;
        let ax_a2 = ax_at * amplt_a2 + gt * tt_a2 + gd * dt_a2 + gu * ut_a2;
        let ax_x2 = gt * tt_x2 + gd * dt_x2 + gu * ut_x2;
        let ax_xf = gt * tt_xf + gd * dt_xf + gu * ut_xf;
        let ax_ms = ax_hkt * st.hk_ms + ax_rtt * st.rt_ms + ax_hk1 * s1.hk_ms + ax_rt1 * s1.rt_ms;
        let ax_re = ax_rtt * st.rt_re + ax_rt1 * s1.rt_re;

        // sensitivities of the residual RES = AMPL2 - AMPL1 - AX*(X2-X1)
        let z_ax = -(x2 - x1);
        let z_a1 = z_ax * ax_a1 - 1.0;
        let z_t1 = z_ax * ax_t1;
        let z_d1 = z_ax * ax_d1;
        let z_u1 = z_ax * ax_u1;
        let z_x1 = z_ax * ax_x1 + a.ax;
        let z_a2 = z_ax * ax_a2 + 1.0;
        let z_t2 = z_ax * ax_t2;
        let z_d2 = z_ax * ax_d2;
        let z_u2 = z_ax * ax_u2;
        let z_x2 = z_ax * ax_x2 - a.ax;
        let _z_xf = z_ax * ax_xf;
        let z_ms = z_ax * ax_ms;
        let z_re = z_ax * ax_re;

        // XT sensitivities with RES stationary under the A2 constraint
        let k = xt_a2 / z_a2;
        self.xt = super::station::TransitionPoint {
            xt,
            xt_a1: xt_a1 - k * z_a1,
            xt_t1: -k * z_t1,
            xt_d1: -k * z_d1,
            xt_u1: -k * z_u1,
            xt_x1: xt_x1 - k * z_x1,
            xt_t2: -k * z_t2,
            xt_d2: -k * z_d2,
            xt_u2: -k * z_u2,
            xt_x2: xt_x2 - k * z_x2,
            xt_ms: -k * z_ms,
            xt_re: -k * z_re,
            xt_xf: 0.0,
        };
    }

    /// Newton system of a transition interval: laminar X1..XT plus turbulent XT..X2.
    /// Port of XFOIL `TRDIF`.
    pub fn trdif(&mut self) {
        let c1sav = self.s1;
        let c2sav = self.s2;
        let (s1, s2) = (c1sav, c2sav);
        let xt = self.xt;
        let bl = self.bl.clone();

        // weighting factors for linear interpolation to the transition point
        let wf2 = (xt.xt - s1.x) / (s2.x - s1.x);
        let wf2_xt = 1.0 / (s2.x - s1.x);
        let wf2_a1 = wf2_xt * xt.xt_a1;
        let wf2_x1 = wf2_xt * xt.xt_x1 + (wf2 - 1.0) / (s2.x - s1.x);
        let wf2_x2 = wf2_xt * xt.xt_x2 - wf2 / (s2.x - s1.x);
        let wf2_t1 = wf2_xt * xt.xt_t1;
        let wf2_t2 = wf2_xt * xt.xt_t2;
        let wf2_d1 = wf2_xt * xt.xt_d1;
        let wf2_d2 = wf2_xt * xt.xt_d2;
        let wf2_u1 = wf2_xt * xt.xt_u1;
        let wf2_u2 = wf2_xt * xt.xt_u2;
        let wf2_ms = wf2_xt * xt.xt_ms;
        let wf2_re = wf2_xt * xt.xt_re;
        let wf2_xf = wf2_xt * xt.xt_xf;
        let wf1 = 1.0 - wf2;
        let w1 = Sens13 {
            a1: -wf2_a1,
            x1: -wf2_x1,
            x2: -wf2_x2,
            t1: -wf2_t1,
            t2: -wf2_t2,
            d1: -wf2_d1,
            d2: -wf2_d2,
            u1: -wf2_u1,
            u2: -wf2_u2,
            ms: -wf2_ms,
            re: -wf2_re,
            xf: -wf2_xf,
        };
        let w2 = Sens13 {
            a1: wf2_a1,
            x1: wf2_x1,
            x2: wf2_x2,
            t1: wf2_t1,
            t2: wf2_t2,
            d1: wf2_d1,
            d2: wf2_d2,
            u1: wf2_u1,
            u2: wf2_u2,
            ms: wf2_ms,
            re: wf2_re,
            xf: wf2_xf,
        };

        // ---- laminar part X1..XT: interpolate primary variables to XT
        let lerp = |v1: f64, v2: f64| Sens13 {
            a1: v1 * w1.a1 + v2 * w2.a1,
            x1: v1 * w1.x1 + v2 * w2.x1,
            x2: v1 * w1.x2 + v2 * w2.x2,
            t1: v1 * w1.t1 + v2 * w2.t1,
            t2: v1 * w1.t2 + v2 * w2.t2,
            d1: v1 * w1.d1 + v2 * w2.d1,
            d2: v1 * w1.d2 + v2 * w2.d2,
            u1: v1 * w1.u1 + v2 * w2.u1,
            u2: v1 * w1.u2 + v2 * w2.u2,
            ms: v1 * w1.ms + v2 * w2.ms,
            re: v1 * w1.re + v2 * w2.re,
            xf: v1 * w1.xf + v2 * w2.xf,
        };
        let tt = s1.t * wf1 + s2.t * wf2;
        let mut ts = lerp(s1.t, s2.t);
        ts.t1 += wf1;
        ts.t2 += wf2;
        let dt = s1.d * wf1 + s2.d * wf2;
        let mut ds = lerp(s1.d, s2.d);
        ds.d1 += wf1;
        ds.d2 += wf2;
        let ut = s1.u * wf1 + s2.u * wf2;
        let mut us = lerp(s1.u, s2.u);
        us.u1 += wf1;
        us.u2 += wf2;
        let xs = Sens13 {
            a1: xt.xt_a1,
            x1: xt.xt_x1,
            x2: xt.xt_x2,
            t1: xt.xt_t1,
            t2: xt.xt_t2,
            d1: xt.xt_d1,
            d2: xt.xt_d2,
            u1: xt.xt_u1,
            u2: xt.xt_u2,
            ms: xt.xt_ms,
            re: xt.xt_re,
            xf: xt.xt_xf,
        };

        // primary "T" variables at XT, placed into station 2
        self.s2.x = xt.xt;
        self.s2.t = tt;
        self.s2.d = dt;
        self.s2.u = ut;
        self.s2.ampl = self.p.amcrit;
        self.s2.s = 0.0;
        self.blkin();
        self.blvar(Regime::Laminar);
        self.blmid(Regime::Laminar);
        self.bldif(Interval::Laminar);

        // convert "T" sensitivities to "1" and "2"; the amplification row is not needed
        let v = self.sys;
        let mut lam = LocalSystem::default();
        for k in 1..3 {
            let r = &v.vs2[k];
            let comb = |f: fn(&Sens13) -> f64| {
                r[1] * f(&ts) + r[2] * f(&ds) + r[3] * f(&us) + r[4] * f(&xs)
            };
            lam.vsrez[k] = v.vsrez[k];
            lam.vsm[k] = v.vsm[k] + r[1] * ts.ms + r[2] * ds.ms + r[3] * us.ms + r[4] * xs.ms;
            lam.vsr[k] = v.vsr[k] + r[1] * ts.re + r[2] * ds.re + r[3] * us.re + r[4] * xs.re;
            lam.vsx[k] = v.vsx[k] + r[1] * ts.xf + r[2] * ds.xf + r[3] * us.xf + r[4] * xs.xf;
            lam.vs1[k][0] = v.vs1[k][0] + r[1] * ts.a1 + r[2] * ds.a1 + r[3] * us.a1 + r[4] * xs.a1;
            lam.vs1[k][1] = v.vs1[k][1] + r[1] * ts.t1 + r[2] * ds.t1 + r[3] * us.t1 + r[4] * xs.t1;
            lam.vs1[k][2] = v.vs1[k][2] + r[1] * ts.d1 + r[2] * ds.d1 + r[3] * us.d1 + r[4] * xs.d1;
            lam.vs1[k][3] = v.vs1[k][3] + r[1] * ts.u1 + r[2] * ds.u1 + r[3] * us.u1 + r[4] * xs.u1;
            lam.vs1[k][4] = v.vs1[k][4] + r[1] * ts.x1 + r[2] * ds.x1 + r[3] * us.x1 + r[4] * xs.x1;
            lam.vs2[k][0] = 0.0;
            lam.vs2[k][1] = comb(|s| s.t2);
            lam.vs2[k][2] = comb(|s| s.d2);
            lam.vs2[k][3] = comb(|s| s.u2);
            lam.vs2[k][4] = comb(|s| s.x2);
        }

        // ---- turbulent part XT..X2
        // equilibrium shear coefficient CQT at the transition point
        self.blvar(Regime::Turbulent);
        // initial shear coefficient ST at the transition point (CQ2.. are really CQT..)
        let st2 = self.s2;
        let ctr = bl.ctrcon * (-bl.ctrcex / (st2.hk - 1.0)).exp();
        let ctr_hk2 = ctr * bl.ctrcex / powi(st2.hk - 1.0, 2);
        let st = ctr * st2.cq;
        let st_tt = ctr * st2.cq_t + st2.cq * ctr_hk2 * st2.hk_t;
        let st_dt = ctr * st2.cq_d + st2.cq * ctr_hk2 * st2.hk_d;
        let st_ut = ctr * st2.cq_u + st2.cq * ctr_hk2 * st2.hk_u;
        let st_ms0 = ctr * st2.cq_ms + st2.cq * ctr_hk2 * st2.hk_ms;
        let st_re0 = ctr * st2.cq_re;
        let chain = |f: fn(&Sens13) -> f64| st_tt * f(&ts) + st_dt * f(&ds) + st_ut * f(&us);
        let ss = Sens13 {
            a1: chain(|s| s.a1),
            x1: chain(|s| s.x1),
            x2: chain(|s| s.x2),
            t1: chain(|s| s.t1),
            t2: chain(|s| s.t2),
            d1: chain(|s| s.d1),
            d2: chain(|s| s.d2),
            u1: chain(|s| s.u1),
            u2: chain(|s| s.u2),
            ms: chain(|s| s.ms) + st_ms0,
            re: chain(|s| s.re) + st_re0,
            xf: chain(|s| s.xf),
        };
        self.s2.ampl = 0.0;
        self.s2.s = st;
        // turbulent secondary "T" variables with the proper Ctau
        self.blvar(Regime::Turbulent);
        // "1" <- "T", "2" <- saved turbulent station
        self.s1 = self.s2;
        self.s2 = c2sav;
        self.blmid(Regime::Turbulent);
        self.bldif(Interval::Turbulent);

        let v = self.sys;
        let mut tur = LocalSystem::default();
        for k in 0..3 {
            let r = &v.vs1[k];
            let comb = |f: fn(&Sens13) -> f64| {
                r[0] * f(&ss) + r[1] * f(&ts) + r[2] * f(&ds) + r[3] * f(&us) + r[4] * f(&xs)
            };
            tur.vsrez[k] = v.vsrez[k];
            tur.vsm[k] =
                v.vsm[k] + r[0] * ss.ms + r[1] * ts.ms + r[2] * ds.ms + r[3] * us.ms + r[4] * xs.ms;
            tur.vsr[k] =
                v.vsr[k] + r[0] * ss.re + r[1] * ts.re + r[2] * ds.re + r[3] * us.re + r[4] * xs.re;
            tur.vsx[k] =
                v.vsx[k] + r[0] * ss.xf + r[1] * ts.xf + r[2] * ds.xf + r[3] * us.xf + r[4] * xs.xf;
            tur.vs1[k][0] = comb(|s| s.a1);
            tur.vs1[k][1] = comb(|s| s.t1);
            tur.vs1[k][2] = comb(|s| s.d1);
            tur.vs1[k][3] = comb(|s| s.u1);
            tur.vs1[k][4] = comb(|s| s.x1);
            tur.vs2[k][0] = v.vs2[k][0];
            tur.vs2[k][1] = v.vs2[k][1]
                + r[0] * ss.t2
                + r[1] * ts.t2
                + r[2] * ds.t2
                + r[3] * us.t2
                + r[4] * xs.t2;
            tur.vs2[k][2] = v.vs2[k][2]
                + r[0] * ss.d2
                + r[1] * ts.d2
                + r[2] * ds.d2
                + r[3] * us.d2
                + r[4] * xs.d2;
            tur.vs2[k][3] = v.vs2[k][3]
                + r[0] * ss.u2
                + r[1] * ts.u2
                + r[2] * ds.u2
                + r[3] * us.u2
                + r[4] * xs.u2;
            tur.vs2[k][4] = v.vs2[k][4]
                + r[0] * ss.x2
                + r[1] * ts.x2
                + r[2] * ds.x2
                + r[3] * us.x2
                + r[4] * xs.x2;
        }

        // sum laminar and turbulent parts
        let sys = &mut self.sys;
        sys.vsrez[0] = tur.vsrez[0];
        sys.vsrez[1] = lam.vsrez[1] + tur.vsrez[1];
        sys.vsrez[2] = lam.vsrez[2] + tur.vsrez[2];
        sys.vsm[0] = tur.vsm[0];
        sys.vsm[1] = lam.vsm[1] + tur.vsm[1];
        sys.vsm[2] = lam.vsm[2] + tur.vsm[2];
        sys.vsr[0] = tur.vsr[0];
        sys.vsr[1] = lam.vsr[1] + tur.vsr[1];
        sys.vsr[2] = lam.vsr[2] + tur.vsr[2];
        sys.vsx[0] = tur.vsx[0];
        sys.vsx[1] = lam.vsx[1] + tur.vsx[1];
        sys.vsx[2] = lam.vsx[2] + tur.vsx[2];
        for l in 0..5 {
            sys.vs1[0][l] = tur.vs1[0][l];
            sys.vs2[0][l] = tur.vs2[0][l];
            sys.vs1[1][l] = lam.vs1[1][l] + tur.vs1[1][l];
            sys.vs2[1][l] = lam.vs2[1][l] + tur.vs2[1][l];
            sys.vs1[2][l] = lam.vs1[2][l] + tur.vs1[2][l];
            sys.vs2[2][l] = lam.vs2[2][l] + tur.vs2[2][l];
        }
        // restore the clobbered "1" quantities
        self.s1 = c1sav;
    }

    /// Newton system of the current interval, in terms of incompressible Ue and Mach.
    /// Port of XFOIL `BLSYS`.
    pub fn blsys(&mut self) {
        let regime = if self.wake {
            Regime::Wake
        } else if self.turb || self.tran {
            Regime::Turbulent
        } else {
            Regime::Laminar
        };
        self.blvar(regime);
        self.blmid(regime);

        // at the similarity station the "1" and "2" variables are the same
        if self.simi {
            self.s1 = self.s2;
        }

        if self.tran {
            self.trdif();
        } else if self.simi {
            self.bldif(Interval::Similarity);
        } else if !self.turb {
            self.bldif(Interval::Laminar);
        } else if self.wake {
            self.bldif(Interval::Wake);
        } else {
            self.bldif(Interval::Turbulent);
        }

        let sys = &mut self.sys;
        if self.simi {
            // "1" variables are really "2" variables here
            for k in 0..4 {
                for l in 0..5 {
                    sys.vs2[k][l] += sys.vs1[k][l];
                    sys.vs1[k][l] = 0.0;
                }
            }
        }
        // convert residual derivatives wrt compressible Uec to incompressible Uei and Mach
        let (s1, s2) = (&self.s1, &self.s2);
        for k in 0..4 {
            let res_u1 = sys.vs1[k][3];
            let res_u2 = sys.vs2[k][3];
            let res_ms = sys.vsm[k];
            sys.vs1[k][3] = res_u1 * s1.u_uei;
            sys.vs2[k][3] = res_u2 * s2.u_uei;
            sys.vsm[k] = res_u1 * s1.u_ms + res_u2 * s2.u_ms + res_ms;
        }
    }

    /// "Dummy" system between the airfoil TE and the first wake point. Port of XFOIL `TESYS`.
    pub fn tesys(&mut self, cte: f64, tte: f64, dte: f64) {
        self.sys.clear();
        self.blvar(Regime::Wake);
        let (sys, s2) = (&mut self.sys, &self.s2);
        sys.vs1[0][0] = -1.0;
        sys.vs2[0][0] = 1.0;
        sys.vsrez[0] = cte - s2.s;
        sys.vs1[1][1] = -1.0;
        sys.vs2[1][1] = 1.0;
        sys.vsrez[1] = tte - s2.t;
        sys.vs1[2][2] = -1.0;
        sys.vs2[2][2] = 1.0;
        sys.vsrez[2] = dte - s2.d - s2.dw;
    }
}

/// Sensitivities of an interpolated "T" quantity w.r.t. the variables of TRDIF.
#[derive(Debug, Clone, Copy, Default)]
struct Sens13 {
    a1: f64,
    x1: f64,
    x2: f64,
    t1: f64,
    t2: f64,
    d1: f64,
    d2: f64,
    u1: f64,
    u2: f64,
    ms: f64,
    re: f64,
    xf: f64,
}

/// State of the last TRCHEK2 iterate, used for the XT sensitivities.
#[derive(Debug, Clone, Copy, Default)]
struct TrIter {
    wf1: f64,
    wf2: f64,
    wf1_a1: f64,
    wf2_a1: f64,
    wf1_x1: f64,
    wf2_x1: f64,
    wf1_x2: f64,
    wf2_x2: f64,
    wf1_xf: f64,
    wf2_xf: f64,
    xt: f64,
    xt_a2: f64,
    tt_a2: f64,
    dt_a2: f64,
    ut_a2: f64,
    amplt_a2: f64,
    st: super::station::Station,
    a: super::transition::IntervalRate,
}
