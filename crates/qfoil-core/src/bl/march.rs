//! Station-by-station BL marching. Port of XFOIL `MRCHUE`, `MRCHDU`, `XIFSET`, `DSLIM`.

use super::closure::hkin;
use super::station::{Kernel, Regime, Rows};
use super::{BoundaryLayer, Side};
use crate::fortran::{pow, powi};
use crate::linalg::gauss;
use crate::paneling::Paneling;
use crate::spline::{EndCondition, sinvrt, splind};

/// Inputs of a march that are not BL state.
#[derive(Debug, Clone, Copy)]
pub struct MarchInputs<'a> {
    /// Airfoil panels (TE gap `ANTE`, geometry for forced transition).
    pub pan: &'a Paneling,
    /// Critical amplification exponent per side (`ACRIT`).
    pub acrit: [f64; 2],
    /// Forced-transition x/c per side (`XSTRIP`; ≥ 1 means free transition).
    pub xstrip: [f64; 2],
}

/// Raises δ* so that Hk ≥ `hklim`. Port of XFOIL `DSLIM`.
pub fn dslim(dstr: f64, thet: f64, msq: f64, hklim: f64) -> f64 {
    let h = dstr / thet;
    let (hk, hk_h, _) = hkin(h, msq);
    let dh = (hklim - hk).max(0.0) / hk_h;
    dstr + dh * thet
}

/// Edge Mach number squared for incompressible edge speed `uei`.
fn edge_msq(k: &Kernel, uei: f64) -> f64 {
    uei * uei * k.p.hstinv / (k.p.gm1bl * (1.0 - 0.5 * uei * uei * k.p.hstinv))
}

impl BoundaryLayer {
    /// BL ξ of the forced-transition point on `side`. Port of XFOIL `XIFSET`.
    pub fn xifset(&self, side: usize, inp: &MarchInputs<'_>) -> f64 {
        let sd = &self.sides[side];
        let xi_te = sd.xssi[self.iblte[side]];
        let xstrip = inp.xstrip[side];
        if xstrip >= 1.0 {
            return xi_te;
        }
        let pan = inp.pan;
        let (x, y, s) = (&pan.nodes.x, &pan.nodes.y, &pan.nodes.s);
        let n = x.len();
        let (xle, yle) = pan.le;
        let (xte, yte) = pan.te;
        let chx = xte - xle;
        let chy = yte - yle;
        let chsq = chx * chx + chy * chy;
        // chord-based x/c (W1) and y/c (W2), splined over the whole airfoil
        let w1: Vec<f64> = (0..n)
            .map(|i| ((x[i] - xle) * chx + (y[i] - yle) * chy) / chsq)
            .collect();
        let w2: Vec<f64> = (0..n)
            .map(|i| ((y[i] - yle) * chx - (x[i] - xle) * chy) / chsq)
            .collect();
        let (mut w3, mut w4) = (vec![0.0; n], vec![0.0; n]);
        splind(
            &w1,
            &mut w3,
            s,
            EndCondition::ZeroThirdDerivative,
            EndCondition::ZeroThirdDerivative,
        );
        splind(
            &w2,
            &mut w4,
            s,
            EndCondition::ZeroThirdDerivative,
            EndCondition::ZeroThirdDerivative,
        );
        let sle = pan.sle;
        let xiforc = if side == 0 {
            let str_ = sle + (s[0] - sle) * xstrip;
            let str_ = sinvrt(str_, xstrip, &w1, &w3, s).s;
            (self.sst - str_).min(xi_te)
        } else {
            let str_ = sle + (s[n - 1] - sle) * xstrip;
            let str_ = sinvrt(str_, xstrip, &w1, &w3, s).s;
            (str_ - self.sst).min(xi_te)
        };
        // stagnation point past the trip: no forcing
        if xiforc < 0.0 { xi_te } else { xiforc }
    }

    /// TE quantities feeding the first wake station: (CTE, TTE, DTE).
    fn te_values(&self, ante: f64) -> (f64, f64, f64) {
        let (a, b) = (&self.sides[0], &self.sides[1]);
        let (i1, i2) = (self.iblte[0], self.iblte[1]);
        let tte = a.thet[i1] + b.thet[i2];
        let dte = a.dstr[i1] + b.dstr[i2] + ante;
        let cte = (a.ctau[i1] * a.thet[i1] + b.ctau[i2] * b.thet[i2]) / tte;
        (cte, tte, dte)
    }

    /// Stores the converged station-2 variables of `ibl` (label 110 in MRCHUE/MRCHDU).
    fn store(
        side: &mut Side,
        ibl: usize,
        k: &Kernel,
        laminar: bool,
        ami: f64,
        cti: f64,
        thi: f64,
        dsi: f64,
        uei: f64,
    ) {
        let s2 = &k.s2;
        side.ctau[ibl] = if laminar { ami } else { cti };
        side.thet[ibl] = thi;
        side.dstr[ibl] = dsi;
        side.uedg[ibl] = uei;
        side.mass[ibl] = dsi * uei;
        side.tau[ibl] = 0.5 * s2.r * s2.u * s2.u * s2.cf;
        side.dis[ibl] = s2.r * s2.u * s2.u * s2.u * s2.di * s2.hs * 0.5;
        side.ctq[ibl] = s2.cq;
        side.delt[ibl] = s2.de;
        side.tstr[ibl] = s2.hs * s2.t;
    }

    /// Recomputes secondary variables for extrapolated station values (label 109).
    fn reset_station(&mut self, k: &mut Kernel, is: usize, ibl: usize, ami: &mut f64, v: [f64; 5]) {
        let [xsi, cti, thi, dsi, uei] = v;
        let dswaki = self.dswak(is, ibl);
        k.blprv(xsi, *ami, cti, thi, dsi, dswaki, uei);
        k.blkin();
        if !k.simi && !k.turb {
            k.trchek();
            *ami = k.s2.ampl;
            self.itran[is] = if k.tran { ibl } else { ibl + 2 };
        }
        let itran = self.itran[is];
        if ibl < itran {
            k.blvar(Regime::Laminar);
        }
        if ibl >= itran {
            k.blvar(Regime::Turbulent);
        }
        if k.wake {
            k.blvar(Regime::Wake);
        }
        if ibl < itran {
            k.blmid(Regime::Laminar);
        }
        if ibl >= itran {
            k.blmid(Regime::Turbulent);
        }
        if k.wake {
            k.blmid(Regime::Wake);
        }
    }

    /// Wake gap δ*_w at station `ibl` (0 on the airfoil).
    fn dswak(&self, is: usize, ibl: usize) -> f64 {
        if ibl > self.iblte[is] {
            self.wgap[ibl - self.iblte[is] - 1]
        } else {
            0.0
        }
    }

    /// Marches both BLs and the wake in direct mode with the current Ue, switching to
    /// inverse (prescribed Hk) mode at separation. Port of XFOIL `MRCHUE`.
    pub fn mrchue(&mut self, k: &mut Kernel, inp: &MarchInputs<'_>) {
        const HLMAX: f64 = 3.8;
        const HTMAX: f64 = 2.5;
        let ante = inp.pan.trailing_edge.ante;
        let hstinv = k.p.hstinv;
        let gm1bl = k.p.gm1bl;
        // Locals that Fortran keeps across stations (declared once per call)
        let (mut cte, mut tte, mut dte) = (0.0, 0.0, 0.0);
        let mut hmax = 0.0;

        for is in 0..2 {
            k.p.amcrit = inp.acrit[is];
            k.p.xiforc = self.xifset(is, inp);

            // similarity station initialised with Thwaites' formula
            let ibl = 1;
            let sd = &self.sides[is];
            let xsi = sd.xssi[ibl];
            let uei = sd.uedg[ibl];
            k.p.bule = 1.0;
            let bule = k.p.bule;
            let ucon = uei / pow(xsi, bule);
            let tsq = 0.45 / (ucon * (5.0 * bule + 1.0) * k.p.reybl) * pow(xsi, 1.0 - bule);
            let mut thi = tsq.sqrt();
            let mut dsi = 2.2 * thi;
            let mut ami = 0.0;
            let mut cti = 0.03;

            k.tran = false;
            k.turb = false;
            self.itran[is] = self.iblte[is];

            let (iblte, nbl) = (self.iblte[is], self.nbl[is]);
            for ibl in 1..nbl {
                let ibm = ibl - 1;
                k.simi = ibl == 1;
                k.wake = ibl > iblte;
                let xsi = self.sides[is].xssi[ibl];
                let mut uei = self.sides[is].uedg[ibl];
                let dswaki = self.dswak(is, ibl);
                let mut direct = true;
                let mut htarg = 0.0;
                let mut dmax = 0.0;
                let mut converged = false;

                for _itbl in 0..25 {
                    // 4x5 linearised system at the current station "2"
                    k.blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
                    k.blkin();

                    if !k.simi && !k.turb {
                        k.trchek();
                        ami = k.s2.ampl;
                        if k.tran {
                            self.itran[is] = ibl;
                            if cti <= 0.0 {
                                cti = 0.03;
                                k.s2.s = cti;
                            }
                        } else {
                            self.itran[is] = ibl + 2;
                        }
                    }
                    let itran = self.itran[is];

                    if ibl == iblte + 1 {
                        (cte, tte, dte) = self.te_values(ante);
                        k.tesys(cte, tte, dte);
                    } else {
                        k.blsys();
                    }

                    if direct {
                        // try direct mode: dUe = 0 in the empty 4th row
                        k.sys.vs2[3][0] = 0.0;
                        k.sys.vs2[3][1] = 0.0;
                        k.sys.vs2[3][2] = 0.0;
                        k.sys.vs2[3][3] = 1.0;
                        k.sys.vsrez[3] = 0.0;
                        gauss(4, &mut Rows(&mut k.sys.vs2), &mut k.sys.vsrez);
                        let r = k.sys.vsrez;

                        dmax = (r[1] / thi).abs().max((r[2] / dsi).abs());
                        if ibl < itran {
                            dmax = dmax.max((r[0] / 10.0).abs());
                        }
                        if ibl >= itran {
                            dmax = dmax.max((r[0] / cti).abs());
                        }
                        let mut rlx = 1.0;
                        if dmax > 0.3 {
                            rlx = 0.3 / dmax;
                        }

                        // is direct mode still applicable?
                        if ibl != iblte + 1 {
                            let msq =
                                uei * uei * hstinv / (gm1bl * (1.0 - 0.5 * uei * uei * hstinv));
                            let htest = (dsi + rlx * r[2]) / (thi + rlx * r[1]);
                            let (hktest, _, _) = hkin(htest, msq);
                            if ibl < itran {
                                hmax = HLMAX;
                            }
                            if ibl >= itran {
                                hmax = HTMAX;
                            }
                            direct = hktest < hmax;
                        }

                        if direct {
                            if ibl >= itran {
                                cti += rlx * r[0];
                            }
                            thi += rlx * r[1];
                            dsi += rlx * r[2];
                        } else {
                            // prescribed Hk for inverse mode at this station
                            let (s1, s2) = (&k.s1, &k.s2);
                            htarg = if ibl < itran {
                                // laminar: slow increase of Hk downstream
                                s1.hk + 0.03 * (s2.x - s1.x) / s1.t
                            } else if ibl == itran {
                                // transition interval: weighted laminar and turbulent
                                s1.hk + (0.03 * (k.xt.xt - s1.x) - 0.15 * (s2.x - k.xt.xt)) / s1.t
                            } else if k.wake {
                                // wake: asymptotic behaviour, three Backward-Euler Newton steps
                                let cnst = 0.03 * (s2.x - s1.x) / s1.t;
                                let mut hk2 = s1.hk;
                                for _ in 0..3 {
                                    // CONST*(HK2-1)**3, 3*CONST*(HK2-1)**2: power first (gfortran)
                                    hk2 = hk2
                                        - (hk2 + cnst * powi(hk2 - 1.0, 3) - s1.hk)
                                            / (1.0 + 3.0 * cnst * powi(hk2 - 1.0, 2));
                                }
                                k.s2.hk = hk2; // Fortran writes HK2 here
                                hk2
                            } else {
                                // turbulent: fast decrease of Hk downstream
                                s1.hk - 0.15 * (s2.x - s1.x) / s1.t
                            };
                            htarg = if k.wake {
                                htarg.max(1.01)
                            } else {
                                htarg.max(hmax)
                            };
                            continue; // GO TO 100: retry this station in inverse mode
                        }
                    } else {
                        // inverse mode: force Hk to HTARG
                        k.sys.vs2[3][0] = 0.0;
                        k.sys.vs2[3][1] = k.s2.hk_t;
                        k.sys.vs2[3][2] = k.s2.hk_d;
                        k.sys.vs2[3][3] = k.s2.hk_u;
                        k.sys.vsrez[3] = htarg - k.s2.hk;
                        gauss(4, &mut Rows(&mut k.sys.vs2), &mut k.sys.vsrez);
                        let r = k.sys.vsrez;
                        // (Ue clamp added MD 3 Apr 03)
                        dmax = (r[1] / thi)
                            .abs()
                            .max((r[2] / dsi).abs())
                            .max((r[3] / uei).abs());
                        if ibl >= itran {
                            dmax = dmax.max((r[0] / cti).abs());
                        }
                        let mut rlx = 1.0;
                        if dmax > 0.3 {
                            rlx = 0.3 / dmax;
                        }
                        if ibl >= itran {
                            cti += rlx * r[0];
                        }
                        thi += rlx * r[1];
                        dsi += rlx * r[2];
                        uei += rlx * r[3];
                    }

                    // eliminate absurd transients
                    if ibl >= itran {
                        cti = cti.min(0.30).max(0.0000001);
                    }
                    let hklim = if ibl <= iblte { 1.02 } else { 1.00005 };
                    let msq = edge_msq(k, uei);
                    let dsw = dsi - dswaki;
                    dsi = dslim(dsw, thi, msq, hklim) + dswaki;

                    if dmax <= 1.0e-5 {
                        converged = true;
                        break;
                    }
                }

                if !converged {
                    // MRCHUE: convergence failed (the Fortran prints a message)
                    if dmax > 0.1 {
                        // the current solution is garbage: extrapolate instead
                        let sd = &self.sides[is];
                        let itran = self.itran[is];
                        // Fortran IBL.GT.3  ->  ibl > 2
                        if ibl > 2 {
                            if ibl <= iblte {
                                thi = sd.thet[ibm] * pow(sd.xssi[ibl] / sd.xssi[ibm], 0.5);
                                dsi = sd.dstr[ibm] * pow(sd.xssi[ibl] / sd.xssi[ibm], 0.5);
                            } else if ibl == iblte + 1 {
                                cti = cte;
                                thi = tte;
                                dsi = dte;
                            } else {
                                thi = sd.thet[ibm];
                                let ratlen = (sd.xssi[ibl] - sd.xssi[ibm]) / (10.0 * sd.dstr[ibm]);
                                dsi = (sd.dstr[ibm] + thi * ratlen) / (1.0 + ratlen);
                            }
                            if ibl == itran {
                                cti = 0.05;
                            }
                            if ibl > itran {
                                cti = sd.ctau[ibm];
                            }
                            uei = sd.uedg[ibl];
                            // Fortran IBL.GT.2 .AND. IBL.LT.NBL
                            if ibl > 1 && ibl < nbl - 1 {
                                uei = 0.5 * (sd.uedg[ibl - 1] + sd.uedg[ibl + 1]);
                            }
                        }
                    }
                    self.reset_station(k, is, ibl, &mut ami, [xsi, cti, thi, dsi, uei]);
                }

                // store primary variables (label 110)
                let laminar = ibl < self.itran[is];
                Self::store(
                    &mut self.sides[is],
                    ibl,
                    k,
                    laminar,
                    ami,
                    cti,
                    thi,
                    dsi,
                    uei,
                );

                // "1" <- "2" for the next station
                k.blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
                k.blkin();
                k.s1 = k.s2;

                // turbulent intervals follow the transition interval or the TE
                if k.tran || ibl == iblte {
                    k.turb = true;
                    self.tforce[is] = k.trforc;
                    self.xssitr[is] = k.xt.xt;
                }
                k.tran = false;

                if ibl == iblte {
                    let (a, b) = (&self.sides[0], &self.sides[1]);
                    thi = a.thet[self.iblte[0]] + b.thet[self.iblte[1]];
                    dsi = a.dstr[self.iblte[0]] + b.dstr[self.iblte[1]] + ante;
                }
            }
        }
    }

    /// Marches both BLs and the wake in mixed mode along the Ue–Hk line quasi-normal to
    /// the local characteristic. Port of XFOIL/QFoil `MRCHDU` (QFoil D5: RLX 0.7).
    pub fn mrchdu(&mut self, k: &mut Kernel, inp: &MarchInputs<'_>) {
        const DEPS: f64 = 5.0e-6;
        // how far Hk may deviate from the specified value
        const SENSWT: f64 = 1000.0;
        let ante = inp.pan.trailing_edge.ante;
        // Locals that Fortran keeps across stations (declared once per call)
        let (mut cte, mut tte, mut dte) = (0.0, 0.0, 0.0);
        let (mut sens, mut sennew) = (0.0, 0.0);
        let mut ami = 0.0;
        let (mut ueref, mut hkref) = (0.0, 0.0);

        for is in 0..2 {
            k.p.amcrit = inp.acrit[is];
            k.p.xiforc = self.xifset(is, inp);
            k.p.bule = 1.0;

            let itrold = self.itran[is];
            k.tran = false;
            k.turb = false;
            self.itran[is] = self.iblte[is];

            let (iblte, nbl) = (self.iblte[is], self.nbl[is]);
            for ibl in 1..nbl {
                let ibm = ibl - 1;
                k.simi = ibl == 1;
                k.wake = ibl > iblte;

                // initialise the station to the existing variables
                let sd = &self.sides[is];
                let xsi = sd.xssi[ibl];
                let mut uei = sd.uedg[ibl];
                let mut thi = sd.thet[ibl];
                let mut dsi = sd.dstr[ibl];
                let mut cti;
                // (fixed bug, MD 7 June 99)
                if ibl < itrold {
                    ami = sd.ctau[ibl];
                    cti = 0.03;
                } else {
                    cti = sd.ctau[ibl];
                    if cti <= 0.0 {
                        cti = 0.03;
                    }
                }
                let dswaki = self.dswak(is, ibl);
                if ibl <= iblte {
                    dsi = (dsi - dswaki).max(1.02000 * thi) + dswaki;
                }
                if ibl > iblte {
                    dsi = (dsi - dswaki).max(1.00005 * thi) + dswaki;
                }

                let mut dmax = 0.0;
                let mut converged = false;
                for itbl in 1..=25 {
                    k.blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
                    k.blkin();

                    if !k.simi && !k.turb {
                        k.trchek();
                        ami = k.s2.ampl;
                        self.itran[is] = if k.tran { ibl } else { ibl + 2 };
                    }

                    if ibl == iblte + 1 {
                        (cte, tte, dte) = self.te_values(ante);
                        k.tesys(cte, tte, dte);
                    } else {
                        k.blsys();
                    }
                    let itran = self.itran[is];

                    if itbl == 1 {
                        // baseline Ue and Hk for the Ue(Hk) relation
                        ueref = k.s2.u;
                        hkref = k.s2.hk;
                        // station was turbulent and is now laminar: extrapolate baseline Hk
                        if ibl < itran && ibl >= itrold {
                            let sd = &self.sides[is];
                            let uem = sd.uedg[ibl - 1];
                            let dsm = sd.dstr[ibl - 1];
                            let thm = sd.thet[ibl - 1];
                            let msq = edge_msq(k, uem);
                            hkref = hkin(dsm / thm, msq).0;
                        }
                        // station was laminar: reinitialise or extrapolate Ctau if now turbulent
                        if ibl < itrold {
                            if k.tran {
                                self.sides[is].ctau[ibl] = 0.03;
                            }
                            if k.turb {
                                self.sides[is].ctau[ibl] = self.sides[is].ctau[ibl - 1];
                            }
                            if k.tran || k.turb {
                                cti = self.sides[is].ctau[ibl];
                                k.s2.s = cti;
                            }
                        }
                    }

                    if k.simi || ibl == iblte + 1 {
                        // similarity station or first wake point: prescribe Ue
                        k.sys.vs2[3][0] = 0.0;
                        k.sys.vs2[3][1] = 0.0;
                        k.sys.vs2[3][2] = 0.0;
                        k.sys.vs2[3][3] = k.s2.u_uei;
                        k.sys.vsrez[3] = ueref - k.s2.u;
                    } else {
                        // Ue-Hk characteristic slope: dUe response to a unit dHk
                        let mut vtmp = k.sys.vs2;
                        let mut vztmp = k.sys.vsrez;
                        vtmp[3][0] = 0.0;
                        vtmp[3][1] = k.s2.hk_t;
                        vtmp[3][2] = k.s2.hk_d;
                        vtmp[3][3] = k.s2.hk_u * k.s2.u_uei;
                        vztmp[3] = 1.0;
                        gauss(4, &mut Rows(&mut vtmp), &mut vztmp);

                        // SENSWT * (normalised dUe/dHk)
                        sennew = SENSWT * vztmp[3] * hkref / ueref;
                        if itbl <= 5 {
                            sens = sennew;
                        } else if itbl <= 15 {
                            sens = 0.5 * (sens + sennew);
                        }

                        // prescribed Ue-Hk combination
                        let s2 = &k.s2;
                        k.sys.vs2[3][0] = 0.0;
                        k.sys.vs2[3][1] = s2.hk_t * hkref;
                        k.sys.vs2[3][2] = s2.hk_d * hkref;
                        k.sys.vs2[3][3] = (s2.hk_u * hkref + sens / ueref) * s2.u_uei;
                        k.sys.vsrez[3] =
                            -(hkref * hkref) * (s2.hk / hkref - 1.0) - sens * (s2.u / ueref - 1.0);
                    }

                    gauss(4, &mut Rows(&mut k.sys.vs2), &mut k.sys.vsrez);
                    let r = k.sys.vsrez;

                    // max changes and under-relaxation (Ue clamp added MD 3 Apr 03)
                    dmax = (r[1] / thi)
                        .abs()
                        .max((r[2] / dsi).abs())
                        .max((r[3] / uei).abs());
                    if ibl >= itran {
                        dmax = dmax.max((r[0] / (10.0 * cti)).abs());
                    }
                    // QFoil D5: RLX = 0.7 (XFOIL 1.0); then 0.3/DMAX may exceed it (S1)
                    let mut rlx = 0.7;
                    if dmax > 0.3 {
                        rlx = 0.3 / dmax;
                    }

                    if ibl < itran {
                        ami += rlx * r[0];
                    }
                    if ibl >= itran {
                        cti += rlx * r[0];
                    }
                    thi += rlx * r[1];
                    dsi += rlx * r[2];
                    uei += rlx * r[3];

                    // eliminate absurd transients
                    if ibl >= itran {
                        cti = cti.min(0.30).max(0.0000001);
                    }
                    let hklim = if ibl <= iblte { 1.02 } else { 1.00005 };
                    let msq = edge_msq(k, uei);
                    let dsw = dsi - dswaki;
                    dsi = dslim(dsw, thi, msq, hklim) + dswaki;

                    if dmax <= DEPS {
                        converged = true;
                        break;
                    }
                }

                if !converged {
                    // MRCHDU: convergence failed (the Fortran prints a message)
                    if dmax > 0.1 {
                        // the current solution is garbage: extrapolate instead
                        let sd = &self.sides[is];
                        let itran = self.itran[is];
                        // Fortran IBL.GT.3  ->  ibl > 2
                        if ibl > 2 {
                            if ibl <= iblte {
                                thi = sd.thet[ibm] * pow(sd.xssi[ibl] / sd.xssi[ibm], 0.5);
                                dsi = sd.dstr[ibm] * pow(sd.xssi[ibl] / sd.xssi[ibm], 0.5);
                                uei = sd.uedg[ibm];
                            } else if ibl == iblte + 1 {
                                cti = cte;
                                thi = tte;
                                dsi = dte;
                                uei = sd.uedg[ibm];
                            } else {
                                thi = sd.thet[ibm];
                                let ratlen = (sd.xssi[ibl] - sd.xssi[ibm]) / (10.0 * sd.dstr[ibm]);
                                dsi = (sd.dstr[ibm] + thi * ratlen) / (1.0 + ratlen);
                                uei = sd.uedg[ibm];
                            }
                            if ibl == itran {
                                cti = 0.05;
                            }
                            if ibl > itran {
                                cti = sd.ctau[ibm];
                            }
                        }
                    }
                    self.reset_station(k, is, ibl, &mut ami, [xsi, cti, thi, dsi, uei]);
                }

                // label 110
                sens = sennew;
                let laminar = ibl < self.itran[is];
                Self::store(
                    &mut self.sides[is],
                    ibl,
                    k,
                    laminar,
                    ami,
                    cti,
                    thi,
                    dsi,
                    uei,
                );

                k.blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
                k.blkin();
                k.s1 = k.s2;

                if k.tran || ibl == iblte {
                    k.turb = true;
                    self.tforce[is] = k.trforc;
                    self.xssitr[is] = k.xt.xt;
                }
                k.tran = false;
            }
        }
        let _ = sens;
    }
}
