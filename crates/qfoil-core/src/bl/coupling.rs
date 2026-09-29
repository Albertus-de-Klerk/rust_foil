//! Viscous–inviscid coupling bookkeeping. Port of XFOIL `STFIND`, `IBLPAN`, `XICALC`,
//! `IBLSYS`, `UICALC`, `QVFUE`, `GAMQV`, `STMOVE`, `UESET`, and QFoil's `VISCAL`
//! initialisation (D9).

use super::{BoundaryLayer, Side};
use crate::fortran::powi;
use crate::inviscid::InviscidSolution;
use crate::paneling::Paneling;
use crate::wake::Wake;

/// Stagnation point between nodes `ist` and `ist+1` (0-based), its arc length and the
/// sensitivities of `sst` to `gam[ist]`, `gam[ist+1]`. Port of XFOIL `STFIND`.
pub fn stfind(gam: &[f64], s: &[f64]) -> (usize, f64, f64, f64) {
    let n = gam.len();
    // Fortran falls back to I = N/2 (1-based) if no sign change is found
    let i = (0..n - 1)
        .find(|&i| gam[i] >= 0.0 && gam[i + 1] < 0.0)
        .unwrap_or(n / 2 - 1);
    let dgam = gam[i + 1] - gam[i];
    let ds = s[i + 1] - s[i];
    // evaluate so as to minimise round-off for very small GAM(I) or GAM(I+1)
    let mut sst = if gam[i] < -gam[i + 1] {
        s[i] - ds * (gam[i] / dgam)
    } else {
        s[i + 1] - ds * (gam[i + 1] / dgam)
    };
    // tweak the stagnation point if it falls right on a node
    if sst <= s[i] {
        sst = s[i] + 1.0e-7;
    }
    if sst >= s[i + 1] {
        sst = s[i + 1] - 1.0e-7;
    }
    let sst_go = (sst - s[i + 1]) / dgam;
    let sst_gp = (s[i] - sst) / dgam;
    (i, sst, sst_go, sst_gp)
}

/// Coordinates of combined node `i` (airfoil `0..n`, then wake).
fn node_xy(pan: &Paneling, wake: &Wake, i: usize) -> (f64, f64) {
    let n = pan.len();
    if i < n {
        (pan.nodes.x[i], pan.nodes.y[i])
    } else {
        (wake.x[i - n], wake.y[i - n])
    }
}

impl BoundaryLayer {
    /// Station capacity per side (Fortran arrays are fixed-size; STMOVE and UPDATE rely on
    /// entries past `NBL`).
    fn capacity(n: usize, nw: usize) -> usize {
        n + nw + 2
    }

    /// Sets up the BL for a fresh operating point: stagnation point, station pointers,
    /// ξ, wake gap, system rows and inviscid edge speed, then QFoil's cold start
    /// (`UEDG = min(UINV, 2)`, `CTAU = 0.01`). The `!LIPAN` block of XFOIL `VISCAL`
    /// followed by `UICALC` and the QFoil initialisation (D9).
    pub fn new(pan: &Paneling, wake: &Wake, sol: &InviscidSolution) -> Self {
        let cap = Self::capacity(pan.len(), wake.len());
        let (ist, sst, sst_go, sst_gp) = stfind(&sol.gam, &pan.nodes.s);
        let mut bl = Self {
            sides: [Side::zeros(cap), Side::zeros(cap)],
            iblte: [0, 0],
            nbl: [0, 0],
            itran: [0, 0],
            xssitr: [0.0, 0.0],
            tforce: [false, false],
            ist,
            sst,
            sst_go,
            sst_gp,
            nsys: 0,
            wgap: vec![0.0; wake.len()],
            xoctr: [1.0, 1.0],
            yoctr: [0.0, 0.0],
        };
        bl.iblpan(pan.len(), wake.len());
        bl.xicalc(pan, wake);
        bl.iblsys();
        bl.uicalc(&sol.qinv, &sol.qinv_a);
        // QFoil D9: always restart from the inviscid Ue, clamping suction peaks at 2.0
        // (includes the dummy station 0, as the Fortran loop runs from IBL = 1)
        for is in 0..2 {
            let sd = &mut bl.sides[is];
            for ibl in 0..bl.nbl[is] {
                sd.uedg[ibl] = if sd.uinv[ibl] > 2.0 {
                    2.0
                } else {
                    sd.uinv[ibl]
                };
                sd.ctau[ibl] = 0.01;
            }
        }
        bl
    }

    /// BL station → panel node pointers. Port of XFOIL `IBLPAN`.
    pub fn iblpan(&mut self, n: usize, nw: usize) {
        let ist = self.ist;
        // upper side: stations 1.. run from the stagnation panel towards node 0
        let mut ibl = 0;
        for i in (0..=ist).rev() {
            ibl += 1;
            self.sides[0].ipan[ibl] = i;
            self.sides[0].vti[ibl] = 1.0;
        }
        self.iblte[0] = ibl;
        self.nbl[0] = ibl + 1;
        // lower side, then the wake
        let mut ibl = 0;
        for i in ist + 1..n {
            ibl += 1;
            self.sides[1].ipan[ibl] = i;
            self.sides[1].vti[ibl] = -1.0;
        }
        self.iblte[1] = ibl;
        for iw in 1..=nw {
            let ibl = self.iblte[1] + iw;
            self.sides[1].ipan[ibl] = n + iw - 1;
            self.sides[1].vti[ibl] = -1.0;
        }
        self.nbl[1] = self.iblte[1] + nw + 1;
        // upper wake pointers (for plotting only)
        for iw in 1..=nw {
            let p = self.sides[1].ipan[self.iblte[1] + iw];
            self.sides[0].ipan[self.iblte[0] + iw] = p;
            self.sides[0].vti[self.iblte[0] + iw] = 1.0;
        }
    }

    /// BL arc length ξ on both sides and the wake, and the TE wake-gap profile.
    /// Port of XFOIL `XICALC`.
    pub fn xicalc(&mut self, pan: &Paneling, wake: &Wake) {
        const XFEPS: f64 = 1.0e-7;
        let s = &pan.nodes.s;
        let n = s.len();
        let xeps = XFEPS * (s[n - 1] - s[0]);
        let sst = self.sst;

        self.sides[0].xssi[0] = 0.0;
        for ibl in 1..=self.iblte[0] {
            let i = self.sides[0].ipan[ibl];
            self.sides[0].xssi[ibl] = (sst - s[i]).max(xeps);
        }
        self.sides[1].xssi[0] = 0.0;
        for ibl in 1..=self.iblte[1] {
            let i = self.sides[1].ipan[ibl];
            self.sides[1].xssi[ibl] = (s[i] - sst).max(xeps);
        }
        let ibl1 = self.iblte[0] + 1;
        self.sides[0].xssi[ibl1] = self.sides[0].xssi[ibl1 - 1];
        let ibl2 = self.iblte[1] + 1;
        self.sides[1].xssi[ibl2] = self.sides[1].xssi[ibl2 - 1];
        for ibl in self.iblte[1] + 2..self.nbl[1] {
            let i = self.sides[1].ipan[ibl];
            let (xa, ya) = node_xy(pan, wake, i);
            let (xb, yb) = node_xy(pan, wake, i - 1);
            let dxssi = ((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb)).sqrt();
            let ibl1 = self.iblte[0] + ibl - self.iblte[1];
            let ibl2 = self.iblte[1] + ibl - self.iblte[1];
            self.sides[0].xssi[ibl1] = self.sides[0].xssi[ibl1 - 1] + dxssi;
            self.sides[1].xssi[ibl2] = self.sides[1].xssi[ibl2 - 1] + dxssi;
        }

        // TE "flap" cubic for the wake gap
        const TELRAT: f64 = 2.50;
        let (xp, yp) = (&pan.nodes.xp, &pan.nodes.yp);
        let crosp = (xp[0] * yp[n - 1] - yp[0] * xp[n - 1])
            / ((powi(xp[0], 2) + powi(yp[0], 2)) * (powi(xp[n - 1], 2) + powi(yp[n - 1], 2)))
                .sqrt();
        let dwdxte = (crosp / (1.0 - powi(crosp, 2)).sqrt())
            .max(-3.0 / TELRAT)
            .min(3.0 / TELRAT);
        let aa = 3.0 + TELRAT * dwdxte;
        let bb = -2.0 - TELRAT * dwdxte;
        let te = &pan.trailing_edge;
        for iw in 1..=self.wgap.len() {
            self.wgap[iw - 1] = if te.sharp {
                0.0
            } else {
                let ibl = self.iblte[1] + iw;
                let sd = &self.sides[1];
                let zn = 1.0 - (sd.xssi[ibl] - sd.xssi[self.iblte[1]]) / (TELRAT * te.ante);
                if zn >= 0.0 {
                    te.ante * (aa + bb * zn) * powi(zn, 2)
                } else {
                    0.0
                }
            };
        }
    }

    /// BL station → Newton-system row pointers. Port of XFOIL `IBLSYS`.
    pub fn iblsys(&mut self) {
        let mut iv = 0;
        for is in 0..2 {
            for ibl in 1..self.nbl[is] {
                self.sides[is].isys[ibl] = iv;
                iv += 1;
            }
        }
        self.nsys = iv;
    }

    /// Inviscid edge speed from the panel speed `qinv` (N+NW). Port of XFOIL `UICALC`.
    pub fn uicalc(&mut self, qinv: &[f64], qinv_a: &[f64]) {
        for is in 0..2 {
            let sd = &mut self.sides[is];
            sd.uinv[0] = 0.0;
            sd.uinv_a[0] = 0.0;
            for ibl in 1..self.nbl[is] {
                let i = sd.ipan[ibl];
                sd.uinv[ibl] = sd.vti[ibl] * qinv[i];
                sd.uinv_a[ibl] = sd.vti[ibl] * qinv_a[i];
            }
        }
    }

    /// Panel viscous speed from BL edge speed. Port of XFOIL `QVFUE`.
    pub fn qvfue(&self, qvis: &mut [f64]) {
        for is in 0..2 {
            let sd = &self.sides[is];
            for ibl in 1..self.nbl[is] {
                qvis[sd.ipan[ibl]] = sd.vti[ibl] * sd.uedg[ibl];
            }
        }
    }

    /// Edge speed from inviscid speed plus all source influence (into `out`, per side).
    /// Port of XFOIL `UESET`.
    pub fn ueset(&self, dij: &crate::linalg::Matrix) -> [Vec<f64>; 2] {
        let mut out = [self.sides[0].uedg.clone(), self.sides[1].uedg.clone()];
        for is in 0..2 {
            let sd = &self.sides[is];
            for ibl in 1..self.nbl[is] {
                let i = sd.ipan[ibl];
                let mut dui = 0.0;
                for js in 0..2 {
                    let sj = &self.sides[js];
                    for jbl in 1..self.nbl[js] {
                        let j = sj.ipan[jbl];
                        let ue_m = -sd.vti[ibl] * sj.vti[jbl] * dij[(i, j)];
                        dui += ue_m * sj.mass[jbl];
                    }
                }
                out[is][ibl] = sd.uinv[ibl] + dui;
            }
        }
        out
    }

    /// Moves the stagnation point to the panel implied by the new vorticity and shifts the
    /// BL arrays accordingly. Port of XFOIL `STMOVE`. Returns `true` if the stagnation
    /// panel changed. `gam`/`qvis` receive XFOIL's tweak of zero edge speeds.
    pub fn stmove(
        &mut self,
        pan: &Paneling,
        wake: &Wake,
        qinv: &[f64],
        qinv_a: &[f64],
        gam: &mut [f64],
        qvis: &mut [f64],
    ) -> bool {
        let istold = self.ist;
        let (ist, sst, sst_go, sst_gp) = stfind(gam, &pan.nodes.s);
        self.ist = ist;
        self.sst = sst;
        self.sst_go = sst_go;
        self.sst_gp = sst_gp;

        let moved = istold != ist;
        if !moved {
            self.xicalc(pan, wake);
        } else {
            self.iblpan(pan.len(), wake.len());
            self.uicalc(qinv, qinv_a);
            self.xicalc(pan, wake);
            self.iblsys();

            let shift = |sd: &mut Side, to: usize, from: usize| {
                sd.ctau[to] = sd.ctau[from];
                sd.thet[to] = sd.thet[from];
                sd.dstr[to] = sd.dstr[from];
                sd.uedg[to] = sd.uedg[from];
            };
            // (grow, shrink): the side gaining stations and the side losing them
            let (grow, shrink, idif) = if ist > istold {
                (0, 1, ist - istold)
            } else {
                (1, 0, istold - ist)
            };
            if ist > istold {
                self.itran[0] += idif;
                self.itran[1] -= idif;
            } else {
                self.itran[0] -= idif;
                self.itran[1] += idif;
            }
            // move the growing side's variables downstream
            {
                let sd = &mut self.sides[grow];
                for ibl in (idif + 1..self.nbl[grow]).rev() {
                    shift(sd, ibl, ibl - idif);
                }
                // fill between the old and new stagnation points
                let dudx = sd.uedg[idif + 1] / sd.xssi[idif + 1];
                for ibl in (1..=idif).rev() {
                    sd.ctau[ibl] = sd.ctau[idif + 1];
                    sd.thet[ibl] = sd.thet[idif + 1];
                    sd.dstr[ibl] = sd.dstr[idif + 1];
                    sd.uedg[ibl] = dudx * sd.xssi[ibl];
                }
            }
            // move the shrinking side's variables upstream
            {
                let sd = &mut self.sides[shrink];
                for ibl in 1..self.nbl[shrink] {
                    shift(sd, ibl, ibl + idif);
                }
            }
            // tweak Ue so it is not zero, in case the stagnation point is right on a node
            const UEPS: f64 = 1.0e-7;
            for is in 0..2 {
                let sd = &mut self.sides[is];
                for ibl in 1..self.nbl[is] {
                    if sd.uedg[ibl] <= UEPS {
                        let i = sd.ipan[ibl];
                        sd.uedg[ibl] = UEPS;
                        qvis[i] = sd.vti[ibl] * UEPS;
                        gam[i] = sd.vti[ibl] * UEPS;
                    }
                }
            }
        }
        // new mass since Ue may have been tweaked
        for is in 0..2 {
            let sd = &mut self.sides[is];
            for ibl in 1..self.nbl[is] {
                sd.mass[ibl] = sd.dstr[ibl] * sd.uedg[ibl];
            }
        }
        moved
    }
}
