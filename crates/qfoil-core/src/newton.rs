//! Global viscous Newton system. Port of XFOIL `SETBL`, `BLSOLV` (xsolve.f) and
//! QFoil `UPDATE` (D6–D8).

use crate::bl::BoundaryLayer;
use crate::bl::march::MarchInputs;
use crate::bl::station::{Kernel, KernelParams, Regime};
use crate::fortran::powi;
use crate::inviscid::forces::{Compressibility, mrcl};
use crate::linalg::Matrix;
use crate::spline::seval;

/// One 3×2 block: `[column][row]`, i.e. Fortran `VA(K, J, IV)` is `va[iv][J-1][K-1]`.
pub type Block = [[f64; 3]; 2];

/// The block-bidiagonal Newton system with a dense mass-defect coupling (COMMON `/VMAT/`).
///
/// Unknowns per row `iv` are (dCτ or dN, dθ, dm). `VA` multiplies row `iv`'s own first
/// two unknowns, `VB` the previous row's, `VM` every row's mass defect, `VZ` the upper-TE
/// row for the first wake point. `VDEL` holds the two right-hand sides (residual, and
/// sensitivity to Re/M or α) and receives the solution.
#[derive(Debug, Clone, PartialEq)]
pub struct NewtonSystem {
    /// Number of rows (`NSYS`).
    pub nsys: usize,
    /// Diagonal blocks.
    pub va: Vec<Block>,
    /// Sub-diagonal blocks.
    pub vb: Vec<Block>,
    /// Right-hand sides / solution.
    pub vdel: Vec<Block>,
    /// Mass coupling: `VM(K, L, IV)` at `vm[(iv*nsys + l)*3 + k]` (Fortran column order).
    pub vm: Vec<f64>,
    /// TE coupling block of the first wake row.
    pub vz: Block,
}

impl NewtonSystem {
    fn zeros(nsys: usize) -> Self {
        Self {
            nsys,
            va: vec![[[0.0; 3]; 2]; nsys],
            vb: vec![[[0.0; 3]; 2]; nsys],
            vdel: vec![[[0.0; 3]; 2]; nsys],
            vm: vec![0.0; 3 * nsys * nsys],
            vz: [[0.0; 3]; 2],
        }
    }

    #[inline(always)]
    fn at(&self, k: usize, l: usize, iv: usize) -> usize {
        (iv * self.nsys + l) * 3 + k
    }

    /// `VM(K, L, IV)` (0-based).
    #[inline(always)]
    pub fn vm(&self, k: usize, l: usize, iv: usize) -> f64 {
        self.vm[self.at(k, l, iv)]
    }

    #[inline(always)]
    fn vm_mut(&mut self, k: usize, l: usize, iv: usize) -> &mut f64 {
        let a = self.at(k, l, iv);
        &mut self.vm[a]
    }

    /// Solves the system in place (solution in `vdel`). Port of XFOIL `BLSOLV`.
    ///
    /// Mass-coupling entries below the drop tolerance `vaccel` (scaled by the airfoil arc
    /// length for rows 2–3) are skipped during elimination, as in XFOIL.
    pub fn blsolv(&mut self, vaccel: f64, arc_length: f64, ivte1: usize, ivz: usize) {
        let nsys = self.nsys;
        let vacc1 = vaccel;
        let vacc2 = vaccel * 2.0 / arc_length;
        let vacc3 = vaccel * 2.0 / arc_length;

        for iv in 0..nsys {
            let ivp = iv + 1;

            // ---- invert the VA(IV) block
            // normalise the first row
            let pivot = 1.0 / self.va[iv][0][0];
            self.va[iv][1][0] *= pivot;
            for l in iv..nsys {
                *self.vm_mut(0, l, iv) *= pivot;
            }
            self.vdel[iv][0][0] *= pivot;
            self.vdel[iv][1][0] *= pivot;

            // eliminate the lower first column
            for k in 1..3 {
                let vtmp = self.va[iv][0][k];
                self.va[iv][1][k] -= vtmp * self.va[iv][1][0];
                for l in iv..nsys {
                    let v = self.vm(0, l, iv);
                    *self.vm_mut(k, l, iv) -= vtmp * v;
                }
                self.vdel[iv][0][k] -= vtmp * self.vdel[iv][0][0];
                self.vdel[iv][1][k] -= vtmp * self.vdel[iv][1][0];
            }

            // normalise the second row
            let pivot = 1.0 / self.va[iv][1][1];
            for l in iv..nsys {
                *self.vm_mut(1, l, iv) *= pivot;
            }
            self.vdel[iv][0][1] *= pivot;
            self.vdel[iv][1][1] *= pivot;

            // eliminate the lower second column
            let k = 2;
            let vtmp = self.va[iv][1][k];
            for l in iv..nsys {
                let v = self.vm(1, l, iv);
                *self.vm_mut(k, l, iv) -= vtmp * v;
            }
            self.vdel[iv][0][k] -= vtmp * self.vdel[iv][0][1];
            self.vdel[iv][1][k] -= vtmp * self.vdel[iv][1][1];

            // normalise the third row
            let pivot = 1.0 / self.vm(2, iv, iv);
            for l in ivp..nsys {
                *self.vm_mut(2, l, iv) *= pivot;
            }
            self.vdel[iv][0][2] *= pivot;
            self.vdel[iv][1][2] *= pivot;

            // eliminate the upper third column
            let vtmp1 = self.vm(0, iv, iv);
            let vtmp2 = self.vm(1, iv, iv);
            for l in ivp..nsys {
                let v3 = self.vm(2, l, iv);
                *self.vm_mut(0, l, iv) -= vtmp1 * v3;
                *self.vm_mut(1, l, iv) -= vtmp2 * v3;
            }
            self.vdel[iv][0][0] -= vtmp1 * self.vdel[iv][0][2];
            self.vdel[iv][0][1] -= vtmp2 * self.vdel[iv][0][2];
            self.vdel[iv][1][0] -= vtmp1 * self.vdel[iv][1][2];
            self.vdel[iv][1][1] -= vtmp2 * self.vdel[iv][1][2];

            // eliminate the upper second column
            let vtmp = self.va[iv][1][0];
            for l in ivp..nsys {
                let v = self.vm(1, l, iv);
                *self.vm_mut(0, l, iv) -= vtmp * v;
            }
            self.vdel[iv][0][0] -= vtmp * self.vdel[iv][0][1];
            self.vdel[iv][1][0] -= vtmp * self.vdel[iv][1][1];

            if iv == nsys - 1 {
                continue;
            }

            // ---- eliminate the VB(IV+1) block, rows 1..3
            for k in 0..3 {
                let vtmp1 = self.vb[ivp][0][k];
                let vtmp2 = self.vb[ivp][1][k];
                let vtmp3 = self.vm(k, iv, ivp);
                for l in ivp..nsys {
                    let d = vtmp1 * self.vm(0, l, iv)
                        + vtmp2 * self.vm(1, l, iv)
                        + vtmp3 * self.vm(2, l, iv);
                    *self.vm_mut(k, l, ivp) -= d;
                }
                for c in 0..2 {
                    let d = vtmp1 * self.vdel[iv][c][0]
                        + vtmp2 * self.vdel[iv][c][1]
                        + vtmp3 * self.vdel[iv][c][2];
                    self.vdel[ivp][c][k] -= d;
                }
            }

            if iv == ivte1 {
                // eliminate the VZ block (upper TE into the first wake row)
                for k in 0..3 {
                    let vtmp1 = self.vz[0][k];
                    let vtmp2 = self.vz[1][k];
                    for l in ivp..nsys {
                        let d = vtmp1 * self.vm(0, l, iv) + vtmp2 * self.vm(1, l, iv);
                        *self.vm_mut(k, l, ivz) -= d;
                    }
                    for c in 0..2 {
                        let d = vtmp1 * self.vdel[iv][c][0] + vtmp2 * self.vdel[iv][c][1];
                        self.vdel[ivz][c][k] -= d;
                    }
                }
            }

            if ivp == nsys - 1 {
                continue;
            }

            // ---- eliminate the lower VM column (skipping entries below the drop tolerance)
            for kv in iv + 2..nsys {
                let vtmp1 = self.vm(0, iv, kv);
                let vtmp2 = self.vm(1, iv, kv);
                let vtmp3 = self.vm(2, iv, kv);
                for (row, vtmp, vacc) in [(0, vtmp1, vacc1), (1, vtmp2, vacc2), (2, vtmp3, vacc3)] {
                    if vtmp.abs() > vacc {
                        for l in ivp..nsys {
                            let v = self.vm(2, l, iv);
                            *self.vm_mut(row, l, kv) -= vtmp * v;
                        }
                        self.vdel[kv][0][row] -= vtmp * self.vdel[iv][0][2];
                        self.vdel[kv][1][row] -= vtmp * self.vdel[iv][1][2];
                    }
                }
            }
        }

        // ---- back substitution: eliminate the upper VM columns
        for iv in (1..nsys).rev() {
            for c in 0..2 {
                let vtmp = self.vdel[iv][c][2];
                for kv in (0..iv).rev() {
                    for row in 0..3 {
                        self.vdel[kv][c][row] -= self.vm(row, iv, kv) * vtmp;
                    }
                }
            }
        }
    }
}

/// Inputs of SETBL that are not BL state.
#[derive(Debug, Clone, Copy)]
pub struct SetblInputs<'a> {
    /// March inputs (panels, Ncrit, trips).
    pub march: MarchInputs<'a>,
    /// Source influence matrix `DIJ`.
    pub dij: &'a Matrix,
    /// Freestream conditions.
    pub flow: &'a crate::settings::FlowConditions,
    /// CL defining M and Re (`CLMR`: current CL for α-prescribed points).
    pub clmr: f64,
    /// α prescribed (`LALFA`); CL-prescribed otherwise.
    pub lalfa: bool,
}

impl BoundaryLayer {
    /// Assembles the global Newton system for the current BL state, marching it first
    /// (MRCHUE on the first call, then MRCHDU). Port of XFOIL `SETBL`.
    ///
    /// `lblini` is XFOIL's "BL initialised" flag; it is set here after MRCHUE.
    pub fn setbl(
        &mut self,
        k: &mut Kernel,
        lblini: &mut bool,
        inp: &SetblInputs<'_>,
    ) -> NewtonSystem {
        let pan = inp.march.pan;
        let dij = inp.dij;
        let qinf = 1.0;

        // M and Re for the current CL, and the freestream constants of the BL equations
        let (minf, ma_clmr, reinf, re_clmr) = mrcl(inp.clmr, inp.flow);
        let msq_clmr = 2.0 * minf * ma_clmr;
        let mut comp = Compressibility {
            minf,
            minf_cl: ma_clmr,
            reinf,
            reinf_cl: re_clmr,
            ..Default::default()
        };
        comp.comset(inp.flow.gamma, qinf);
        k.p = KernelParams::new(
            minf,
            reinf,
            comp.tklam,
            comp.tkl_msq,
            inp.flow.gamma,
            qinf,
            self.wgap[0],
        );

        if !*lblini {
            // initialise the BL by marching with Ue (inverse mode at separation)
            self.mrchue(k, &inp.march);
            *lblini = true;
        }
        // march with the current Ue and Ds to establish transition
        self.mrchdu(k, &inp.march);

        // USAV = UINV + DIJ*MASS; UEDG keeps the marched values
        let usav = self.ueset(dij);

        let nsys = self.nsys;
        let (s0, s1) = (&self.sides[0], &self.sides[1]);
        let (iblte0, iblte1) = (self.iblte[0], self.iblte[1]);
        let ile1 = s0.ipan[1];
        let ile2 = s1.ipan[1];
        let ite1 = s0.ipan[iblte0];
        let ite2 = s1.ipan[iblte1];
        let jvte1 = s0.isys[iblte0];
        let jvte2 = s1.isys[iblte1];
        let dule1 = s0.uedg[1] - usav[0][1];
        let dule2 = s1.uedg[1] - usav[1][1];

        // LE and TE Ue sensitivities w.r.t. all mass values
        let (mut ule1_m, mut ule2_m, mut ute1_m, mut ute2_m) = (
            vec![0.0; nsys],
            vec![0.0; nsys],
            vec![0.0; nsys],
            vec![0.0; nsys],
        );
        for js in 0..2 {
            let sj = &self.sides[js];
            for jbl in 1..self.nbl[js] {
                let j = sj.ipan[jbl];
                let jv = sj.isys[jbl];
                ule1_m[jv] = -s0.vti[1] * sj.vti[jbl] * dij[(ile1, j)];
                ule2_m[jv] = -s1.vti[1] * sj.vti[jbl] * dij[(ile2, j)];
                ute1_m[jv] = -s0.vti[iblte0] * sj.vti[jbl] * dij[(ite1, j)];
                ute2_m[jv] = -s1.vti[iblte1] * sj.vti[jbl] * dij[(ite2, j)];
            }
        }
        let ule1_a = s0.uinv_a[1];
        let ule2_a = s1.uinv_a[1];

        let mut sys = NewtonSystem::zeros(nsys);
        let (mut u1_m, mut d1_m, mut u2_m, mut d2_m) = (
            vec![0.0; nsys],
            vec![0.0; nsys],
            vec![0.0; nsys],
            vec![0.0; nsys],
        );
        // Fortran locals that persist across stations
        let (mut ami, mut cti) = (0.0, 0.0);

        for is in 0..2 {
            // no station "1" at the similarity station
            u1_m.fill(0.0);
            d1_m.fill(0.0);
            let (mut u1_a, mut d1_a) = (0.0, 0.0);
            let (mut due1, mut dds1) = (0.0, 0.0);
            k.p.bule = 1.0;
            k.p.amcrit = inp.march.acrit[is];
            k.p.xiforc = self.xifset(is, &inp.march);
            k.tran = false;
            k.turb = false;

            let (iblte, nbl) = (self.iblte[is], self.nbl[is]);
            let (mut cte_cte1, mut cte_cte2, mut cte_tte1, mut cte_tte2) = (0.0, 0.0, 0.0, 0.0);
            let (tte_tte1, tte_tte2) = (1.0, 1.0);

            for ibl in 1..nbl {
                let iv = self.sides[is].isys[ibl];
                k.simi = ibl == 1;
                k.wake = ibl > iblte;
                k.tran = ibl == self.itran[is];
                k.turb = ibl > self.itran[is];
                let sd = &self.sides[is];
                let i = sd.ipan[ibl];

                // primary variables of the current station
                let xsi = sd.xssi[ibl];
                if ibl < self.itran[is] {
                    ami = sd.ctau[ibl];
                }
                if ibl >= self.itran[is] {
                    cti = sd.ctau[ibl];
                }
                let uei = sd.uedg[ibl];
                let thi = sd.thet[ibl];
                let mdi = sd.mass[ibl];
                let dsi = mdi / uei;
                let dswaki = if k.wake {
                    self.wgap[ibl - iblte - 1]
                } else {
                    0.0
                };

                // derivatives of DSI (= D2)
                let d2_m2 = 1.0 / uei;
                let d2_u2 = -dsi / uei;
                for js in 0..2 {
                    let sj = &self.sides[js];
                    for jbl in 1..self.nbl[js] {
                        let j = sj.ipan[jbl];
                        let jv = sj.isys[jbl];
                        u2_m[jv] = -sd.vti[ibl] * sj.vti[jbl] * dij[(i, j)];
                        d2_m[jv] = d2_u2 * u2_m[jv];
                    }
                }
                d2_m[iv] += d2_m2;
                let u2_a = sd.uinv_a[ibl];
                let d2_a = d2_u2 * u2_a;
                // "forced" changes from the UEDG - USAV mismatch
                let due2 = sd.uedg[ibl] - usav[is][ibl];
                let dds2 = d2_u2 * due2;

                k.blprv(xsi, ami, cti, thi, dsi, dswaki, uei);
                k.blkin();
                // transition interval: set XT etc.
                if k.tran {
                    k.trchek();
                    ami = k.s2.ampl;
                }

                if ibl == iblte + 1 {
                    // start of the wake: TE base thickness added to Dstar
                    let (a, b) = (&self.sides[0], &self.sides[1]);
                    let tte = a.thet[iblte0] + b.thet[iblte1];
                    let dte = a.dstr[iblte0] + b.dstr[iblte1] + pan.trailing_edge.ante;
                    let cte =
                        (a.ctau[iblte0] * a.thet[iblte0] + b.ctau[iblte1] * b.thet[iblte1]) / tte;
                    k.tesys(cte, tte, dte);
                    let dte_mte1 = 1.0 / a.uedg[iblte0];
                    let dte_ute1 = -a.dstr[iblte0] / a.uedg[iblte0];
                    let dte_mte2 = 1.0 / b.uedg[iblte1];
                    let dte_ute2 = -b.dstr[iblte1] / b.uedg[iblte1];
                    cte_cte1 = a.thet[iblte0] / tte;
                    cte_cte2 = b.thet[iblte1] / tte;
                    cte_tte1 = (a.ctau[iblte0] - cte) / tte;
                    cte_tte2 = (b.ctau[iblte1] - cte) / tte;
                    // D1 depends on both TE Ds values
                    for js in 0..2 {
                        let sj = &self.sides[js];
                        for jbl in 1..self.nbl[js] {
                            let jv = sj.isys[jbl];
                            d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
                        }
                    }
                    d1_m[jvte1] += dte_mte1;
                    d1_m[jvte2] += dte_mte2;
                    due1 = 0.0;
                    dds1 = dte_ute1 * (a.uedg[iblte0] - usav[0][iblte0])
                        + dte_ute2 * (b.uedg[iblte1] - usav[1][iblte1]);
                } else {
                    k.blsys();
                }

                // wall shear and equilibrium shear for output
                {
                    let s2 = &k.s2;
                    let sd = &mut self.sides[is];
                    sd.tau[ibl] = 0.5 * s2.r * s2.u * s2.u * s2.cf;
                    sd.dis[ibl] = s2.r * s2.u * s2.u * s2.u * s2.di * s2.hs * 0.5;
                    sd.ctq[ibl] = s2.cq;
                    sd.delt[ibl] = s2.de;
                    sd.uslp[ibl] = 1.60 / (1.0 + s2.us);
                }

                // XI sensitivities w.r.t. LE Ue changes
                let (xi_ule1, xi_ule2) = if is == 0 {
                    (self.sst_go, -self.sst_gp)
                } else {
                    (-self.sst_go, self.sst_gp)
                };

                // BL coefficients into the global Jacobian
                let v = &k.sys;
                for row in 0..3 {
                    let xt = v.vs1[row][4] + v.vs2[row][4] + v.vsx[row];
                    for jv in 0..nsys {
                        *sys.vm_mut(row, jv, iv) = v.vs1[row][2] * d1_m[jv]
                            + v.vs1[row][3] * u1_m[jv]
                            + v.vs2[row][2] * d2_m[jv]
                            + v.vs2[row][3] * u2_m[jv]
                            + (v.vs1[row][4] + v.vs2[row][4] + v.vsx[row])
                                * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
                    }
                    sys.vb[iv][0][row] = v.vs1[row][0];
                    sys.vb[iv][1][row] = v.vs1[row][1];
                    sys.va[iv][0][row] = v.vs2[row][0];
                    sys.va[iv][1][row] = v.vs2[row][1];
                    sys.vdel[iv][1][row] = if inp.lalfa {
                        v.vsr[row] * re_clmr + v.vsm[row] * msq_clmr
                    } else {
                        (v.vs1[row][3] * u1_a + v.vs1[row][2] * d1_a)
                            + (v.vs2[row][3] * u2_a + v.vs2[row][2] * d2_a)
                            + xt * (xi_ule1 * ule1_a + xi_ule2 * ule2_a)
                    };
                    sys.vdel[iv][0][row] = v.vsrez[row]
                        + (v.vs1[row][3] * due1 + v.vs1[row][2] * dds1)
                        + (v.vs2[row][3] * due2 + v.vs2[row][2] * dds2)
                        + xt * (xi_ule1 * dule1 + xi_ule2 * dule2);
                }

                if ibl == iblte + 1 {
                    // TTE, DTE, CTE coefficients
                    for row in 0..3 {
                        sys.vz[0][row] = v.vs1[row][0] * cte_cte1;
                        sys.vz[1][row] = v.vs1[row][0] * cte_tte1 + v.vs1[row][1] * tte_tte1;
                        sys.vb[iv][0][row] = v.vs1[row][0] * cte_cte2;
                        sys.vb[iv][1][row] = v.vs1[row][0] * cte_tte2 + v.vs1[row][1] * tte_tte2;
                    }
                }

                // turbulent intervals follow the transition interval
                if k.tran {
                    k.turb = true;
                    self.itran[is] = ibl;
                    self.tforce[is] = k.trforc;
                    self.xssitr[is] = k.xt.xt;
                    // transition x/c for output
                    let str_ = if is == 0 {
                        self.sst - k.xt.xt
                    } else {
                        self.sst + k.xt.xt
                    };
                    let (xle, yle) = pan.le;
                    let (xte, yte) = pan.te;
                    let chx = xte - xle;
                    let chy = yte - yle;
                    let chsq = powi(chx, 2) + powi(chy, 2);
                    let nd = &pan.nodes;
                    let xtr = seval(str_, &nd.x, &nd.xp, &nd.s);
                    let ytr = seval(str_, &nd.y, &nd.yp, &nd.s);
                    self.xoctr[is] = ((xtr - xle) * chx + (ytr - yle) * chy) / chsq;
                    self.yoctr[is] = ((ytr - yle) * chx - (xtr - xle) * chy) / chsq;
                }
                k.tran = false;

                if ibl == iblte {
                    // "2" variables at the TE become wake correlations for the next station
                    k.turb = true;
                    k.wake = true;
                    k.blvar(Regime::Wake);
                    k.blmid(Regime::Wake);
                }

                u1_m.copy_from_slice(&u2_m);
                d1_m.copy_from_slice(&d2_m);
                u1_a = u2_a;
                d1_a = d2_a;
                due1 = due2;
                dds1 = dds2;

                // next station
                k.s1 = k.s2;
            }
        }
        sys
    }
}

/// Convergence information from [`BoundaryLayer::update`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UpdateResult {
    /// RMS normalised change (`RMSBL`).
    pub rmsbl: f64,
    /// Largest normalised change (`RMXBL`).
    pub rmxbl: f64,
    /// Under-relaxation factor used (`RLX`).
    pub rlx: f64,
    /// Variable with the largest change: `n`, `C`, `T`, `D` or `U` (`VMXBL`).
    pub vmxbl: char,
    /// Station of the largest change (0-based).
    pub imxbl: usize,
    /// Side of the largest change.
    pub ismxbl: usize,
}

/// Inputs of UPDATE.
#[derive(Debug, Clone, Copy)]
pub struct UpdateInputs<'a> {
    /// Airfoil node x.
    pub x: &'a [f64],
    /// Airfoil node y.
    pub y: &'a [f64],
    /// Source influence matrix.
    pub dij: &'a Matrix,
    /// Mach number and dM/dCL.
    pub minf: f64,
    /// dM/dCL.
    pub minf_cl: f64,
    /// `GAMMA - 1`.
    pub gamm1: f64,
    /// α prescribed (`LALFA`).
    pub lalfa: bool,
    /// Prescribed CL for CL-prescribed points (`CLSPEC`).
    pub clspec: f64,
}

impl BoundaryLayer {
    /// Adds the Newton deltas with under-relaxation and returns the RMS/max change.
    /// Updates `cl` (α prescribed) or `alfa` (CL prescribed). Port of QFoil `UPDATE`
    /// (D6: RLX 0.7, D7: CL change ±0.15 without the MATYP guard, D8: DHI 1.0, DLO -0.4).
    pub fn update(
        &mut self,
        sys: &NewtonSystem,
        inp: &UpdateInputs<'_>,
        cl: &mut f64,
        alfa: &mut f64,
    ) -> UpdateResult {
        let qinf = 1.0;
        let dij = inp.dij;
        let minf = inp.minf;
        let gamm1 = inp.gamm1;
        let (x, y) = (inp.x, inp.y);
        let n = x.len();

        // max allowable α and CL changes per iteration (QFoil D7)
        let dalmax = 0.5 * crate::fortran::PI / 180.0;
        let dalmin = -0.5 * crate::fortran::PI / 180.0;
        let dclmax = 0.15;
        let dclmin = -0.15;

        let hstinv = gamm1 * powi(minf / qinf, 2) / (1.0 + 0.5 * gamm1 * powi(minf, 2));

        // new Ue assuming no under-relaxation, and its sensitivity to "AC"
        let cap = self.sides[0].uedg.len();
        let mut unew = [vec![0.0; cap], vec![0.0; cap]];
        let mut u_ac = [vec![0.0; cap], vec![0.0; cap]];
        for is in 0..2 {
            let sd = &self.sides[is];
            for ibl in 1..self.nbl[is] {
                let i = sd.ipan[ibl];
                let mut dui = 0.0;
                let mut dui_ac = 0.0;
                for js in 0..2 {
                    let sj = &self.sides[js];
                    for jbl in 1..self.nbl[js] {
                        let j = sj.ipan[jbl];
                        let jv = sj.isys[jbl];
                        let ue_m = -sd.vti[ibl] * sj.vti[jbl] * dij[(i, j)];
                        dui += ue_m * (sj.mass[jbl] + sys.vdel[jv][0][2]);
                        dui_ac += ue_m * (-sys.vdel[jv][1][2]);
                    }
                }
                let uinv_ac = if inp.lalfa { 0.0 } else { sd.uinv_a[ibl] };
                unew[is][ibl] = sd.uinv[ibl] + dui;
                u_ac[is][ibl] = uinv_ac + dui_ac;
            }
        }

        // new Qtan from the new Ue
        let mut qnew = vec![0.0; n];
        let mut q_ac = vec![0.0; n];
        for is in 0..2 {
            let sd = &self.sides[is];
            for ibl in 1..=self.iblte[is] {
                let i = sd.ipan[ibl];
                qnew[i] = sd.vti[ibl] * unew[is][ibl];
                q_ac[i] = sd.vti[ibl] * u_ac[is][ibl];
            }
        }

        // new CL from the new Qtan
        let sa = alfa.sin();
        let ca = alfa.cos();
        let beta = (1.0 - powi(minf, 2)).sqrt();
        let beta_msq = -0.5 / beta;
        let bfac = 0.5 * powi(minf, 2) / (1.0 + beta);
        let bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;
        let (mut clnew, mut cl_a, mut cl_ms, mut cl_ac) = (0.0, 0.0, 0.0, 0.0);
        let cginc = 1.0 - powi(qnew[0] / qinf, 2);
        let mut cpg1 = cginc / (beta + bfac * cginc);
        let mut cpg1_ms = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
        let cpi_q = -2.0 * qnew[0] / powi(qinf, 2);
        let cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
        let mut cpg1_ac = cpc_cpi * cpi_q * q_ac[0];
        for i in 0..n {
            let ip = if i == n - 1 { 0 } else { i + 1 };
            let cginc = 1.0 - powi(qnew[ip] / qinf, 2);
            let cpg2 = cginc / (beta + bfac * cginc);
            let cpg2_ms = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
            let cpi_q = -2.0 * qnew[ip] / powi(qinf, 2);
            let cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
            let cpg2_ac = cpc_cpi * cpi_q * q_ac[ip];
            let dx = (x[ip] - x[i]) * ca + (y[ip] - y[i]) * sa;
            let dx_a = -(x[ip] - x[i]) * sa + (y[ip] - y[i]) * ca;
            let ag = 0.5 * (cpg2 + cpg1);
            let ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
            let ag_ac = 0.5 * (cpg2_ac + cpg1_ac);
            clnew += dx * ag;
            cl_a += dx_a * ag;
            cl_ms += dx * ag_ms;
            cl_ac += dx * ag_ac;
            cpg1 = cpg2;
            cpg1_ms = cpg2_ms;
            cpg1_ac = cpg2_ac;
        }

        // QFoil D6: initial under-relaxation 0.7 (XFOIL 1.0)
        let mut rlx = 0.7;
        let dac = if inp.lalfa {
            // α prescribed: AC is CL; Re may depend on CL
            let dac = (clnew - *cl) / (1.0 - cl_ac - cl_ms * 2.0 * minf * inp.minf_cl);
            if rlx * dac > dclmax {
                rlx = dclmax / dac;
            }
            if rlx * dac < dclmin {
                rlx = dclmin / dac;
            }
            dac
        } else {
            // CL prescribed: AC is α
            let dac = (clnew - inp.clspec) / (0.0 - cl_ac - cl_a);
            if rlx * dac > dalmax {
                rlx = dalmax / dac;
            }
            if rlx * dac < dalmin {
                rlx = dalmin / dac;
            }
            dac
        };

        let mut rmsbl = 0.0;
        let mut rmxbl: f64 = 0.0;
        let (mut vmxbl, mut imxbl, mut ismxbl) = (' ', 0, 0);
        // QFoil D8: DHI 1.0, DLO -0.4 (XFOIL 1.5, -0.5)
        let dhi = 1.0;
        let dlo = -0.4;

        // changes in BL variables and under-relaxation
        for is in 0..2 {
            let sd = &self.sides[is];
            let itran = self.itran[is];
            for ibl in 1..self.nbl[is] {
                let iv = sd.isys[ibl];
                let d = &sys.vdel[iv];
                let dctau = d[0][0] - dac * d[1][0];
                let dthet = d[0][1] - dac * d[1][1];
                let dmass = d[0][2] - dac * d[1][2];
                let duedg = unew[is][ibl] + dac * u_ac[is][ibl] - sd.uedg[ibl];
                let ddstr = (dmass - sd.dstr[ibl] * duedg) / sd.uedg[ibl];

                let dn1 = if ibl < itran {
                    dctau / 10.0
                } else {
                    dctau / sd.ctau[ibl]
                };
                let dn2 = dthet / sd.thet[ibl];
                let dn3 = ddstr / sd.dstr[ibl];
                let dn4 = duedg.abs() / 0.25;
                rmsbl = rmsbl + powi(dn1, 2) + powi(dn2, 2) + powi(dn3, 2) + powi(dn4, 2);

                let mut track = |dn: f64, val: f64, var: char, rlx: &mut f64| {
                    let rdn = *rlx * dn;
                    if dn.abs() > rmxbl.abs() {
                        rmxbl = val;
                        vmxbl = var;
                        imxbl = ibl;
                        ismxbl = is;
                    }
                    if rdn > dhi {
                        *rlx = dhi / dn;
                    }
                    if rdn < dlo {
                        *rlx = dlo / dn;
                    }
                };
                track(dn1, dn1, if ibl < itran { 'n' } else { 'C' }, &mut rlx);
                track(dn2, dn2, 'T', &mut rlx);
                track(dn3, dn3, 'D', &mut rlx);
                track(dn4, duedg, 'U', &mut rlx);
            }
        }
        rmsbl = (rmsbl / (4.0 * (self.nbl[0] + self.nbl[1]) as f64)).sqrt();

        if inp.lalfa {
            *cl += rlx * dac;
        } else {
            *alfa += rlx * dac;
        }

        // update BL variables with the under-relaxed changes
        for is in 0..2 {
            let (iblte, itran, nbl) = (self.iblte[is], self.itran[is], self.nbl[is]);
            for ibl in 1..nbl {
                let sd = &mut self.sides[is];
                let iv = sd.isys[ibl];
                let d = &sys.vdel[iv];
                let dctau = d[0][0] - dac * d[1][0];
                let dthet = d[0][1] - dac * d[1][1];
                let dmass = d[0][2] - dac * d[1][2];
                let duedg = unew[is][ibl] + dac * u_ac[is][ibl] - sd.uedg[ibl];
                let ddstr = (dmass - sd.dstr[ibl] * duedg) / sd.uedg[ibl];
                sd.ctau[ibl] += rlx * dctau;
                sd.thet[ibl] += rlx * dthet;
                sd.dstr[ibl] += rlx * ddstr;
                sd.uedg[ibl] += rlx * duedg;
                let dswaki = if ibl > iblte {
                    self.wgap[ibl - iblte - 1]
                } else {
                    0.0
                };
                let sd = &mut self.sides[is];
                // eliminate absurd transients
                if ibl >= itran {
                    sd.ctau[ibl] = sd.ctau[ibl].min(0.25);
                }
                let hklim = if ibl <= iblte { 1.02 } else { 1.00005 };
                let msq = powi(sd.uedg[ibl], 2) * hstinv
                    / (gamm1 * (1.0 - 0.5 * powi(sd.uedg[ibl], 2) * hstinv));
                let dsw = sd.dstr[ibl] - dswaki;
                sd.dstr[ibl] = crate::bl::march::dslim(dsw, sd.thet[ibl], msq, hklim) + dswaki;
                // new mass defect (nonlinear update)
                sd.mass[ibl] = sd.dstr[ibl] * sd.uedg[ibl];
            }
            // no "islands" of negative Ue (Fortran IBL = 3..IBLTE)
            let sd = &mut self.sides[is];
            for ibl in 2..=iblte {
                if sd.uedg[ibl - 1] > 0.0 && sd.uedg[ibl] <= 0.0 {
                    sd.uedg[ibl] = sd.uedg[ibl - 1];
                    sd.mass[ibl] = sd.dstr[ibl] * sd.uedg[ibl];
                }
            }
        }

        // upper wake arrays mirror the lower wake arrays
        let (i0, i1) = (self.iblte[0], self.iblte[1]);
        for kbl in 1..self.nbl[1] - i1 {
            let src = self.sides[1].clone_station(i1 + kbl);
            self.sides[0].set_station(i0 + kbl, &src);
        }

        UpdateResult {
            rmsbl,
            rmxbl,
            rlx,
            vmxbl,
            imxbl,
            ismxbl,
        }
    }
}
