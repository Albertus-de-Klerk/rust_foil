//! Inviscid panel solution. Port of XFOIL `GGCALC`, `SPECAL`, `QISET`, `TECALC`
//! (strengths), `QDCALC`, with the kernels in [`influence`] and forces in [`forces`].

pub mod forces;
pub mod influence;

use forces::{Compressibility, clcalc, cpcalc, mrcl};
use influence::{FieldPoint, Influence, NodeRef, Strengths, psilin, pswlin};

use crate::linalg::{LinalgError, LuFactors, Matrix, ludcmp};
use crate::operating::OperatingPoint;
use crate::paneling::Paneling;
use crate::settings::{MachType, Settings};
use crate::wake::Wake;

/// Distance of the sharp-TE internal control point, as a fraction of the smaller TE panel
/// (`BWT`).
const BWT: f64 = 0.1;

/// The inviscid solution (COMMON `/CR03/`, `/CR04/`, `/CR06/` strengths).
#[derive(Debug, Clone, PartialEq)]
pub struct InviscidSolution {
    /// LU-factored influence matrix `AIJ` ((N+1)², incl. the Kutta row and ψ₀ column).
    pub aij: LuFactors,
    /// Source right-hand sides `BIJ` for the airfoil sources ((N+1) × N), not back-solved.
    pub bij: Matrix,
    /// Unit vorticity for α = 0°, 90° (N+1 each; row N is ψ₀).
    pub gamu: [Vec<f64>; 2],
    /// Angle of attack of `gam` (rad).
    pub alfa: f64,
    /// Surface vorticity `GAM` (N).
    pub gam: Vec<f64>,
    /// dGAM/dα (N).
    pub gam_a: Vec<f64>,
    /// Stream function inside the airfoil (`PSIO`).
    pub psio: f64,
    /// TE panel source strength (`SIGTE`).
    pub sigte: f64,
    /// TE panel vortex strength (`GAMTE`).
    pub gamte: f64,
    /// Inviscid tangential speed for α = 0°, 90° on airfoil and wake (`QINVU`, N+NW).
    pub qinvu: [Vec<f64>; 2],
    /// Inviscid tangential speed at `alfa` (`QINV`, N+NW).
    pub qinv: Vec<f64>,
    /// dQINV/dα (N+NW).
    pub qinv_a: Vec<f64>,
    /// Inviscid Cp (`CPI`; N after SPECAL, N+NW after a viscous solve).
    pub cpi: Vec<f64>,
    /// Viscous tangential speed on airfoil and wake (`QVIS`, N+NW).
    pub qvis: Vec<f64>,
    /// Source strengths `SIG` (N+NW). XFOIL never sets them in analysis; they stay zero.
    pub sig: Vec<f64>,
    /// Source influence matrix `DIJ` ((N+NW)²), once [`Self::qdcalc`] has run.
    pub dij: Option<Matrix>,
}

/// The unfactored GGCALC system: `AIJ`, `BIJ` and the right-hand sides `GAMU`.
#[derive(Debug, Clone, PartialEq)]
pub struct GgcalcSystem {
    /// dψ/dγ with the ψ₀ column and Kutta row ((N+1)²).
    pub aij: Matrix,
    /// −dψ/dσ for the airfoil sources ((N+1) × N).
    pub bij: Matrix,
    /// Right-hand sides for α = 0°, 90° (N+1 each).
    pub gamu: [Vec<f64>; 2],
}

impl GgcalcSystem {
    /// Assembles the system. First part of XFOIL `GGCALC`.
    pub fn assemble(pan: &Paneling) -> Self {
        let n = pan.len();
        let (x, y) = (&pan.nodes.x, &pan.nodes.y);
        let qinf = 1.0;
        let gam = vec![0.0; n];
        let sig = vec![0.0; n];
        let mut gamu = [vec![0.0; n + 1], vec![0.0; n + 1]];
        let mut aij = Matrix::zeros(n + 1, n + 1);
        let mut bij = Matrix::zeros(n + 1, n);
        let mut inf = Influence::new(n, 0);

        // Psi = Psio on the surface; unknowns are gamma_i and Psio
        for i in 0..n {
            {
                let st = Strengths {
                    gam: &gam,
                    gamu: [&gamu[0], &gamu[1]],
                    sig: &sig,
                    alfa: 0.0,
                    qinf,
                };
                let fp = FieldPoint {
                    node: NodeRef::Airfoil(i),
                    x: x[i],
                    y: y[i],
                    nx: pan.nx[i],
                    ny: pan.ny[i],
                };
                psilin(pan, &st, fp, true, &mut inf);
            }
            let res1 = qinf * y[i];
            let res2 = -qinf * x[i];
            for j in 0..n {
                aij[(i, j)] = inf.dzdg[j];
                bij[(i, j)] = -inf.dzdm[j];
            }
            aij[(i, n)] = -1.0;
            gamu[0][i] = -res1;
            gamu[1][i] = -res2;
        }

        // Kutta condition: gamma_1 + gamma_N = 0
        let res = 0.0_f64;
        for j in 0..=n {
            aij[(n, j)] = 0.0;
        }
        aij[(n, 0)] = 1.0;
        aij[(n, n - 1)] = 1.0;
        gamu[0][n] = -res;
        gamu[1][n] = -res;
        for j in 0..n {
            bij[(n, j)] = 0.0;
        }

        if pan.trailing_edge.sharp {
            // zero internal velocity in the TE corner: replace the equation at node N
            let (xp, yp) = (&pan.nodes.xp, &pan.nodes.yp);
            let ag1 = (-yp[0]).atan2(-xp[0]);
            let ag2 = crate::fortran::atanc(yp[n - 1], xp[n - 1], ag1);
            let abis = 0.5 * (ag1 + ag2);
            let cbis = abis.cos();
            let sbis = abis.sin();
            let ds1 = ((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1])).sqrt();
            let ds2 = ((x[n - 1] - x[n - 2]) * (x[n - 1] - x[n - 2])
                + (y[n - 1] - y[n - 2]) * (y[n - 1] - y[n - 2]))
                .sqrt();
            let dsmin = ds1.min(ds2);
            let xbis = pan.te.0 - BWT * dsmin * cbis;
            let ybis = pan.te.1 - BWT * dsmin * sbis;
            let st = Strengths {
                gam: &gam,
                gamu: [&gamu[0], &gamu[1]],
                sig: &sig,
                alfa: 0.0,
                qinf,
            };
            let fp = FieldPoint {
                node: NodeRef::Off,
                x: xbis,
                y: ybis,
                nx: -sbis,
                ny: cbis,
            };
            psilin(pan, &st, fp, true, &mut inf);
            for j in 0..n {
                aij[(n - 1, j)] = inf.dqdg[j];
                bij[(n - 1, j)] = -inf.dqdm[j];
            }
            aij[(n - 1, n)] = 0.0;
            gamu[0][n - 1] = -cbis;
            gamu[1][n - 1] = -sbis;
        }

        Self { aij, bij, gamu }
    }
}

impl InviscidSolution {
    /// Builds and factors the vortex influence system and solves for the α = 0°, 90°
    /// vorticity. Port of XFOIL `GGCALC`.
    pub fn new(pan: &Paneling) -> Result<Self, LinalgError> {
        let n = pan.len();
        let GgcalcSystem { aij, bij, mut gamu } = GgcalcSystem::assemble(pan);
        let aij = ludcmp(aij)?;
        aij.baksub(&mut gamu[0]);
        aij.baksub(&mut gamu[1]);
        let qinvu = [gamu[0][..n].to_vec(), gamu[1][..n].to_vec()];
        let gam = vec![0.0; n];

        Ok(Self {
            aij,
            bij,
            gamu,
            alfa: 0.0,
            gam,
            gam_a: vec![0.0; n],
            psio: 0.0,
            sigte: 0.0,
            gamte: 0.0,
            qinv: vec![0.0; n],
            qinv_a: vec![0.0; n],
            qinvu,
            cpi: vec![0.0; n],
            qvis: vec![0.0; n],
            sig: vec![0.0; n],
            dij: None,
        })
    }

    /// Number of airfoil nodes.
    pub fn n(&self) -> usize {
        self.gam.len()
    }

    /// Current singularity strengths for the kernels.
    pub fn strengths(&self, qinf: f64) -> Strengths<'_> {
        Strengths {
            gam: &self.gam,
            gamu: [&self.gamu[0], &self.gamu[1]],
            sig: &self.sig[..self.n()],
            alfa: self.alfa,
            qinf,
        }
    }

    /// Superimposes the unit solutions for `alfa`, including TE strengths. The first part of
    /// XFOIL `SPECAL` plus `TECALC`.
    fn superimpose(&mut self, pan: &Paneling, alfa: f64) {
        let n = self.n();
        let cosa = alfa.cos();
        let sina = alfa.sin();
        for i in 0..n {
            self.gam[i] = cosa * self.gamu[0][i] + sina * self.gamu[1][i];
            self.gam_a[i] = -sina * self.gamu[0][i] + cosa * self.gamu[1][i];
        }
        self.psio = cosa * self.gamu[0][n] + sina * self.gamu[1][n];
        self.alfa = alfa;

        // TECALC, strength part
        let te = &pan.trailing_edge;
        let (scs, sds) = if te.sharp {
            (1.0, 0.0)
        } else {
            (te.ante / te.dste, te.aste / te.dste)
        };
        self.sigte = 0.5 * (self.gam[0] - self.gam[n - 1]) * scs;
        self.gamte = -0.5 * (self.gam[0] - self.gam[n - 1]) * sds;
    }

    /// Inviscid speed on airfoil and wake at the current α. Port of XFOIL `QISET`.
    pub fn qiset(&mut self) {
        let cosa = self.alfa.cos();
        let sina = self.alfa.sin();
        let m = self.qinvu[0].len();
        self.qinv.resize(m, 0.0);
        self.qinv_a.resize(m, 0.0);
        for i in 0..m {
            self.qinv[i] = cosa * self.qinvu[0][i] + sina * self.qinvu[1][i];
            self.qinv_a[i] = -sina * self.qinvu[0][i] + cosa * self.qinvu[1][i];
        }
    }

    /// Converges the inviscid solution at angle `alfa` (rad), iterating on the lift-dependent
    /// Mach number. Port of XFOIL `SPECAL`.
    pub fn specal(&mut self, pan: &Paneling, alfa: f64, settings: &Settings) -> OperatingPoint {
        let (x, y) = (&pan.nodes.x, &pan.nodes.y);
        let flow = &settings.flow;
        let qinf = 1.0;
        self.superimpose(pan, alfa);
        self.qiset();

        let forces = |minf: f64, s: &Self| {
            clcalc(
                x,
                y,
                &s.gam,
                &s.gam_a,
                alfa,
                minf,
                qinf,
                settings.moment_ref,
            )
        };

        let mut clm = 1.0;
        let mut comp = Compressibility::at_cl(clm, flow, qinf);
        let mut f = forces(comp.minf, self);

        for _ in 0..20 {
            let msq_clm = 2.0 * comp.minf * comp.minf_cl;
            let dclm = (f.cl - clm) / (1.0 - f.cl_msq * msq_clm);
            let clm1 = clm;
            let mut rlx = 1.0;
            // under-relax so that M(CL) stays below 1
            for _ in 0..12 {
                clm = clm1 + rlx * dclm;
                let (minf, minf_clm, reinf, reinf_cl) = mrcl(clm, flow);
                comp.minf = minf;
                comp.minf_cl = minf_clm;
                comp.reinf = reinf;
                comp.reinf_cl = reinf_cl;
                if flow.mach_type == MachType::Fixed || minf == 0.0 || minf_clm != 0.0 {
                    break;
                }
                rlx *= 0.5;
            }
            comp.comset(flow.gamma, qinf);
            f = forces(comp.minf, self);
            if dclm.abs() <= 1.0e-6 {
                break;
            }
        }

        // final Mach, CL and Cp
        let comp = Compressibility::at_cl(f.cl, flow, qinf);
        let f = forces(comp.minf, self);
        let n = self.n();
        self.cpi.resize(n, 0.0);
        cpcalc(&self.qinv[..n], qinf, comp.minf, &mut self.cpi);
        OperatingPoint::inviscid(alfa, qinf, f, comp)
    }

    /// Adds the wake part of the inviscid speed (QWCALC) and resets `QINV` (QISET).
    pub fn set_wake_speeds(&mut self, pan: &Paneling, wake: &Wake) {
        let n = self.n();
        let te = [self.qinvu[0][n - 1], self.qinvu[1][n - 1]];
        let qw = wake.qwcalc(pan, &self.strengths(1.0), te);
        for c in 0..2 {
            self.qinvu[c].truncate(n);
            self.qinvu[c].extend_from_slice(&qw[c]);
        }
        self.sig.resize(n + wake.len(), 0.0);
        self.qvis.resize(n + wake.len(), 0.0);
        self.qiset();
    }

    /// Surface vorticity from the viscous speed. Port of XFOIL `GAMQV`.
    pub fn gamqv(&mut self) {
        let n = self.n();
        self.gam[..n].copy_from_slice(&self.qvis[..n]);
        self.gam_a[..n].copy_from_slice(&self.qinv_a[..n]);
    }

    /// Source influence matrix dQtan/dσ for airfoil and wake. Port of XFOIL `QDCALC`.
    pub fn qdcalc(&mut self, pan: &Paneling, wake: &Wake) {
        let n = self.n();
        let nw = wake.len();
        let (x, y) = (&pan.nodes.x, &pan.nodes.y);
        let mut dij = Matrix::zeros(n + nw, n + nw);
        let mut inf = Influence::new(n, nw);
        let wake_sig = vec![0.0; nw];

        // airfoil sources: back-solve each BIJ column
        let mut col = vec![0.0; n + 1];
        for j in 0..n {
            col.copy_from_slice(self.bij.col(j));
            self.aij.baksub(&mut col);
            dij.col_mut(j)[..n].copy_from_slice(&col[..n]);
        }

        // wake sources: dPsi/dm on the airfoil, Kutta row has no source influence
        let mut bijw = Matrix::zeros(n + 1, nw);
        for i in 0..n {
            let fp = FieldPoint {
                node: NodeRef::Airfoil(i),
                x: x[i],
                y: y[i],
                nx: pan.nx[i],
                ny: pan.ny[i],
            };
            pswlin(wake, &wake_sig, fp, &mut inf);
            for k in 0..nw {
                bijw[(i, k)] = -inf.dzdm[n + k];
            }
        }
        for k in 0..nw {
            bijw[(n, k)] = 0.0;
        }
        if pan.trailing_edge.sharp {
            for k in 0..nw {
                bijw[(n - 1, k)] = 0.0;
            }
        }
        for k in 0..nw {
            self.aij.baksub(bijw.col_mut(k));
        }
        for i in 0..n {
            for k in 0..nw {
                dij[(i, n + k)] = bijw[(i, k)];
            }
        }

        // velocities at wake points: direct influences, plus dQtan/dGam (CIJ) for later
        let mut cij = Matrix::zeros(nw, n);
        {
            let st = self.strengths(1.0);
            for k in 0..nw {
                let i = n + k;
                let fp = FieldPoint {
                    node: NodeRef::Wake(k),
                    x: wake.x[k],
                    y: wake.y[k],
                    nx: wake.nx[k],
                    ny: wake.ny[k],
                };
                psilin(pan, &st, fp, true, &mut inf);
                for j in 0..n {
                    cij[(k, j)] = inf.dqdg[j];
                }
                for j in 0..n {
                    dij[(i, j)] = inf.dqdm[j];
                }
                pswlin(wake, &wake_sig, fp, &mut inf);
                for kk in 0..nw {
                    dij[(i, n + kk)] = inf.dqdm[n + kk];
                }
            }
        }

        // all sources also change the airfoil vorticity, which changes the wake Qtan
        for k in 0..nw {
            let i = n + k;
            for j in 0..n {
                let mut sum = 0.0;
                for kk in 0..n {
                    sum += cij[(k, kk)] * dij[(kk, j)];
                }
                dij[(i, j)] += sum;
            }
            for j in 0..nw {
                let mut sum = 0.0;
                for kk in 0..n {
                    sum += cij[(k, kk)] * bijw[(kk, j)];
                }
                dij[(i, n + j)] += sum;
            }
        }

        // first wake point has the TE velocity
        for j in 0..n + nw {
            dij[(n, j)] = dij[(n - 1, j)];
        }
        self.dij = Some(dij);
    }
}
