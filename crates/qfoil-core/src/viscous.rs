//! Viscous operating point: the Newton iteration of XFOIL/QFoil `VISCAL` and the drag
//! integration `CDCALC` (QFoil D10).

use crate::bl::BoundaryLayer;
use crate::bl::march::MarchInputs;
use crate::bl::station::{Kernel, KernelParams};
use crate::bl::transition::AmplificationModel;
use crate::fortran::{pow, powi};
use crate::inviscid::InviscidSolution;
use crate::inviscid::forces::{Compressibility, clcalc, cpcalc};
use crate::newton::{NewtonSystem, SetblInputs, UpdateInputs, UpdateResult};
use crate::operating::OperatingPoint;
use crate::paneling::Paneling;
use crate::settings::{BlParams, Settings};
use crate::wake::Wake;

/// Viscous-analysis parameters (OPER/VPAR: `ACRIT`, `XSTRIP`, `ITER`, `VACC`, `IDAMP`).
#[derive(Debug, Clone, PartialEq)]
pub struct ViscousSettings {
    /// Critical amplification exponent per side (`N`, `NT`/`NB`).
    pub ncrit: [f64; 2],
    /// Forced-transition x/c per side (`XTR`); ≥ 1 means free transition.
    pub xtrip: [f64; 2],
    /// Newton iteration limit per point (`ITER`). QFoil's default is 20; 100 is used by
    /// the golden reference data.
    pub max_iterations: usize,
    /// BLSOLV drop tolerance (`VACCEL`).
    pub vaccel: f64,
    /// Amplification model (`IDAMP`).
    pub amplification: AmplificationModel,
    /// BL model constants (`/BLPAR/`).
    pub bl: BlParams,
}

impl Default for ViscousSettings {
    fn default() -> Self {
        Self {
            ncrit: [9.0, 9.0],
            xtrip: [1.0, 1.0],
            max_iterations: 100,
            vaccel: 0.01,
            amplification: AmplificationModel::Envelope,
            bl: BlParams::default(),
        }
    }
}

/// VISCAL convergence tolerance on the RMS Newton change (`EPS1`).
pub const EPS1: f64 = 1.0e-4;

/// Stage of the Newton loop at which the observer of [`viscal`] is called.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage {
    /// After SETBL (system assembled).
    Assembled,
    /// After BLSOLV (solution in `vdel`).
    Solved,
    /// After UPDATE.
    Updated,
    /// End of the iteration (after CLCALC/CDCALC).
    Iterated,
}

/// A view of the solver state passed to the [`viscal`] observer.
#[derive(Debug)]
pub struct IterationView<'a> {
    /// 1-based iteration number.
    pub iteration: usize,
    /// Current stage.
    pub stage: Stage,
    /// Boundary layer.
    pub bl: &'a BoundaryLayer,
    /// Operating point (CL, CD, ...).
    pub op: &'a OperatingPoint,
    /// Newton system (assembled or solved).
    pub sys: &'a NewtonSystem,
    /// UPDATE result, from [`Stage::Updated`] on.
    pub update: Option<&'a UpdateResult>,
}

/// Result of [`viscal`].
#[derive(Debug, Clone)]
pub struct ViscalResult {
    /// `RMSBL < EPS1` was reached (`LVCONV`).
    pub converged: bool,
    /// Newton iterations performed.
    pub iterations: usize,
    /// Last UPDATE result.
    pub last: Option<UpdateResult>,
    /// Final boundary layer.
    pub bl: BoundaryLayer,
    /// Viscous Cp on airfoil and wake (`CPV`, N+NW).
    pub cpv: Vec<f64>,
}

/// Total and friction drag. Port of XFOIL/QFoil `CDCALC`.
///
/// QFoil D10: if the wake-end speed ratio `u = Ue/Q∞ ≤ 1`, the Squire–Young momentum
/// thickness is corrected by `θ(1 + (1-u)(u(GWAKE·H1 − 1) − 1))` with Head's
/// `H1 = 3.15 + 1.72/(H−1)`; otherwise the original XFOIL formula is used.
pub fn cdcalc(
    bl: &BoundaryLayer,
    pan: &Paneling,
    alfa: f64,
    qinf: f64,
    tklam: f64,
    gwake: f64,
    bl_initialised: bool,
) -> (f64, f64) {
    let sa = alfa.sin();
    let ca = alfa.cos();
    let cd = if bl_initialised {
        // variables at the end of the wake
        let w = &bl.sides[1];
        let last = bl.nbl[1] - 1;
        let thwake = w.thet[last];
        let urat = w.uedg[last] / qinf;
        let uewake = w.uedg[last] * (1.0 - tklam) / (1.0 - tklam * powi(urat, 2));
        let shwake = w.dstr[last] / w.thet[last];
        // Squire-Young extrapolation to downstream infinity
        let uetinf = uewake / qinf;
        if uetinf <= 1.0 {
            let g_corr = gwake;
            let h1_val = 3.15 + (1.72 / (shwake - 1.0));
            let factor = (1.0 - uetinf) * (uetinf * (g_corr * h1_val - 1.0) - 1.0);
            let dtheta = thwake * factor;
            let th_total = thwake + dtheta;
            2.0 * th_total * pow(uetinf, 0.5 * (5.0 + shwake))
        } else {
            // deep stall / untrustworthy wake: original XFOIL CD
            2.0 * thwake * pow(uetinf, 0.5 * (5.0 + shwake))
        }
    } else {
        0.0
    };

    // friction drag
    let (x, y) = (&pan.nodes.x, &pan.nodes.y);
    let mut cdf = 0.0;
    for is in 0..2 {
        let sd = &bl.sides[is];
        // Fortran IBL = 3..IBLTE
        for ibl in 2..=bl.iblte[is] {
            let i = sd.ipan[ibl];
            let im = sd.ipan[ibl - 1];
            let dx = (x[i] - x[im]) * ca + (y[i] - y[im]) * sa;
            cdf += 0.5 * (sd.tau[ibl] + sd.tau[ibl - 1]) * dx * 2.0 / powi(qinf, 2);
        }
    }
    (cd, cdf)
}

/// Converges a viscous operating point at the current angle of attack.
/// Port of XFOIL/QFoil `VISCAL` for a fresh point (QBlade usage; D9).
///
/// `sol` must hold the inviscid solution at `op.alfa` (from [`InviscidSolution::specal`]).
/// `observe` is called at every [`Stage`] of every iteration; pass `|_| {}` if unused.
pub fn viscal(
    pan: &Paneling,
    sol: &mut InviscidSolution,
    op: &mut OperatingPoint,
    settings: &Settings,
    vs: &ViscousSettings,
    mut observe: impl FnMut(&IterationView<'_>),
) -> ViscalResult {
    let flow = &settings.flow;
    let qinf = op.qinf;

    // wake trajectory, wake inviscid speeds, source influence (XYWAKE, QWCALC, QISET, QDCALC)
    let wake = Wake::trace(pan, &sol.strengths(qinf), settings.wake_length);
    sol.set_wake_speeds(pan, &wake);

    // BL pointers, xi, UINV and QFoil's cold start
    let mut bl = BoundaryLayer::new(pan, &wake, sol);
    sol.qdcalc(pan, &wake);

    let march = MarchInputs {
        pan,
        acrit: vs.ncrit,
        xstrip: vs.xtrip,
    };
    let mut k = Kernel::new(KernelParams::default(), vs.bl.clone(), vs.amplification);
    let mut lblini = false;
    let s = &pan.nodes.s;
    let arc_length = s[s.len() - 1] - s[0];

    let mut converged = false;
    let mut iterations = 0;
    let mut last = None;
    for iter in 1..=vs.max_iterations {
        iterations = iter;
        let dij = sol.dij.take().expect("QDCALC ran");

        // SETBL, BLSOLV, UPDATE
        let inp = SetblInputs {
            march,
            dij: &dij,
            flow,
            clmr: op.forces.cl,
            lalfa: true,
        };
        let mut sys = bl.setbl(&mut k, &mut lblini, &inp);
        observe(&IterationView {
            iteration: iter,
            stage: Stage::Assembled,
            bl: &bl,
            op,
            sys: &sys,
            update: None,
        });
        let ivte1 = bl.sides[0].isys[bl.iblte[0]];
        let ivz = bl.sides[1].isys[bl.iblte[1] + 1];
        sys.blsolv(vs.vaccel, arc_length, ivte1, ivz);
        observe(&IterationView {
            iteration: iter,
            stage: Stage::Solved,
            bl: &bl,
            op,
            sys: &sys,
            update: None,
        });
        let upd = UpdateInputs {
            x: &pan.nodes.x,
            y: &pan.nodes.y,
            dij: &dij,
            minf: op.comp.minf,
            minf_cl: op.comp.minf_cl,
            gamm1: flow.gamma - 1.0,
            lalfa: true,
            clspec: 0.0,
        };
        let r = bl.update(&sys, &upd, &mut op.forces.cl, &mut op.alfa);
        observe(&IterationView {
            iteration: iter,
            stage: Stage::Updated,
            bl: &bl,
            op,
            sys: &sys,
            update: Some(&r),
        });
        sol.dij = Some(dij);

        // new Mach and Re from the new CL (MRCL, COMSET)
        op.comp = Compressibility::at_cl(op.forces.cl, flow, qinf);

        // edge velocities -> panel speed -> vorticity; relocate the stagnation point
        bl.qvfue(&mut sol.qvis);
        sol.gamqv();
        let (qinv, qinv_a) = (sol.qinv.clone(), sol.qinv_a.clone());
        bl.stmove(pan, &wake, &qinv, &qinv_a, &mut sol.gam, &mut sol.qvis);

        // updated CL, CM, CD
        op.forces = clcalc(
            &pan.nodes.x,
            &pan.nodes.y,
            &sol.gam,
            &sol.gam_a,
            op.alfa,
            op.comp.minf,
            qinf,
            settings.moment_ref,
        );
        let (cd, cdf) = cdcalc(&bl, pan, op.alfa, qinf, op.comp.tklam, vs.bl.gwake, lblini);
        op.cd = cd;
        op.cdf = cdf;
        observe(&IterationView {
            iteration: iter,
            stage: Stage::Iterated,
            bl: &bl,
            op,
            sys: &sys,
            update: Some(&r),
        });
        last = Some(r);

        if r.rmsbl < EPS1 {
            converged = true;
            break;
        }
    }

    // final Cp distributions on airfoil and wake
    let m = sol.qinv.len();
    sol.cpi.resize(m, 0.0);
    cpcalc(&sol.qinv, qinf, op.comp.minf, &mut sol.cpi);
    let mut cpv = vec![0.0; m];
    cpcalc(&sol.qvis, qinf, op.comp.minf, &mut cpv);

    ViscalResult {
        converged,
        iterations,
        last,
        bl,
        cpv,
    }
}
