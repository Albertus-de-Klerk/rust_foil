//! Panel influence kernels. Port of XFOIL `PSILIN` and `PSWLIN`.
//!
//! Not ported: the `GEOLIN` branch (geometric sensitivities, used only by inverse design)
//! and the ground-effect image system (`LIMAGE`, off by default and not reachable from
//! QFoil's CLI).

use crate::fortran::PI;
use crate::paneling::Paneling;
use crate::wake::Wake;

/// `1/(2π)` (`HOPI`).
const HOPI: f64 = 0.50 / PI;
/// `1/(4π)` (`QOPI`).
const QOPI: f64 = 0.25 / PI;

/// Identity of the field point, which replaces the overloaded Fortran argument `I`
/// (0 = off-body, 1..N airfoil, N+1..N+NW wake).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeRef {
    /// A point on neither airfoil nor wake (Fortran `I = 0`).
    Off,
    /// Airfoil node, 0-based.
    Airfoil(usize),
    /// Wake node, 0-based (Fortran `N+1+k`).
    Wake(usize),
}

/// Field point for [`psilin`] / [`pswlin`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FieldPoint {
    /// Which node, if any.
    pub node: NodeRef,
    /// Position.
    pub x: f64,
    /// Position.
    pub y: f64,
    /// Direction along which `psi_ni` / `qtan` are taken.
    pub nx: f64,
    /// Direction along which `psi_ni` / `qtan` are taken.
    pub ny: f64,
}

/// Singularity strengths seen by [`psilin`].
#[derive(Debug, Clone, Copy)]
pub struct Strengths<'a> {
    /// Current vorticity `GAM` (N).
    pub gam: &'a [f64],
    /// Unit vorticity for α = 0° and 90° (`GAMU`, at least N rows each).
    pub gamu: [&'a [f64]; 2],
    /// Source strengths `SIG` (airfoil part, N). Always zero in analysis, kept for parity.
    pub sig: &'a [f64],
    /// Angle of attack (rad).
    pub alfa: f64,
    /// Freestream speed.
    pub qinf: f64,
}

/// Results and sensitivity vectors of [`psilin`]/[`pswlin`] (XFOIL `/QMAT/` outputs).
///
/// Reused across calls to avoid allocation; `psilin` writes the airfoil part `[..n]` of
/// `dzdm`/`dqdm`, `pswlin` the wake part `[n..]`.
#[derive(Debug, Clone, PartialEq)]
pub struct Influence {
    /// Stream function ψ.
    pub psi: f64,
    /// dψ/dn at the field point.
    pub psi_ni: f64,
    /// dψ/dγ (N).
    pub dzdg: Vec<f64>,
    /// dψ/dσ (N + NW).
    pub dzdm: Vec<f64>,
    /// dQtan/dγ (N).
    pub dqdg: Vec<f64>,
    /// dQtan/dσ (N + NW).
    pub dqdm: Vec<f64>,
    /// Qtan for α = 0°.
    pub qtan1: f64,
    /// Qtan for α = 90°.
    pub qtan2: f64,
    /// dψ/dQ∞.
    pub z_qinf: f64,
    /// dψ/dα.
    pub z_alfa: f64,
}

impl Influence {
    /// Workspace for `n` airfoil and `nw` wake nodes.
    pub fn new(n: usize, nw: usize) -> Self {
        Self {
            psi: 0.0,
            psi_ni: 0.0,
            dzdg: vec![0.0; n],
            dzdm: vec![0.0; n + nw],
            dqdg: vec![0.0; n],
            dqdm: vec![0.0; n + nw],
            qtan1: 0.0,
            qtan2: 0.0,
            z_qinf: 0.0,
            z_alfa: 0.0,
        }
    }
}

/// Stream function at a field point due to the freestream and the airfoil vortex (and,
/// with `siglin`, source) panels, plus sensitivities. Port of XFOIL `PSILIN`.
pub fn psilin(
    pan: &Paneling,
    st: &Strengths<'_>,
    fp: FieldPoint,
    siglin: bool,
    out: &mut Influence,
) {
    let (x, y, s) = (&pan.nodes.x, &pan.nodes.y, &pan.nodes.s);
    let n = x.len();
    let (xi, yi, nxi, nyi) = (fp.x, fp.y, fp.nx, fp.ny);

    let seps = (s[n - 1] - s[0]) * 1.0e-5;
    let cosa = st.alfa.cos();
    let sina = st.alfa.sin();

    out.dzdg.fill(0.0);
    out.dqdg.fill(0.0);
    out.dzdm[..n].fill(0.0);
    out.dqdm[..n].fill(0.0);
    out.z_qinf = 0.0;
    out.z_alfa = 0.0;
    let mut psi = 0.0;
    let mut psi_ni = 0.0;
    let mut qtan1 = 0.0;
    let mut qtan2 = 0.0;

    let te = &pan.trailing_edge;
    let (scs, sds) = if te.sharp {
        (1.0, 0.0)
    } else {
        (te.ante / te.dste, te.aste / te.dste)
    };

    let on_airfoil = matches!(fp.node, NodeRef::Airfoil(_));
    let is_node = |j: usize| fp.node == NodeRef::Airfoil(j);

    // Values of the last panel visited; the TE panel (jo = n-1) reuses them after the loop.
    let (mut jo_te, mut jp_te) = (0, 0);
    let (mut x1, mut x2, mut yy, mut g1, mut g2, mut t1, mut t2, mut apan) =
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut x1i, mut x2i, mut yyi) = (0.0, 0.0, 0.0);
    let mut reached_te = false;

    // Fortran JO = 1..N  ->  jo = 0..n-1
    for jo in 0..n {
        let mut jp = jo + 1;
        let mut jm = jo.wrapping_sub(1);
        let mut jq = jp + 1;
        if jo == 0 {
            jm = jo;
        } else if jo == n - 2 {
            jq = jp;
        } else if jo == n - 1 {
            jp = 0;
            let (dx, dy) = (x[jo] - x[jp], y[jo] - y[jp]);
            if dx * dx + dy * dy < seps * seps {
                // closed TE: no TE panel (GO TO 12)
                break;
            }
        }

        let dso = ((x[jo] - x[jp]) * (x[jo] - x[jp]) + (y[jo] - y[jp]) * (y[jo] - y[jp])).sqrt();
        if dso == 0.0 {
            continue; // null panel
        }
        let dsio = 1.0 / dso;
        apan = pan.apanel[jo];

        let rx1 = xi - x[jo];
        let ry1 = yi - y[jo];
        let rx2 = xi - x[jp];
        let ry2 = yi - y[jp];
        let sx = (x[jp] - x[jo]) * dsio;
        let sy = (y[jp] - y[jo]) * dsio;
        x1 = sx * rx1 + sy * ry1;
        x2 = sx * rx2 + sy * ry2;
        yy = sx * ry1 - sy * rx1;
        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        // reflection flag keeps atan2 within ±π/2 off the surface
        let sgn = if on_airfoil {
            1.0
        } else {
            1.0_f64.copysign(yy)
        };

        if !is_node(jo) && rs1 > 0.0 {
            g1 = rs1.ln();
            t1 = (sgn * x1).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI;
        } else {
            g1 = 0.0;
            t1 = 0.0;
        }
        if !is_node(jp) && rs2 > 0.0 {
            g2 = rs2.ln();
            t2 = (sgn * x2).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI;
        } else {
            g2 = 0.0;
            t2 = 0.0;
        }

        x1i = sx * nxi + sy * nyi;
        x2i = sx * nxi + sy * nyi;
        yyi = sx * nyi - sy * nxi;

        if jo == n - 1 {
            jo_te = jo;
            jp_te = jp;
            reached_te = true;
            break; // GO TO 11: TE panel
        }

        if siglin {
            let sig = st.sig;
            // midpoint quantities
            let x0 = 0.5 * (x1 + x2);
            let rs0 = x0 * x0 + yy * yy;
            let g0 = rs0.ln();
            let t0 = (sgn * x0).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI;

            // 1-0 half-panel
            let dxinv = 1.0 / (x1 - x0);
            let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
            let pdif =
                ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv;
            let psx1 = -(t1 - apan);
            let psx0 = t0 - apan;
            let psyy = 0.5 * (g1 - g0);
            let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
            let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
            let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

            let dsm =
                ((x[jp] - x[jm]) * (x[jp] - x[jm]) + (y[jp] - y[jm]) * (y[jp] - y[jm])).sqrt();
            let dsim = 1.0 / dsm;
            let ssum = (sig[jp] - sig[jo]) * dsio + (sig[jp] - sig[jm]) * dsim;
            let sdif = (sig[jp] - sig[jo]) * dsio - (sig[jp] - sig[jm]) * dsim;
            psi += QOPI * (psum * ssum + pdif * sdif);

            out.dzdm[jm] += QOPI * (-psum * dsim + pdif * dsim);
            out.dzdm[jo] += QOPI * (-psum * dsio - pdif * dsio);
            out.dzdm[jp] += QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim));

            let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
            let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
            psi_ni += QOPI * (psni * ssum + pdni * sdif);
            out.dqdm[jm] += QOPI * (-psni * dsim + pdni * dsim);
            out.dqdm[jo] += QOPI * (-psni * dsio - pdni * dsio);
            out.dqdm[jp] += QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim));

            // 0-2 half-panel
            let dxinv = 1.0 / (x0 - x2);
            let psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
            let pdif =
                ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv;
            let psx0 = -(t0 - apan);
            let psx2 = t2 - apan;
            let psyy = 0.5 * (g0 - g2);
            let pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
            let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
            let pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

            let dsp =
                ((x[jq] - x[jo]) * (x[jq] - x[jo]) + (y[jq] - y[jo]) * (y[jq] - y[jo])).sqrt();
            let dsip = 1.0 / dsp;
            let ssum = (sig[jq] - sig[jo]) * dsip + (sig[jp] - sig[jo]) * dsio;
            let sdif = (sig[jq] - sig[jo]) * dsip - (sig[jp] - sig[jo]) * dsio;
            psi += QOPI * (psum * ssum + pdif * sdif);

            out.dzdm[jo] += QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
            out.dzdm[jp] += QOPI * (psum * dsio - pdif * dsio);
            out.dzdm[jq] += QOPI * (psum * dsip + pdif * dsip);

            let psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
            let pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
            psi_ni += QOPI * (psni * ssum + pdni * sdif);
            out.dqdm[jo] += QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
            out.dqdm[jp] += QOPI * (psni * dsio - pdni * dsio);
            out.dqdm[jq] += QOPI * (psni * dsip + pdni * dsip);
        }

        // vortex panel contribution
        let dxinv = 1.0 / (x1 - x2);
        let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
        let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;
        let psx1 = 0.5 * g1;
        let psx2 = -0.5 * g2;
        let psyy = t1 - t2;
        let pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv;
        let pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv;
        let pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv;

        let gsum1 = st.gamu[0][jp] + st.gamu[0][jo];
        let gsum2 = st.gamu[1][jp] + st.gamu[1][jo];
        let gdif1 = st.gamu[0][jp] - st.gamu[0][jo];
        let gdif2 = st.gamu[1][jp] - st.gamu[1][jo];
        let gsum = st.gam[jp] + st.gam[jo];
        let gdif = st.gam[jp] - st.gam[jo];

        psi += QOPI * (psis * gsum + psid * gdif);
        out.dzdg[jo] += QOPI * (psis - psid);
        out.dzdg[jp] += QOPI * (psis + psid);

        let psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
        let pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
        psi_ni += QOPI * (gsum * psni + gdif * pdni);
        qtan1 += QOPI * (gsum1 * psni + gdif1 * pdni);
        qtan2 += QOPI * (gsum2 * psni + gdif2 * pdni);
        out.dqdg[jo] += QOPI * (psni - pdni);
        out.dqdg[jp] += QOPI * (psni + pdni);
    }

    if reached_te {
        // TE panel (label 11): uniform source and vortex of strengths set by the TE gap
        let (jo, jp) = (jo_te, jp_te);
        let psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
        let pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
        let psigx1 = -(t1 - apan);
        let psigx2 = t2 - apan;
        let psigyy = 0.5 * (g1 - g2);
        let pgamx1 = 0.5 * g1;
        let pgamx2 = -0.5 * g2;
        let pgamyy = t1 - t2;
        let psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
        let pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

        let sigte1 = 0.5 * scs * (st.gamu[0][jp] - st.gamu[0][jo]);
        let sigte2 = 0.5 * scs * (st.gamu[1][jp] - st.gamu[1][jo]);
        let gamte1 = -0.5 * sds * (st.gamu[0][jp] - st.gamu[0][jo]);
        let gamte2 = -0.5 * sds * (st.gamu[1][jp] - st.gamu[1][jo]);
        let sigte = 0.5 * scs * (st.gam[jp] - st.gam[jo]);
        let gamte = -0.5 * sds * (st.gam[jp] - st.gam[jo]);

        psi += HOPI * (psig * sigte + pgam * gamte);
        out.dzdg[jo] -= HOPI * psig * scs * 0.5;
        out.dzdg[jp] += HOPI * psig * scs * 0.5;
        out.dzdg[jo] += HOPI * pgam * sds * 0.5;
        out.dzdg[jp] -= HOPI * pgam * sds * 0.5;

        psi_ni += HOPI * (psigni * sigte + pgamni * gamte);
        qtan1 += HOPI * (psigni * sigte1 + pgamni * gamte1);
        qtan2 += HOPI * (psigni * sigte2 + pgamni * gamte2);
        out.dqdg[jo] -= HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
        out.dqdg[jp] += HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
    }

    // freestream (label 12)
    let qinf = st.qinf;
    psi += qinf * (cosa * yi - sina * xi);
    psi_ni += qinf * (cosa * nyi - sina * nxi);
    qtan1 += qinf * nyi;
    qtan2 -= qinf * nxi;
    out.z_qinf += cosa * yi - sina * xi;
    out.z_alfa -= qinf * (sina * yi + cosa * xi);

    out.psi = psi;
    out.psi_ni = psi_ni;
    out.qtan1 = qtan1;
    out.qtan2 = qtan2;
}

/// Stream function and its sensitivities due to the wake source panels.
/// Port of XFOIL `PSWLIN`. Writes `dzdm[n..]`, `dqdm[n..]`, `psi`, `psi_ni`.
pub fn pswlin(wake: &Wake, wake_sig: &[f64], fp: FieldPoint, out: &mut Influence) {
    let (x, y) = (&wake.x, &wake.y);
    let nw = x.len();
    let n = out.dzdg.len();
    let (xi, yi, nxi, nyi) = (fp.x, fp.y, fp.nx, fp.ny);
    let (dzdm, dqdm) = (&mut out.dzdm[n..], &mut out.dqdm[n..]);
    dzdm.fill(0.0);
    dqdm.fill(0.0);
    let mut psi = 0.0;
    let mut psi_ni = 0.0;
    let sig = wake_sig;
    let on_wake = matches!(fp.node, NodeRef::Wake(_));
    let is_node = |j: usize| fp.node == NodeRef::Wake(j);

    // Fortran JO = N+1 .. N+NW-1  ->  wake-local jo = 0..nw-2
    for jo in 0..nw - 1 {
        let jp = jo + 1;
        let mut jm = jo.wrapping_sub(1);
        let mut jq = jp + 1;
        if jo == 0 {
            jm = jo;
        } else if jo == nw - 2 {
            jq = jp;
        }
        let dso = ((x[jo] - x[jp]) * (x[jo] - x[jp]) + (y[jo] - y[jp]) * (y[jo] - y[jp])).sqrt();
        let dsio = 1.0 / dso;
        let apan = wake.apanel[jo];

        let rx1 = xi - x[jo];
        let ry1 = yi - y[jo];
        let rx2 = xi - x[jp];
        let ry2 = yi - y[jp];
        let sx = (x[jp] - x[jo]) * dsio;
        let sy = (y[jp] - y[jo]) * dsio;
        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;
        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        let sgn = if on_wake { 1.0 } else { 1.0_f64.copysign(yy) };
        // note the minus sign on the reflection term, unlike PSILIN
        let (g1, t1) = if !is_node(jo) && rs1 > 0.0 {
            (
                rs1.ln(),
                (sgn * x1).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI,
            )
        } else {
            (0.0, 0.0)
        };
        let (g2, t2) = if !is_node(jp) && rs2 > 0.0 {
            (
                rs2.ln(),
                (sgn * x2).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI,
            )
        } else {
            (0.0, 0.0)
        };
        let x1i = sx * nxi + sy * nyi;
        let x2i = sx * nxi + sy * nyi;
        let yyi = sx * nyi - sy * nxi;

        let x0 = 0.5 * (x1 + x2);
        let rs0 = x0 * x0 + yy * yy;
        let g0 = rs0.ln();
        let t0 = (sgn * x0).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI;

        // 1-0 half-panel
        let dxinv = 1.0 / (x1 - x0);
        let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
        let pdif =
            ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv;
        let psx1 = -(t1 - apan);
        let psx0 = t0 - apan;
        let psyy = 0.5 * (g1 - g0);
        let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
        let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
        let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

        let dsm = ((x[jp] - x[jm]) * (x[jp] - x[jm]) + (y[jp] - y[jm]) * (y[jp] - y[jm])).sqrt();
        let dsim = 1.0 / dsm;
        let ssum = (sig[jp] - sig[jo]) * dsio + (sig[jp] - sig[jm]) * dsim;
        let sdif = (sig[jp] - sig[jo]) * dsio - (sig[jp] - sig[jm]) * dsim;
        psi += QOPI * (psum * ssum + pdif * sdif);
        dzdm[jm] += QOPI * (-psum * dsim + pdif * dsim);
        dzdm[jo] += QOPI * (-psum * dsio - pdif * dsio);
        dzdm[jp] += QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim));

        let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
        let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
        psi_ni += QOPI * (psni * ssum + pdni * sdif);
        dqdm[jm] += QOPI * (-psni * dsim + pdni * dsim);
        dqdm[jo] += QOPI * (-psni * dsio - pdni * dsio);
        dqdm[jp] += QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim));

        // 0-2 half-panel
        let dxinv = 1.0 / (x0 - x2);
        let psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
        let pdif =
            ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv;
        let psx0 = -(t0 - apan);
        let psx2 = t2 - apan;
        let psyy = 0.5 * (g0 - g2);
        let pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
        let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
        let pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

        let dsp = ((x[jq] - x[jo]) * (x[jq] - x[jo]) + (y[jq] - y[jo]) * (y[jq] - y[jo])).sqrt();
        let dsip = 1.0 / dsp;
        let ssum = (sig[jq] - sig[jo]) * dsip + (sig[jp] - sig[jo]) * dsio;
        let sdif = (sig[jq] - sig[jo]) * dsip - (sig[jp] - sig[jo]) * dsio;
        psi += QOPI * (psum * ssum + pdif * sdif);
        dzdm[jo] += QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
        dzdm[jp] += QOPI * (psum * dsio - pdif * dsio);
        dzdm[jq] += QOPI * (psum * dsip + pdif * dsip);

        let psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
        let pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
        psi_ni += QOPI * (psni * ssum + pdni * sdif);
        dqdm[jo] += QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
        dqdm[jp] += QOPI * (psni * dsio - pdni * dsio);
        dqdm[jq] += QOPI * (psni * dsip + pdni * dsip);
    }
    out.psi = psi;
    out.psi_ni = psi_ni;
    let _ = n;
    let _ = nw;
}
