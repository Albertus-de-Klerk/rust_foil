//! Airfoil panel nodes. Port of XFOIL `PANGEN`, `ABCOPY`, `TECALC` (geometry part),
//! `NCALC` and `APCALC`.

use crate::error::GeometryError;
use crate::fortran::PI;
use crate::geometry::{AirfoilSource, Geometry, lefind};
use crate::limits::IQX;
use crate::spline::{Curve, deval, segspl, seval, trisol};

/// Curvature-based panelling parameters (XFOIL `/CR12/`, `NPAN`; defaults from `INIT`).
#[derive(Debug, Clone, PartialEq)]
pub struct PanelingSettings {
    /// Number of panel nodes (`NPAN`).
    pub npan: usize,
    /// Curvature attraction (`CVPAR`): 0 uniform, ~1 strongly bunched.
    pub cvpar: f64,
    /// TE / LE panel density ratio (`CTERAT`).
    pub cterat: f64,
    /// Refined-area / LE panel density ratio (`CTRRAT`).
    pub ctrrat: f64,
    /// Top-side refinement x/c limits (`XSREF1`, `XSREF2`).
    pub xsref: (f64, f64),
    /// Bottom-side refinement x/c limits (`XPREF1`, `XPREF2`).
    pub xpref: (f64, f64),
}

impl Default for PanelingSettings {
    fn default() -> Self {
        Self {
            npan: 160,
            cvpar: 1.0,
            cterat: 0.15,
            ctrrat: 0.2,
            xsref: (1.0, 1.0),
            xpref: (1.0, 1.0),
        }
    }
}

/// How the current (panel) airfoil is derived from the buffer airfoil.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum PanelingMode {
    /// QFoil's behaviour for the input: files use their nodes (ABCOPY), NACA
    /// sections are re-panelled with default settings (PANGEN).
    #[default]
    Auto,
    /// Use the buffer points as panel nodes (`ABCOPY`, what `LOAD` does).
    InputNodes,
    /// Re-panel by curvature (`PANE` / `PANGEN`).
    Generate(PanelingSettings),
}

/// Trailing-edge geometry (`TECALC`, geometric part).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrailingEdge {
    /// TE gap area projected normal to the TE bisector (`ANTE`).
    pub ante: f64,
    /// TE gap area projected along the bisector (`ASTE`).
    pub aste: f64,
    /// TE gap length (`DSTE`).
    pub dste: f64,
    /// `DSTE < 1e-4 * CHORD` (`SHARP`).
    pub sharp: bool,
}

/// The current airfoil: panel nodes and derived geometry (COMMON `/CR05/`, `/CR06/` part).
#[derive(Debug, Clone, PartialEq)]
pub struct Paneling {
    /// Nodes `X, Y`, arc length `S`, spline derivatives `XP, YP`.
    pub nodes: Curve,
    /// Unit normal x components (`NX`).
    pub nx: Vec<f64>,
    /// Unit normal y components (`NY`).
    pub ny: Vec<f64>,
    /// Panel angles (`APANEL`); element `n-1` is the TE panel.
    pub apanel: Vec<f64>,
    /// LE arc length (`SLE`).
    pub sle: f64,
    /// LE point.
    pub le: (f64, f64),
    /// TE midpoint.
    pub te: (f64, f64),
    /// LE–TE distance (`CHORD`).
    pub chord: f64,
    /// TE gap geometry.
    pub trailing_edge: TrailingEdge,
}

impl Paneling {
    /// Builds panels according to `mode`.
    pub fn new(geom: &Geometry, mode: &PanelingMode) -> Result<Self, GeometryError> {
        match mode {
            PanelingMode::Auto => match geom.source {
                AirfoilSource::File => Self::from_buffer(geom),
                AirfoilSource::Naca => Self::generate(geom, &PanelingSettings::default()),
            },
            PanelingMode::InputNodes => Self::from_buffer(geom),
            PanelingMode::Generate(settings) => Self::generate(geom, settings),
        }
    }

    /// Number of panel nodes `N`.
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Whether there are no nodes.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Uses the buffer points as panel nodes. Port of XFOIL `ABCOPY`.
    ///
    /// Successive identical points are removed pairwise, as in the Fortran (a run of
    /// three identical points keeps two).
    pub fn from_buffer(geom: &Geometry) -> Result<Self, GeometryError> {
        let b = &geom.buffer;
        let max = IQX - 5;
        if b.len() > max {
            return Err(GeometryError::TooManyPoints {
                what: "panel nodes (ABCOPY)",
                got: b.len(),
                max,
            });
        }
        let (mut x, mut y) = (b.x.clone(), b.y.clone());
        // Fortran: I=1; loop { I=I+1; if doubled, shift down and N=N-1 ; until I >= N }.
        // After a removal the index still advances, so only pairs are stripped.
        let mut i = 0;
        loop {
            i += 1;
            if x[i - 1] == x[i] && y[i - 1] == y[i] {
                x.remove(i);
                y.remove(i);
            }
            if i + 1 >= x.len() {
                break;
            }
        }
        Self::finish(Curve::new(x, y)?)
    }

    /// Curvature-based node distribution. Port of XFOIL `PANGEN`.
    pub fn generate(geom: &Geometry, set: &PanelingSettings) -> Result<Self, GeometryError> {
        let b = &geom.buffer;
        let (xb, yb, sb, xbp, ybp) = (&b.x, &b.y, &b.s, &b.xp, &b.yp);
        let nb = b.len();
        let n = set.npan;
        const IPFAC: usize = 5;
        let nn = IPFAC * (n - 1) + 1;
        let len = nb.max(nn);
        let (mut w1, mut w2, mut w3, mut w4) = (
            vec![0.0; len],
            vec![0.0; len],
            vec![0.0; len],
            vec![0.0; len],
        );
        let (mut w5, mut w6) = (vec![0.0; len], vec![0.0; len]);

        // normalising length (~ chord) and curvature array
        let sbref = 0.5 * (sb[nb - 1] - sb[0]);
        for i in 0..nb {
            w5[i] = b.curvature(sb[i]).abs() * sbref;
        }

        let sble = lefind(xb, xbp, yb, ybp, sb);
        let cvle = b.curvature(sble).abs() * sbref;

        // Doubled point (sharp corner) at the LE: `ible` is 0-based (Fortran IBLE = ible+1).
        let ible = (0..nb - 1).find(|&i| sble == sb[i] && sble == sb[i + 1]);

        let xble = seval(sble, xb, xbp, sb);
        let yble = seval(sble, yb, ybp, sb);
        let xbte = 0.5 * (xb[0] + xb[nb - 1]);
        let ybte = 0.5 * (yb[0] + yb[nb - 1]);
        let chbsq = (xbte - xble) * (xbte - xble) + (ybte - yble) * (ybte - yble);

        // average curvature over 2*NK+1 points near the LE
        const NK: i32 = 3;
        let mut cvsum = 0.0;
        for k in -NK..=NK {
            let frac = f64::from(k) / f64::from(NK);
            let sbk = sble + frac * sbref / cvle.max(20.0);
            cvsum += b.curvature(sbk).abs() * sbref;
        }
        let mut cvavg = cvsum / f64::from(2 * NK + 1);
        if ible.is_some() {
            cvavg = 10.0;
        }

        let cc = 6.0 * set.cvpar;
        let cvte = cvavg * set.cterat;
        w5[0] = cvte;
        w5[nb - 1] = cvte;

        // smoothing length: 1 / averaged LE curvature, bounded (NPAN/2 is integer division)
        let smool = (1.0 / cvavg.max(20.0)).max(0.25 / (set.npan / 2) as f64);
        let smoosq = (smool * sbref) * (smool * sbref);

        // tridiagonal system for the smoothed curvature
        w2[0] = 1.0;
        w3[0] = 0.0;
        for i in 1..nb - 1 {
            let dsm = sb[i] - sb[i - 1];
            let dsp = sb[i + 1] - sb[i];
            let dso = 0.5 * (sb[i + 1] - sb[i - 1]);
            if dsm == 0.0 || dsp == 0.0 {
                // leave curvature at a corner unchanged
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
            } else {
                w1[i] = smoosq * (-1.0 / dsm) / dso;
                w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
                w3[i] = smoosq * (-1.0 / dsp) / dso;
            }
        }
        w1[nb - 1] = 0.0;
        w2[nb - 1] = 1.0;

        // fix curvature at the LE by modifying the equations next to it
        for i in 1..nb - 1 {
            let at_le = sb[i] == sble || ible.is_some_and(|ib| i == ib || i == ib + 1);
            if at_le {
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w5[i] = cvle;
            } else if sb[i - 1] < sble && sb[i] > sble {
                // Fortran reads SB(I-2); with I=2 that is outside the array.
                if i < 2 {
                    return Err(GeometryError::LeadingEdgeAtEnd);
                }
                // equation at the node just before the LE point
                let dsm = sb[i - 1] - sb[i - 2];
                let dsp = sble - sb[i - 1];
                let dso = 0.5 * (sble - sb[i - 2]);
                w1[i - 1] = smoosq * (-1.0 / dsm) / dso;
                w2[i - 1] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
                w3[i - 1] = 0.0;
                w5[i - 1] += smoosq * cvle / (dsp * dso);
                // equation at the node just after the LE point
                let dsm = sb[i] - sble;
                let dsp = sb[i + 1] - sb[i];
                let dso = 0.5 * (sb[i + 1] - sble);
                w1[i] = 0.0;
                w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
                w3[i] = smoosq * (-1.0 / dsp) / dso;
                w5[i] += smoosq * cvle / (dsm * dso);
                break;
            }
        }

        // artificial curvature at refinement points
        for i in 1..nb - 1 {
            let xoc = ((xb[i] - xble) * (xbte - xble) + (yb[i] - yble) * (ybte - yble)) / chbsq;
            let (lo, hi) = if sb[i] < sble { set.xsref } else { set.xpref };
            if xoc > lo && xoc < hi {
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w5[i] = cvle * set.ctrrat;
            }
        }

        // solve for the smoothed curvature (separately either side of a sharp LE)
        match ible {
            None => trisol(&mut w2[..nb], &w1[..nb], &mut w3[..nb], &mut w5[..nb]),
            Some(ib) => {
                // Fortran: TRISOL(..., IBLE) then TRISOL(W(IBLE+1), ..., NB-IBLE)
                let k = ib + 1;
                trisol(&mut w2[..k], &w1[..k], &mut w3[..k], &mut w5[..k]);
                trisol(&mut w2[k..nb], &w1[k..nb], &mut w3[k..nb], &mut w5[k..nb]);
            }
        }

        let cvmax = w5[..nb].iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        for v in &mut w5[..nb] {
            *v /= cvmax;
        }
        segspl(&w5[..nb], &mut w6[..nb], sb)?;
        let (cv, cvp) = (&w5[..nb], &w6[..nb]);

        // initial node positions, uniform in s, with IPFAC times more nodes than requested
        let rdste = 0.667;
        let rtf = (rdste - 1.0) * 2.0 + 1.0;
        let mut snew = vec![0.0; nn];
        let mut nn1 = 0; // Fortran NN1 (1-based count of nodes up to a sharp LE)
        match ible {
            None => {
                let dsavg = (sb[nb - 1] - sb[0]) / ((nn - 3) as f64 + 2.0 * rtf);
                snew[0] = sb[0];
                for i in 1..nn - 1 {
                    // Fortran I = i+1: SB(1) + DSAVG*(FLOAT(I-2) + RTF)
                    snew[i] = sb[0] + dsavg * ((i - 1) as f64 + rtf);
                }
                snew[nn - 1] = sb[nb - 1];
            }
            Some(ib) => {
                let nfrac1 = (n * (ib + 1)) / nb;
                nn1 = IPFAC * (nfrac1 - 1) + 1;
                let dsavg1 = (sble - sb[0]) / ((nn1 - 2) as f64 + rtf);
                snew[0] = sb[0];
                for i in 1..nn1 {
                    snew[i] = sb[0] + dsavg1 * ((i - 1) as f64 + rtf);
                }
                let nn2 = nn - nn1 + 1;
                let dsavg2 = (sb[nb - 1] - sble) / ((nn2 - 2) as f64 + rtf);
                // Fortran: DO I=2, NN2-1 ; SNEW(I-1+NN1) = SBLE + DSAVG2*(FLOAT(I-2)+RTF)
                for i in 2..nn2 {
                    snew[i - 2 + nn1] = sble + dsavg2 * ((i - 2) as f64 + rtf);
                }
                snew[nn - 1] = sb[nb - 1];
            }
        }

        // Newton iteration for node positions:
        // (1 + C*curvature)*ds is made equal on both sides of every node
        for _ in 0..20 {
            let mut cv1 = seval(snew[0], cv, cvp, sb);
            let mut cv2 = seval(snew[1], cv, cvp, sb);
            let mut cvs1 = deval(snew[0], cv, cvp, sb);
            let mut cvs2 = deval(snew[1], cv, cvp, sb);
            let mut cavm = (cv1 * cv1 + cv2 * cv2).sqrt();
            let (mut cavm_s1, mut cavm_s2) = if cavm == 0.0 {
                (0.0, 0.0)
            } else {
                (cvs1 * cv1 / cavm, cvs2 * cv2 / cavm)
            };

            for i in 1..nn - 1 {
                let dsm = snew[i] - snew[i - 1];
                let dsp = snew[i] - snew[i + 1];
                let cv3 = seval(snew[i + 1], cv, cvp, sb);
                let cvs3 = deval(snew[i + 1], cv, cvp, sb);
                let cavp = (cv3 * cv3 + cv2 * cv2).sqrt();
                let (cavp_s2, cavp_s3) = if cavp == 0.0 {
                    (0.0, 0.0)
                } else {
                    (cvs2 * cv2 / cavp, cvs3 * cv3 / cavp)
                };
                let fm = cc * cavm + 1.0;
                let fp = cc * cavp + 1.0;
                let rez = dsp * fp + dsm * fm;
                w1[i] = -fm + cc * dsm * cavm_s1;
                w2[i] = fp + fm + cc * (dsp * cavp_s2 + dsm * cavm_s2);
                w3[i] = -fp + cc * dsp * cavp_s3;
                w4[i] = -rez;
                cv1 = cv2;
                cv2 = cv3;
                cvs1 = cvs2;
                cvs2 = cvs3;
                cavm = cavp;
                cavm_s1 = cavp_s2;
                cavm_s2 = cavp_s3;
            }
            let _ = (cv1, cvs1); // rolled like the Fortran; not read after the loop

            // fix endpoints (TE)
            w2[0] = 1.0;
            w3[0] = 0.0;
            w4[0] = 0.0;
            w1[nn - 1] = 0.0;
            w2[nn - 1] = 1.0;
            w4[nn - 1] = 0.0;

            if rtf != 1.0 {
                // TE panel length ratio RTF on the panels next to the TE
                let i = 1;
                w4[i] = -((snew[i] - snew[i - 1]) + rtf * (snew[i] - snew[i + 1]));
                w1[i] = -1.0;
                w2[i] = 1.0 + rtf;
                w3[i] = -rtf;
                let i = nn - 2;
                w4[i] = -((snew[i] - snew[i + 1]) + rtf * (snew[i] - snew[i - 1]));
                w3[i] = -1.0;
                w2[i] = 1.0 + rtf;
                w1[i] = -rtf;
            }
            if ible.is_some() {
                // pin the sharp LE node (Fortran I = NN1)
                let i = nn1 - 1;
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w4[i] = sble - snew[i];
            }

            trisol(&mut w2[..nn], &w1[..nn], &mut w3[..nn], &mut w4[..nn]);

            // under-relax to keep nodes in order
            let mut rlx = 1.0;
            let mut dmax = 0.0_f64;
            for i in 0..nn - 1 {
                let ds = snew[i + 1] - snew[i];
                let dds = w4[i + 1] - w4[i];
                let dsrat = 1.0 + rlx * dds / ds;
                if dsrat > 4.0 {
                    rlx = (4.0 - 1.0) * ds / dds;
                }
                if dsrat < 0.2 {
                    rlx = (0.2 - 1.0) * ds / dds;
                }
                dmax = w4[i].abs().max(dmax);
            }
            for i in 1..nn - 1 {
                snew[i] += rlx * w4[i];
            }
            if dmax.abs() < 1.0e-3 {
                break;
            }
        }

        // panel nodes: every IPFAC-th temporary node
        let mut s: Vec<f64> = (0..n).map(|i| snew[IPFAC * i]).collect();
        let mut x: Vec<f64> = s.iter().map(|&v| seval(v, xb, xbp, sb)).collect();
        let mut y: Vec<f64> = s.iter().map(|&v| seval(v, yb, ybp, sb)).collect();

        // insert a node at every buffer corner (doubled point)
        for ib in 0..nb - 1 {
            if sb[ib] != sb[ib + 1] {
                continue;
            }
            let (xbc, ybc, sbc) = (xb[ib], yb[ib], sb[ib]);
            let Some(i) = s.iter().position(|&v| v > sbc) else {
                continue;
            };
            x.insert(i, xbc);
            y.insert(i, ybc);
            s.insert(i, sbc);
            if x.len() > IQX - 1 {
                return Err(GeometryError::TooManyPoints {
                    what: "panel nodes (PANGEN)",
                    got: x.len(),
                    max: IQX - 1,
                });
            }
            // shift neighbours to keep panel sizes comparable (Fortran I-2 >= 1, I+2 <= N)
            if i >= 2 {
                s[i - 1] = 0.5 * (s[i] + s[i - 2]);
                x[i - 1] = seval(s[i - 1], xb, xbp, sb);
                y[i - 1] = seval(s[i - 1], yb, ybp, sb);
            }
            if i + 2 < x.len() {
                s[i + 1] = 0.5 * (s[i] + s[i + 2]);
                x[i + 1] = seval(s[i + 1], xb, xbp, sb);
                y[i + 1] = seval(s[i + 1], yb, ybp, sb);
            }
        }

        Self::finish(Curve::new(x, y)?)
    }

    /// Common tail of PANGEN and ABCOPY: LE/TE, chord, TECALC, NCALC, APCALC.
    fn finish(nodes: Curve) -> Result<Self, GeometryError> {
        let n = nodes.len();
        let sle = nodes.leading_edge();
        let le = (nodes.x_at(sle), nodes.y_at(sle));
        let te = (
            0.5 * (nodes.x[0] + nodes.x[n - 1]),
            0.5 * (nodes.y[0] + nodes.y[n - 1]),
        );
        let chord = ((te.0 - le.0) * (te.0 - le.0) + (te.1 - le.1) * (te.1 - le.1)).sqrt();
        let trailing_edge = tecalc(&nodes, chord);
        let (nx, ny) = ncalc(&nodes);
        let apanel = apcalc(&nodes, &nx, &ny, trailing_edge.sharp);
        Ok(Self {
            nodes,
            nx,
            ny,
            apanel,
            sle,
            le,
            te,
            chord,
            trailing_edge,
        })
    }
}

/// TE gap geometry. Geometric part of XFOIL `TECALC` (the TE panel strengths belong to
/// the inviscid solution).
fn tecalc(p: &Curve, chord: f64) -> TrailingEdge {
    let n = p.len();
    let dxte = p.x[0] - p.x[n - 1];
    let dyte = p.y[0] - p.y[n - 1];
    let dxs = 0.5 * (-p.xp[0] + p.xp[n - 1]);
    let dys = 0.5 * (-p.yp[0] + p.yp[n - 1]);
    let ante = dxs * dyte - dys * dxte;
    let aste = dxs * dxte + dys * dyte;
    let dste = (dxte * dxte + dyte * dyte).sqrt();
    TrailingEdge {
        ante,
        aste,
        dste,
        sharp: dste < 0.0001 * chord,
    }
}

/// Unit normals at the nodes, averaged at corners. Port of XFOIL `NCALC`.
///
/// NCALC re-splines X and Y with SEGSPL; the result equals the node derivatives
/// `xp`, `yp` already held by the curve, so those are used directly.
fn ncalc(p: &Curve) -> (Vec<f64>, Vec<f64>) {
    let n = p.len();
    let (mut xn, mut yn) = (vec![0.0; n], vec![0.0; n]);
    for i in 0..n {
        let sx = p.yp[i];
        let sy = -p.xp[i];
        let smod = (sx * sx + sy * sy).sqrt();
        if smod == 0.0 {
            xn[i] = -1.0;
            yn[i] = 0.0;
        } else {
            xn[i] = sx / smod;
            yn[i] = sy / smod;
        }
    }
    for i in 0..n - 1 {
        if p.s[i] == p.s[i + 1] {
            let sx = 0.5 * (xn[i] + xn[i + 1]);
            let sy = 0.5 * (yn[i] + yn[i + 1]);
            let smod = (sx * sx + sy * sy).sqrt();
            let (a, b) = if smod == 0.0 {
                (-1.0, 0.0)
            } else {
                (sx / smod, sy / smod)
            };
            xn[i] = a;
            yn[i] = b;
            xn[i + 1] = a;
            yn[i + 1] = b;
        }
    }
    (xn, yn)
}

/// Panel angles. Port of XFOIL `APCALC`.
fn apcalc(p: &Curve, nx: &[f64], ny: &[f64], sharp: bool) -> Vec<f64> {
    let n = p.len();
    let mut ap = vec![0.0; n];
    for i in 0..n - 1 {
        let sx = p.x[i + 1] - p.x[i];
        let sy = p.y[i + 1] - p.y[i];
        ap[i] = if sx == 0.0 && sy == 0.0 {
            (-ny[i]).atan2(-nx[i])
        } else {
            sx.atan2(-sy)
        };
    }
    // TE panel, from node N to node 1
    ap[n - 1] = if sharp {
        PI
    } else {
        let sx = p.x[0] - p.x[n - 1];
        let sy = p.y[0] - p.y[n - 1];
        (-sx).atan2(sy) + PI
    };
    ap
}
