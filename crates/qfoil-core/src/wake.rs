//! Wake geometry. Port of XFOIL `XYWAKE`, `SETEXP` (xutils.f) and `QWCALC`.

use crate::inviscid::influence::{FieldPoint, Influence, NodeRef, Strengths, psilin};
use crate::limits::IWX;
use crate::paneling::Paneling;

/// Wake nodes behind the TE (Fortran indices `N+1 .. N+NW`; wake-local index `k` is
/// Fortran `N+1+k`).
#[derive(Debug, Clone, PartialEq)]
pub struct Wake {
    /// x coordinates.
    pub x: Vec<f64>,
    /// y coordinates.
    pub y: Vec<f64>,
    /// Arc length, continuing from the airfoil's `S(N)`.
    pub s: Vec<f64>,
    /// Unit normal x components.
    pub nx: Vec<f64>,
    /// Unit normal y components.
    pub ny: Vec<f64>,
    /// Wake panel normal angles (`APANEL(N+1+k)`); the last entry is unused.
    pub apanel: Vec<f64>,
    /// Angle of attack the trajectory was traced for (`AWAKE`).
    pub alfa: f64,
}

impl Wake {
    /// Number of wake nodes `NW`.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Whether the wake has no nodes.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Number of wake nodes for `n` airfoil nodes: `N/12 + 10*INT(WAKLEN)`, at most `IWX`.
    pub fn node_count(n: usize, wake_length: f64) -> usize {
        (n / 12 + 10 * wake_length.trunc() as usize).min(IWX)
    }

    /// Traces the wake as a streamline of the current inviscid flow. Port of XFOIL `XYWAKE`.
    pub fn trace(pan: &Paneling, st: &Strengths<'_>, wake_length: f64) -> Self {
        let (x, y, s, xp, yp) = (
            &pan.nodes.x,
            &pan.nodes.y,
            &pan.nodes.s,
            &pan.nodes.xp,
            &pan.nodes.yp,
        );
        let n = x.len();
        let nw = Self::node_count(n, wake_length);

        let ds1 = 0.5 * (s[1] - s[0] + s[n - 1] - s[n - 2]);
        let snew = setexp(ds1, wake_length * pan.chord, nw);

        let xte = 0.5 * (x[0] + x[n - 1]);
        let yte = 0.5 * (y[0] + y[n - 1]);

        let mut w = Self {
            x: vec![0.0; nw],
            y: vec![0.0; nw],
            s: vec![0.0; nw],
            nx: vec![0.0; nw],
            ny: vec![0.0; nw],
            apanel: vec![0.0; nw],
            alfa: st.alfa,
        };
        let mut inf = Influence::new(n, 0);

        // first wake point a tiny distance behind the TE, along the TE bisector
        let sx = 0.5 * (yp[n - 1] - yp[0]);
        let sy = 0.5 * (xp[0] - xp[n - 1]);
        let smod = (sx * sx + sy * sy).sqrt();
        w.nx[0] = sx / smod;
        w.ny[0] = sy / smod;
        w.x[0] = xte - 0.0001 * w.ny[0];
        w.y[0] = yte + 0.0001 * w.nx[0];
        w.s[0] = s[n - 1];

        // Each point sets the normal of the next from the local stream-function gradient.
        let gradient = |k: usize, xk: f64, yk: f64, inf: &mut Influence| -> (f64, f64) {
            let fp = |nx, ny| FieldPoint {
                node: NodeRef::Wake(k),
                x: xk,
                y: yk,
                nx,
                ny,
            };
            psilin(pan, st, fp(1.0, 0.0), false, inf);
            let psi_x = inf.psi_ni;
            psilin(pan, st, fp(0.0, 1.0), false, inf);
            (psi_x, inf.psi_ni)
        };

        let (psi_x, psi_y) = gradient(0, w.x[0], w.y[0], &mut inf);
        let g = (psi_x * psi_x + psi_y * psi_y).sqrt();
        w.nx[1] = -psi_x / g;
        w.ny[1] = -psi_y / g;
        w.apanel[0] = psi_y.atan2(psi_x);

        for k in 1..nw {
            let ds = snew[k] - snew[k - 1];
            w.x[k] = w.x[k - 1] - ds * w.ny[k];
            w.y[k] = w.y[k - 1] + ds * w.nx[k];
            w.s[k] = w.s[k - 1] + ds;
            if k == nw - 1 {
                break;
            }
            let (psi_x, psi_y) = gradient(k, w.x[k], w.y[k], &mut inf);
            let g = (psi_x * psi_x + psi_y * psi_y).sqrt();
            w.nx[k + 1] = -psi_x / g;
            w.ny[k + 1] = -psi_y / g;
            w.apanel[k] = psi_y.atan2(psi_x);
        }
        w
    }

    /// Inviscid tangential velocity on the wake for α = 0° and 90°. Port of XFOIL `QWCALC`.
    ///
    /// Returns `NW` values for each case; entry 0 (the point at the TE) copies the airfoil's
    /// last node, `qinvu_te`.
    pub fn qwcalc(&self, pan: &Paneling, st: &Strengths<'_>, qinvu_te: [f64; 2]) -> [Vec<f64>; 2] {
        let nw = self.len();
        let mut q = [vec![0.0; nw], vec![0.0; nw]];
        q[0][0] = qinvu_te[0];
        q[1][0] = qinvu_te[1];
        let mut inf = Influence::new(pan.len(), 0);
        for k in 1..nw {
            let fp = FieldPoint {
                node: NodeRef::Wake(k),
                x: self.x[k],
                y: self.y[k],
                nx: self.nx[k],
                ny: self.ny[k],
            };
            psilin(pan, st, fp, false, &mut inf);
            q[0][k] = inf.qtan1;
            q[1][k] = inf.qtan2;
        }
        q
    }
}

/// Geometrically stretched spacing `s[0] = 0`, `s[1] = ds1`, `s[nn-1] = smax`.
/// Port of XFOIL `SETEXP`.
///
/// `nn < 3` stops QFoil (`SETEXP: Cannot fill array`); it cannot occur for a valid wake
/// (`NW >= 10`), so it is an invariant here.
pub fn setexp(ds1: f64, smax: f64, nn: usize) -> Vec<f64> {
    assert!(nn >= 3, "SETEXP: cannot fill array, n too small");
    let sigma = smax / ds1;
    let nex = nn - 1;
    let rnex = nex as f64;
    let rni = 1.0 / rnex;

    // quadratic for the initial geometric-ratio guess
    let aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
    let bbb = rnex * (rnex - 1.0) / 2.0;
    let ccc = rnex - sigma;
    let disc = (bbb * bbb - 4.0 * aaa * ccc).max(0.0);
    let mut ratio = if nex == 2 {
        -ccc / bbb + 1.0
    } else {
        (-bbb + disc.sqrt()) / (2.0 * aaa) + 1.0
    };

    if ratio != 1.0 {
        // Newton iteration for the ratio. RATIO**NEX has a variable integer exponent, which
        // gfortran evaluates with libgcc __powidf2; f64::powi uses the same algorithm.
        let nexi = nex as i32;
        for _ in 0..100 {
            let sigman = (ratio.powi(nexi) - 1.0) / (ratio - 1.0);
            let res = sigman.powf(rni) - sigma.powf(rni);
            let dresdr = rni * sigman.powf(rni) * (rnex * ratio.powi(nexi - 1) - sigman)
                / (ratio.powi(nexi) - 1.0);
            let dratio = -res / dresdr;
            ratio += dratio;
            if dratio.abs() < 1.0e-5 {
                break;
            }
        }
    }
    let mut s = vec![0.0; nn];
    let mut ds = ds1;
    for i in 1..nn {
        s[i] = s[i - 1] + ds;
        ds *= ratio;
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn setexp_hits_end_point() {
        let s = setexp(0.01, 1.0, 23);
        assert_eq!(s[0], 0.0);
        assert_relative_eq!(s[1], 0.01);
        assert_relative_eq!(s[22], 1.0, epsilon = 1e-4);
        let r = (s[2] - s[1]) / (s[1] - s[0]);
        assert_relative_eq!((s[3] - s[2]) / (s[2] - s[1]), r, epsilon = 1e-12);
    }
}
