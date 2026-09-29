//! Cubic splines parameterised by arc length. Port of XFOIL `spline.f`.
//!
//! All routines work on slices. Index conventions: Fortran element `I` is Rust `[I-1]`.
//! The floating-point operation order follows the Fortran source so results are
//! bit-identical to the reference build (checked against the golden dumps).

use crate::error::SplineError;

/// End condition of a spline (`XS1`/`XS2` arguments of SPLIND).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EndCondition {
    /// Zero second derivative (`999.0` in the Fortran).
    ZeroSecondDerivative,
    /// Zero third derivative (`-999.0` in the Fortran).
    ZeroThirdDerivative,
    /// Prescribed first derivative.
    Slope(f64),
}

/// Solves a tridiagonal system in place.
///
/// Port of XFOIL `TRISOL`. `a` is the diagonal, `b` the sub-diagonal (`b[k]` multiplies
/// unknown `k-1` in row `k`), `c` the super-diagonal. The right-hand side `d` is replaced
/// by the solution. `a` and `c` are destroyed.
pub fn trisol(a: &mut [f64], b: &[f64], c: &mut [f64], d: &mut [f64]) {
    let kk = d.len();
    for k in 1..kk {
        let km = k - 1;
        c[km] /= a[km];
        d[km] /= a[km];
        a[k] -= b[k] * c[km];
        d[k] -= b[k] * d[km];
    }
    d[kk - 1] /= a[kk - 1];
    for k in (0..kk - 1).rev() {
        d[k] -= c[k] * d[k + 1];
    }
}

/// Interior rows shared by SPLINE and SPLIND: returns (A, B, C) with XS filled for rows 2..N-1.
fn interior_rows(x: &[f64], xs: &mut [f64], s: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = x.len();
    let (mut a, mut b, mut c) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    for i in 1..n.saturating_sub(1) {
        let dsm = s[i] - s[i - 1];
        let dsp = s[i + 1] - s[i];
        b[i] = dsp;
        a[i] = 2.0 * (dsm + dsp);
        c[i] = dsm;
        xs[i] = 3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm);
    }
    (a, b, c)
}

/// Spline derivatives `xs = dx/ds` with zero second derivative ends. Port of XFOIL `SPLINE`.
pub fn spline(x: &[f64], xs: &mut [f64], s: &[f64]) {
    let n = x.len();
    let (mut a, mut b, mut c) = interior_rows(x, xs, s);
    a[0] = 2.0;
    c[0] = 1.0;
    xs[0] = 3.0 * (x[1] - x[0]) / (s[1] - s[0]);
    b[n - 1] = 1.0;
    a[n - 1] = 2.0;
    xs[n - 1] = 3.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);
    trisol(&mut a, &b, &mut c, &mut xs[..n]);
}

/// Spline derivatives with the given end conditions. Port of XFOIL `SPLIND`.
pub fn splind(x: &[f64], xs: &mut [f64], s: &[f64], start: EndCondition, end: EndCondition) {
    let n = x.len();
    let (mut a, mut b, mut c) = interior_rows(x, xs, s);

    match start {
        EndCondition::ZeroSecondDerivative => {
            a[0] = 2.0;
            c[0] = 1.0;
            xs[0] = 3.0 * (x[1] - x[0]) / (s[1] - s[0]);
        }
        EndCondition::ZeroThirdDerivative => {
            a[0] = 1.0;
            c[0] = 1.0;
            xs[0] = 2.0 * (x[1] - x[0]) / (s[1] - s[0]);
        }
        EndCondition::Slope(v) => {
            a[0] = 1.0;
            c[0] = 0.0;
            xs[0] = v;
        }
    }
    let last = n - 1;
    match end {
        EndCondition::ZeroSecondDerivative => {
            b[last] = 1.0;
            a[last] = 2.0;
            xs[last] = 3.0 * (x[last] - x[last - 1]) / (s[last] - s[last - 1]);
        }
        EndCondition::ZeroThirdDerivative => {
            b[last] = 1.0;
            a[last] = 1.0;
            xs[last] = 2.0 * (x[last] - x[last - 1]) / (s[last] - s[last - 1]);
        }
        EndCondition::Slope(v) => {
            a[last] = 1.0;
            b[last] = 0.0;
            xs[last] = v;
        }
    }
    // A two-point segment with zero-third-derivative ends would be singular.
    if n == 2
        && start == EndCondition::ZeroThirdDerivative
        && end == EndCondition::ZeroThirdDerivative
    {
        b[last] = 1.0;
        a[last] = 2.0;
        xs[last] = 3.0 * (x[last] - x[last - 1]) / (s[last] - s[last - 1]);
    }
    trisol(&mut a, &b, &mut c, &mut xs[..n]);
}

/// Segmented spline: derivative discontinuities where successive `s` values are equal.
///
/// Port of XFOIL `SEGSPL` (zero-third-derivative ends on each segment).
pub fn segspl(x: &[f64], xs: &mut [f64], s: &[f64]) -> Result<(), SplineError> {
    segspld(
        x,
        xs,
        s,
        EndCondition::ZeroThirdDerivative,
        EndCondition::ZeroThirdDerivative,
    )
}

/// Segmented spline with end conditions applied to every segment. Port of XFOIL `SEGSPLD`.
pub fn segspld(
    x: &[f64],
    xs: &mut [f64],
    s: &[f64],
    start: EndCondition,
    end: EndCondition,
) -> Result<(), SplineError> {
    let n = x.len();
    if n < 2 {
        return Err(SplineError::TooFewPoints(n));
    }
    if s[0] == s[1] {
        return Err(SplineError::DuplicatedEndPoint { first: true });
    }
    if s[n - 1] == s[n - 2] {
        return Err(SplineError::DuplicatedEndPoint { first: false });
    }
    // Fortran: DO ISEG=2, N-2 ; IF(S(ISEG).EQ.S(ISEG+1)) -> 0-based iseg = 1..=n-3
    let mut iseg0 = 0;
    for iseg in 1..n.saturating_sub(2) {
        if s[iseg] == s[iseg + 1] {
            let r = iseg0..=iseg;
            splind(&x[r.clone()], &mut xs[r.clone()], &s[r], start, end);
            iseg0 = iseg + 1;
        }
    }
    let r = iseg0..n;
    splind(&x[r.clone()], &mut xs[r.clone()], &s[r], start, end);
    Ok(())
}

/// Bisection for the interval containing `ss`. Returns the 0-based index of its upper node.
///
/// Port of the search loop shared by SEVAL, DEVAL, D2VAL, CURV and CURVS.
#[inline]
fn upper_node(ss: f64, s: &[f64]) -> usize {
    let mut ilow = 0;
    let mut i = s.len() - 1;
    while i - ilow > 1 {
        // 0-based midpoint equals Fortran (I+ILOW)/2 - 1 for the shifted indices.
        let imid = (i + ilow) / 2;
        if ss < s[imid] {
            i = imid;
        } else {
            ilow = imid;
        }
    }
    i
}

/// Cubic coefficients of the interval ending at node `i`: `(ds, t, cx1, cx2)`.
#[inline]
fn local(ss: f64, x: &[f64], xs: &[f64], s: &[f64], i: usize) -> (f64, f64, f64, f64) {
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;
    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    (ds, t, cx1, cx2)
}

/// Spline value `x(ss)`. Port of XFOIL `SEVAL`.
pub fn seval(ss: f64, x: &[f64], xs: &[f64], s: &[f64]) -> f64 {
    let i = upper_node(ss, s);
    let (_, t, cx1, cx2) = local(ss, x, xs, s, i);
    t * x[i] + (1.0 - t) * x[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2)
}

/// Spline derivative `dx/ds(ss)`. Port of XFOIL `DEVAL`.
pub fn deval(ss: f64, x: &[f64], xs: &[f64], s: &[f64]) -> f64 {
    let i = upper_node(ss, s);
    let (ds, t, cx1, cx2) = local(ss, x, xs, s, i);
    let d = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    d / ds
}

/// Spline second derivative `d²x/ds²(ss)`. Port of XFOIL `D2VAL`.
pub fn d2val(ss: f64, x: &[f64], xs: &[f64], s: &[f64]) -> f64 {
    let i = upper_node(ss, s);
    let (ds, t, cx1, cx2) = local(ss, x, xs, s, i);
    let d = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
    d / (ds * ds)
}

/// First and second parametric derivatives of both coordinates on one interval.
struct Derivs {
    xd: f64,
    xdd: f64,
    yd: f64,
    ydd: f64,
    cx: (f64, f64),
    cy: (f64, f64),
    ds: f64,
}

fn derivs(ss: f64, x: &[f64], xs: &[f64], y: &[f64], ys: &[f64], s: &[f64]) -> Derivs {
    let i = upper_node(ss, s);
    let (ds, t, cx1, cx2) = local(ss, x, xs, s, i);
    let xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    let xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
    let cy1 = ds * ys[i - 1] - y[i] + y[i - 1];
    let cy2 = ds * ys[i] - y[i] + y[i - 1];
    let yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
    let ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;
    Derivs {
        xd,
        xdd,
        yd,
        ydd,
        cx: (cx1, cx2),
        cy: (cy1, cy2),
        ds,
    }
}

/// Curvature of the splined curve at `ss`. Port of XFOIL `CURV`.
pub fn curv(ss: f64, x: &[f64], xs: &[f64], y: &[f64], ys: &[f64], s: &[f64]) -> f64 {
    let d = derivs(ss, x, xs, y, ys, s);
    let sd = (d.xd * d.xd + d.yd * d.yd).sqrt().max(0.001 * d.ds);
    (d.xd * d.ydd - d.yd * d.xdd) / (sd * sd * sd)
}

/// Derivative of curvature along the curve at `ss`. Port of XFOIL `CURVS`.
pub fn curvs(ss: f64, x: &[f64], xs: &[f64], y: &[f64], ys: &[f64], s: &[f64]) -> f64 {
    let d = derivs(ss, x, xs, y, ys, s);
    let xddd = 6.0 * d.cx.0 + 6.0 * d.cx.1;
    let yddd = 6.0 * d.cy.0 + 6.0 * d.cy.1;
    let sd = (d.xd * d.xd + d.yd * d.yd).sqrt().max(0.001 * d.ds);
    let bot = sd * sd * sd;
    let dbotdt = 3.0 * sd * (d.xd * d.xdd + d.yd * d.ydd);
    let top = d.xd * d.ydd - d.yd * d.xdd;
    let dtopdt = d.xd * yddd - d.yd * xddd;
    (dtopdt * bot - dbotdt * top) / (bot * bot)
}

/// Result of [`sinvrt`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Inversion {
    /// `s` such that `x(s) = xi`, or the initial guess if the iteration failed.
    pub s: f64,
    /// Whether the Newton iteration converged.
    pub converged: bool,
}

/// Inverts the spline: finds `s` with `x(s) = xi` from the initial guess `si`.
///
/// Port of XFOIL `SINVRT`. On failure the initial guess is returned, as in the Fortran.
pub fn sinvrt(si: f64, xi: f64, x: &[f64], xs: &[f64], s: &[f64]) -> Inversion {
    let n = s.len();
    let mut sv = si;
    for _ in 0..10 {
        let res = seval(sv, x, xs, s) - xi;
        let resp = deval(sv, x, xs, s);
        let ds = -res / resp;
        sv += ds;
        if (ds / (s[n - 1] - s[0])).abs() < 1.0e-5 {
            return Inversion {
                s: sv,
                converged: true,
            };
        }
    }
    Inversion {
        s: si,
        converged: false,
    }
}

/// Cumulative arc length of a polyline. Port of XFOIL `SCALC`.
pub fn scalc(x: &[f64], y: &[f64], s: &mut [f64]) {
    s[0] = 0.0;
    for i in 1..x.len() {
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }
}

/// A 2-D curve splined in arc length: SCALC followed by SEGSPL of x and y.
#[derive(Debug, Clone, PartialEq)]
pub struct Curve {
    /// x coordinates.
    pub x: Vec<f64>,
    /// y coordinates.
    pub y: Vec<f64>,
    /// Arc length at each node.
    pub s: Vec<f64>,
    /// dx/ds at each node.
    pub xp: Vec<f64>,
    /// dy/ds at each node.
    pub yp: Vec<f64>,
}

impl Curve {
    /// Computes arc length and segmented splines for the given nodes.
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Result<Self, SplineError> {
        let n = x.len();
        let mut s = vec![0.0; n];
        scalc(&x, &y, &mut s);
        Self::with_arc_length(x, y, s)
    }

    /// Splines nodes with a given arc-length parameter (no SCALC).
    pub fn with_arc_length(x: Vec<f64>, y: Vec<f64>, s: Vec<f64>) -> Result<Self, SplineError> {
        let n = x.len();
        let mut xp = vec![0.0; n];
        let mut yp = vec![0.0; n];
        segspl(&x, &mut xp, &s)?;
        segspl(&y, &mut yp, &s)?;
        Ok(Self { x, y, s, xp, yp })
    }

    /// Number of nodes.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Whether the curve has no nodes.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// `x(s)`.
    pub fn x_at(&self, s: f64) -> f64 {
        seval(s, &self.x, &self.xp, &self.s)
    }

    /// `y(s)`.
    pub fn y_at(&self, s: f64) -> f64 {
        seval(s, &self.y, &self.yp, &self.s)
    }

    /// Signed curvature at `s`.
    pub fn curvature(&self, s: f64) -> f64 {
        curv(s, &self.x, &self.xp, &self.y, &self.yp, &self.s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn trisol_solves_system() {
        // [2 1 0; 1 2 1; 0 1 2] u = [3 4 3] -> u = [1 1 1]
        let mut a = vec![2.0, 2.0, 2.0];
        let b = vec![0.0, 1.0, 1.0];
        let mut c = vec![1.0, 1.0, 0.0];
        let mut d = vec![3.0, 4.0, 3.0];
        trisol(&mut a, &b, &mut c, &mut d);
        for v in d {
            assert_relative_eq!(v, 1.0, epsilon = 1e-15);
        }
    }

    #[test]
    fn spline_reproduces_cubic_interpolant_on_line() {
        let s: Vec<f64> = (0..6).map(|i| i as f64 * 0.3).collect();
        let x: Vec<f64> = s.iter().map(|v| 2.0 * v + 1.0).collect();
        let mut xs = vec![0.0; 6];
        segspl(&x, &mut xs, &s).unwrap();
        for &d in &xs {
            assert_relative_eq!(d, 2.0, epsilon = 1e-13);
        }
        assert_relative_eq!(seval(0.75, &x, &xs, &s), 2.5, epsilon = 1e-13);
        assert_relative_eq!(deval(0.75, &x, &xs, &s), 2.0, epsilon = 1e-13);
        assert_relative_eq!(d2val(0.75, &x, &xs, &s), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn circle_curvature() {
        let n = 181;
        let (x, y): (Vec<f64>, Vec<f64>) = (0..n)
            .map(|i| {
                let t = std::f64::consts::PI * i as f64 / (n - 1) as f64;
                (t.cos(), t.sin())
            })
            .unzip();
        let c = Curve::new(x, y).unwrap();
        let mid = c.s[n / 2];
        // chord-length parameterisation: ~2.5e-5 curvature error at 1 degree spacing
        assert_relative_eq!(c.curvature(mid), 1.0, epsilon = 1e-4);
        let inv = sinvrt(mid + 0.01, c.x_at(mid), &c.x, &c.xp, &c.s);
        assert!(inv.converged);
        assert_relative_eq!(inv.s, mid, epsilon = 1e-6);
    }

    #[test]
    fn segments_split_at_duplicate_arc_length() {
        // Corner at node 2 (0-based) duplicated: two independent straight segments.
        let x = vec![0.0, 1.0, 2.0, 2.0, 2.0, 2.0];
        let y = vec![0.0, 0.0, 0.0, 0.0, 1.0, 2.0];
        let c = Curve::new(x, y).unwrap();
        assert_eq!(c.s[2], c.s[3]);
        assert_relative_eq!(c.xp[1], 1.0, epsilon = 1e-14);
        assert_relative_eq!(c.yp[4], 1.0, epsilon = 1e-14);
    }

    #[test]
    fn duplicated_end_point_is_an_error() {
        let s = vec![0.0, 0.0, 1.0];
        let mut xs = vec![0.0; 3];
        assert!(segspl(&[0.0, 0.0, 1.0], &mut xs, &s).is_err());
    }
}
