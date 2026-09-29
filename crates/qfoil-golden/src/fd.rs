//! Finite-difference verification of analytic derivatives.

/// Central difference of `g` with respect to argument `i` at `x0`.
///
/// Step `h = 1e-6 * max(|x_i|, 1e-2)`; the result is O(h²) accurate for smooth `g`.
pub fn central_difference<const N: usize>(
    g: impl Fn([f64; N]) -> f64,
    x0: [f64; N],
    i: usize,
) -> f64 {
    let h = 1e-6 * x0[i].abs().max(1e-2);
    let (mut xp, mut xm) = (x0, x0);
    xp[i] += h;
    xm[i] -= h;
    (g(xp) - g(xm)) / (xp[i] - xm[i])
}

/// Panics unless `analytic` matches the finite difference `fd` to `rel` relative accuracy
/// (with an absolute floor of `rel * scale` for derivatives near zero).
pub fn check_derivative(
    what: &str,
    arg: usize,
    x0: &[f64],
    analytic: f64,
    fd: f64,
    rel: f64,
    scale: f64,
) {
    let err = (analytic - fd).abs();
    let tol = rel * analytic.abs().max(fd.abs()).max(scale);
    assert!(
        err <= tol,
        "{what}: d/dx[{arg}] at {x0:?}: analytic {analytic:e} vs finite difference {fd:e} (err {err:e} > {tol:e})"
    );
}

/// Checks analytic partial derivatives of a closure against central differences.
///
/// ```ignore
/// assert_derivatives!(
///     |x: [f64; 3]| hst(x[0], x[1], x[2]), at [2.5, 1000.0, 0.1],
///     value = |r: &Correlation| r.v, rel = 1e-6,
///     { 0 => |r: &Correlation| r.hk, 1 => |r: &Correlation| r.rt, 2 => |r: &Correlation| r.msq }
/// );
/// ```
/// `scale` (optional, default 1e-6) is the absolute floor for derivatives near zero,
/// relative to which `rel` is applied.
#[macro_export]
macro_rules! assert_derivatives {
    ($f:expr, at $x0:expr, value = $val:expr, rel = $rel:expr $(, scale = $scale:expr)?,
     { $($i:literal => $der:expr),+ $(,)? }) => {{
        let f = $f;
        let val = $val;
        let x0 = $x0;
        let scale: f64 = 1e-6 $(* 0.0 + $scale)?;
        let r0 = f(x0);
        $(
            let analytic: f64 = ($der)(&r0);
            let fd = $crate::fd::central_difference(|x| val(&f(x)), x0, $i);
            $crate::fd::check_derivative(stringify!($der), $i, &x0, analytic, fd, $rel, scale);
        )+
    }};
}
