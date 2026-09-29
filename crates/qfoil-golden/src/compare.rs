//! Element-wise comparison of ported results against golden values.

use std::fmt;

/// Distance between two finite `f64` values in units in the last place.
///
/// `0` means bit-identical (with `+0.0 == -0.0`). Values of different sign count
/// the ULPs through zero. Returns `u64::MAX` if either value is NaN.
pub fn ulp_distance(a: f64, b: f64) -> u64 {
    if a.is_nan() || b.is_nan() {
        return u64::MAX;
    }
    if a == b {
        return 0;
    }
    // Map the IEEE bit pattern onto a monotonic integer line.
    let key = |x: f64| {
        let bits = x.to_bits() as i64;
        if bits < 0 { i64::MIN - bits } else { bits }
    };
    key(a).abs_diff(key(b))
}

/// The worst element of a slice comparison.
#[derive(Debug, Clone, PartialEq)]
pub struct Mismatch {
    /// Record name, for messages.
    pub name: String,
    /// 0-based index of the worst element (Fortran index is `index + 1`).
    pub index: usize,
    /// Ported value.
    pub got: f64,
    /// Golden value.
    pub want: f64,
    /// ULP distance of the worst element.
    pub ulps: u64,
    /// Absolute difference of the worst element.
    pub abs: f64,
    /// Number of elements exceeding the tolerance.
    pub count: usize,
}

impl fmt::Display for Mismatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: {} element(s) out of tolerance; worst at [{}] (Fortran {}): got {:e}, want {:e} ({} ulp, abs {:e})",
            self.name,
            self.count,
            self.index,
            self.index + 1,
            self.got,
            self.want,
            self.ulps,
            self.abs
        )
    }
}

/// Compares `got` with `want` element-wise.
///
/// An element passes if it is within `max_ulps` ULPs **or** within `abs_tol` absolutely.
/// Returns `Ok(max_ulps_seen)` or the worst offending element.
pub fn compare_slices(
    name: &str,
    got: &[f64],
    want: &[f64],
    max_ulps: u64,
    abs_tol: f64,
) -> Result<u64, Mismatch> {
    assert_eq!(
        got.len(),
        want.len(),
        "{name}: length {} vs golden {}",
        got.len(),
        want.len()
    );
    let mut worst: Option<Mismatch> = None;
    let mut max_seen = 0;
    let mut count = 0;
    for (i, (&g, &w)) in got.iter().zip(want).enumerate() {
        let ulps = ulp_distance(g, w);
        let abs = (g - w).abs();
        max_seen = max_seen.max(ulps);
        if ulps > max_ulps && (abs.is_nan() || abs > abs_tol) {
            count += 1;
            if worst.as_ref().is_none_or(|m| ulps > m.ulps) {
                worst = Some(Mismatch {
                    name: name.to_owned(),
                    index: i,
                    got: g,
                    want: w,
                    ulps,
                    abs,
                    count: 0,
                });
            }
        }
    }
    match worst {
        Some(mut m) => {
            m.count = count;
            Err(m)
        }
        None => Ok(max_seen),
    }
}

/// Asserts that ported slices or scalars match golden values.
///
/// ```ignore
/// assert_golden!(ulps = 0; "X" => &pan.x, dump.reals("X"));
/// assert_golden!(ulps = 4, abs = 1e-15;
///     "Y" => &pan.y, dump.reals("Y");
///     "SLE" => pan.sle, dump.real("SLE"));
/// ```
///
/// Every pair is checked before failing, so one run reports all mismatches.
#[macro_export]
macro_rules! assert_golden {
    (ulps = $ulps:expr $(, abs = $abs:expr)?; $($name:literal => $got:expr, $want:expr);+ $(;)?) => {{
        let abs_tol: f64 = 0.0 $(+ $abs)?;
        let mut failures: Vec<String> = Vec::new();
        $(
            let (got_ref, want_ref) = (&$got, &$want);
            let got = $crate::AsSlice::as_f64_slice(got_ref);
            let want = $crate::AsSlice::as_f64_slice(want_ref);
            match $crate::compare_slices($name, got.as_ref(), want.as_ref(), $ulps, abs_tol) {
                Ok(_) => {}
                Err(m) => failures.push(m.to_string()),
            }
        )+
        if !failures.is_empty() {
            panic!("golden mismatch:\n  {}", failures.join("\n  "));
        }
    }};
}

/// Scalars and slices accepted by [`assert_golden!`].
pub trait AsSlice {
    /// View as a slice of `f64`.
    fn as_f64_slice(&self) -> std::borrow::Cow<'_, [f64]>;
}

impl AsSlice for f64 {
    fn as_f64_slice(&self) -> std::borrow::Cow<'_, [f64]> {
        std::borrow::Cow::Owned(vec![*self])
    }
}
impl AsSlice for [f64] {
    fn as_f64_slice(&self) -> std::borrow::Cow<'_, [f64]> {
        std::borrow::Cow::Borrowed(self)
    }
}
impl AsSlice for Vec<f64> {
    fn as_f64_slice(&self) -> std::borrow::Cow<'_, [f64]> {
        std::borrow::Cow::Borrowed(self)
    }
}
impl<T: AsSlice + ?Sized> AsSlice for &T {
    fn as_f64_slice(&self) -> std::borrow::Cow<'_, [f64]> {
        (**self).as_f64_slice()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ulps() {
        assert_eq!(ulp_distance(1.0, 1.0), 0);
        assert_eq!(ulp_distance(0.0, -0.0), 0);
        assert_eq!(ulp_distance(1.0, f64::from_bits(1.0f64.to_bits() + 3)), 3);
        assert_eq!(
            ulp_distance(-f64::MIN_POSITIVE, f64::MIN_POSITIVE),
            2 * f64::MIN_POSITIVE.to_bits()
        );
        assert_eq!(ulp_distance(f64::NAN, 1.0), u64::MAX);
    }

    #[test]
    fn slice_report() {
        let e = compare_slices("v", &[1.0, 2.0, 3.5], &[1.0, 2.0, 3.0], 0, 0.0).unwrap_err();
        assert_eq!((e.index, e.count), (2, 1));
        assert!(compare_slices("v", &[1.0], &[1.0 + 1e-15], 0, 1e-14).is_ok());
    }
}
