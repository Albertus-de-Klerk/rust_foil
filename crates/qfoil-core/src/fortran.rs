//! Helpers that reproduce gfortran's arithmetic exactly.

/// `x**n` for a small integer literal `n`, expanded the way gfortran (GCC `powi`) does at `-O2`.
///
/// GCC expands constant integer powers with multiplication chains; `f64::powi` goes through
/// a runtime routine whose rounding can differ. Only the exponents used in the Fortran
/// source are provided.
#[inline(always)]
pub(crate) fn powi(x: f64, n: u32) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        2 => x * x,
        3 => x * x * x,
        4 => {
            let x2 = x * x;
            x2 * x2
        }
        _ => unreachable!("powi({n}) not used in the Fortran source"),
    }
}

/// π as XFOIL computes it: `PI = 4.0*ATAN(1.0)` (INIT).
pub(crate) const PI: f64 = std::f64::consts::PI;

#[cfg(test)]
mod tests {
    #[test]
    fn pi_matches_fortran_definition() {
        assert_eq!(4.0 * 1.0f64.atan(), super::PI);
    }
}
