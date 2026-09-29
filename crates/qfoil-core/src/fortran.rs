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

/// `ATAN2(y, x)` continued from a nearby angle `thold` without jumping across the branch
/// cut. Port of XFOIL `ATANC` (xutils.f).
pub(crate) fn atanc(y: f64, x: f64, thold: f64) -> f64 {
    // Fortran DATA TPI /6.2831853071795864769/ rounds to TAU (checked in the tests).
    const TPI: f64 = std::f64::consts::TAU;
    let thnew = y.atan2(x);
    let dthet = thnew - thold;
    // angle change cannot exceed ±π: remove multiples of 2π (INT truncates toward zero)
    let dtcorr = dthet - TPI * ((dthet + PI.copysign(dthet)) / TPI).trunc();
    thold + dtcorr
}

#[cfg(test)]
mod tests {
    #[test]
    fn pi_matches_fortran_definition() {
        assert_eq!(4.0 * 1.0f64.atan(), super::PI);
    }

    #[test]
    #[allow(clippy::approx_constant, clippy::excessive_precision)]
    fn tpi_literal_is_tau() {
        assert_eq!(6.2831853071795864769_f64, std::f64::consts::TAU);
    }

    #[test]
    fn atanc_continues_across_branch_cut() {
        let pi = super::PI;
        assert!((super::atanc(-1.0, -1.0, 0.75 * pi) - 1.25 * pi).abs() < 1e-15);
    }
}
