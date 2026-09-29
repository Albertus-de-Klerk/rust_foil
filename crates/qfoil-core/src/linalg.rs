//! Dense linear algebra. Port of XFOIL `xsolve.f` (`GAUSS`, `LUDCMP`, `BAKSUB`).
//!
//! These are hand-ported rather than delegated to a library so that the floating-point
//! operation order matches the reference exactly.

use std::ops::{Index, IndexMut};

/// Linear-algebra errors (Fortran `STOP`s).
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum LinalgError {
    /// LUDCMP's scaling array holds at most `NVX = 500` rows.
    #[error("LUDCMP: {0} unknowns exceed the QFoil limit of 500 (at most 499 panel nodes)")]
    TooLarge(usize),
}

/// Maximum LUDCMP system size (`NVX`).
pub const LUDCMP_MAX: usize = 500;

/// A dense matrix stored column-major, like Fortran arrays.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    /// A zero matrix.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Column-major data.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Column `j` as a slice.
    pub fn col(&self, j: usize) -> &[f64] {
        &self.data[j * self.rows..(j + 1) * self.rows]
    }

    /// Column `j` as a mutable slice.
    pub fn col_mut(&mut self, j: usize) -> &mut [f64] {
        &mut self.data[j * self.rows..(j + 1) * self.rows]
    }

    /// The leading `rows × cols` block.
    pub fn block(&self, rows: usize, cols: usize) -> Vec<f64> {
        (0..cols)
            .flat_map(|j| self.col(j)[..rows].to_vec())
            .collect()
    }
}

impl Index<(usize, usize)> for Matrix {
    type Output = f64;
    #[inline(always)]
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        debug_assert!(i < self.rows && j < self.cols);
        &self.data[j * self.rows + i]
    }
}

impl IndexMut<(usize, usize)> for Matrix {
    #[inline(always)]
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut f64 {
        debug_assert!(i < self.rows && j < self.cols);
        &mut self.data[j * self.rows + i]
    }
}

/// LU factors with row pivots, from [`ludcmp`].
#[derive(Debug, Clone, PartialEq)]
pub struct LuFactors {
    /// L (unit diagonal, below) and U (on and above the diagonal) in one matrix.
    pub lu: Matrix,
    /// Row interchanged with row `j` at step `j` (0-based; Fortran `INDX(J)-1`).
    pub pivots: Vec<usize>,
}

/// Crout LU factorisation with implicit-scaling partial pivoting. Port of XFOIL `LUDCMP`.
///
/// A singular matrix is not detected (division by zero), as in the Fortran.
pub fn ludcmp(mut a: Matrix) -> Result<LuFactors, LinalgError> {
    let n = a.rows();
    if n > LUDCMP_MAX {
        return Err(LinalgError::TooLarge(n));
    }
    let mut vv = vec![0.0; n];
    for i in 0..n {
        let mut aamax = 0.0_f64;
        for j in 0..n {
            aamax = a[(i, j)].abs().max(aamax);
        }
        vv[i] = 1.0 / aamax;
    }
    let mut pivots = vec![0; n];
    for j in 0..n {
        for i in 0..j {
            let mut sum = a[(i, j)];
            for k in 0..i {
                sum -= a[(i, k)] * a[(k, j)];
            }
            a[(i, j)] = sum;
        }
        let mut aamax = 0.0;
        let mut imax = j;
        for i in j..n {
            let mut sum = a[(i, j)];
            for k in 0..j {
                sum -= a[(i, k)] * a[(k, j)];
            }
            a[(i, j)] = sum;
            let dum = vv[i] * sum.abs();
            if dum >= aamax {
                imax = i;
                aamax = dum;
            }
        }
        if j != imax {
            for k in 0..n {
                let dum = a[(imax, k)];
                a[(imax, k)] = a[(j, k)];
                a[(j, k)] = dum;
            }
            vv[imax] = vv[j];
        }
        pivots[j] = imax;
        if j != n - 1 {
            let dum = 1.0 / a[(j, j)];
            for i in j + 1..n {
                a[(i, j)] *= dum;
            }
        }
    }
    Ok(LuFactors { lu: a, pivots })
}

impl LuFactors {
    /// Solves `A x = b` in place. Port of XFOIL `BAKSUB`.
    pub fn baksub(&self, b: &mut [f64]) {
        let a = &self.lu;
        let n = a.rows();
        // Fortran II: first row with a non-zero right-hand side (0 = none yet).
        let mut ii: Option<usize> = None;
        for i in 0..n {
            let ll = self.pivots[i];
            let mut sum = b[ll];
            b[ll] = b[i];
            if let Some(ii) = ii {
                for j in ii..i {
                    sum -= a[(i, j)] * b[j];
                }
            } else if sum != 0.0 {
                ii = Some(i);
            }
            b[i] = sum;
        }
        for i in (0..n).rev() {
            let mut sum = b[i];
            for j in i + 1..n {
                sum -= a[(i, j)] * b[j];
            }
            b[i] = sum / a[(i, i)];
        }
    }
}

/// Gaussian elimination with partial pivoting on the leading `n × n` block of `z`,
/// one right-hand side. Port of XFOIL `GAUSS` (NRHS = 1). `z` is destroyed and `r`
/// receives the solution.
pub fn gauss(n: usize, z: &mut impl IndexMut<(usize, usize), Output = f64>, r: &mut [f64]) {
    for np in 0..n - 1 {
        let np1 = np + 1;
        // max pivot; a later row wins only if strictly larger
        let mut nx = np;
        for k in np1..n {
            if z[(k, np)].abs() - z[(nx, np)].abs() > 0.0 {
                nx = k;
            }
        }
        let pivot = 1.0 / z[(nx, np)];
        z[(nx, np)] = z[(np, np)];
        for l in np1..n {
            let temp = z[(nx, l)] * pivot;
            z[(nx, l)] = z[(np, l)];
            z[(np, l)] = temp;
        }
        let temp = r[nx] * pivot;
        r[nx] = r[np];
        r[np] = temp;
        for k in np1..n {
            let ztmp = z[(k, np)];
            for l in np1..n {
                z[(k, l)] -= ztmp * z[(np, l)];
            }
            r[k] -= ztmp * r[np];
        }
    }
    r[n - 1] /= z[(n - 1, n - 1)];
    for np in (0..n - 1).rev() {
        for k in np + 1..n {
            r[np] -= z[(np, k)] * r[k];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn sample() -> (Matrix, Vec<f64>) {
        let mut a = Matrix::zeros(3, 3);
        let v = [[2.0, 1.0, 1.0], [4.0, -6.0, 0.0], [-2.0, 7.0, 2.0]];
        for i in 0..3 {
            for j in 0..3 {
                a[(i, j)] = v[i][j];
            }
        }
        (a, vec![5.0, -2.0, 9.0]) // solution (1, 1, 2)
    }

    #[test]
    fn lu_solves() {
        let (a, mut b) = sample();
        let lu = ludcmp(a).unwrap();
        lu.baksub(&mut b);
        for (x, e) in b.iter().zip([1.0, 1.0, 2.0]) {
            assert_relative_eq!(*x, e, epsilon = 1e-14);
        }
    }

    #[test]
    fn gauss_solves() {
        let (mut a, mut b) = sample();
        gauss(3, &mut a, &mut b);
        for (x, e) in b.iter().zip([1.0, 1.0, 2.0]) {
            assert_relative_eq!(*x, e, epsilon = 1e-14);
        }
    }

    #[test]
    fn ludcmp_limit() {
        assert_eq!(
            ludcmp(Matrix::zeros(501, 501)).unwrap_err(),
            LinalgError::TooLarge(501)
        );
    }
}
