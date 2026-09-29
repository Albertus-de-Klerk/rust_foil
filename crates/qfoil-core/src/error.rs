//! Error types. Non-convergence of a viscous point is **not** an error: it is reported
//! as a status on the operating point. Errors cover invalid input and failures that
//! would stop QFoil (`STOP` statements, rejected commands).

/// Spline construction errors (XFOIL `SEGSPL` stops).
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SplineError {
    /// Fewer than two points.
    #[error("a spline needs at least 2 points, got {0}")]
    TooFewPoints(usize),
    /// First or last point duplicated (zero-length end segment).
    #[error("{} input point duplicated", if *first { "first" } else { "last" })]
    DuplicatedEndPoint {
        /// `true` for the first point, `false` for the last.
        first: bool,
    },
}

/// Errors reading airfoil coordinates (XFOIL `AREAD`).
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ParseError {
    /// Fewer than two data lines.
    #[error("unexpected end of coordinate data")]
    UnexpectedEnd,
    /// The second line of a labelled file does not start with two numbers.
    #[error("unrecognisable file format (line {line})")]
    UnrecognisedFormat {
        /// 1-based line number.
        line: usize,
    },
    /// A coordinate line could not be read.
    #[error("bad coordinate line {line}: `{text}`")]
    BadLine {
        /// 1-based line number.
        line: usize,
        /// Line content.
        text: String,
    },
    /// Multi-element MSES files are not supported (QFoil asks interactively).
    #[error("multi-element coordinate files are not supported")]
    MultiElement,
}

/// Errors building geometry or panels.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum GeometryError {
    /// Unsupported NACA designation.
    #[error("NACA designation {0} not implemented (4-digit, or 5-digit 210xx..250xx)")]
    NacaDesignation(u32),
    /// Too few points to define an airfoil.
    #[error("airfoil needs at least 4 points, got {0}")]
    TooFewPoints(usize),
    /// More points than QFoil's arrays allow (`IQX` limits, see PORTING_PLAN D13).
    #[error("{what}: {got} points exceed the QFoil limit of {max}")]
    TooManyPoints {
        /// Which limit.
        what: &'static str,
        /// Requested count.
        got: usize,
        /// Maximum allowed.
        max: usize,
    },
    /// The leading-edge point lies at the first buffer interval, where PANGEN would read
    /// outside its arrays (Fortran `SB(I-2)` with `I = 2`).
    #[error("leading edge too close to the first buffer point")]
    LeadingEdgeAtEnd,
    /// Spline failure.
    #[error(transparent)]
    Spline(#[from] SplineError),
}
