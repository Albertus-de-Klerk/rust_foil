//! Airfoil polar analysis: a Rust port of QFOIL 0.9, an XFOIL 6.99 derivative.
//!
//! The numerical code performs no I/O. Each item names the Fortran routine it ports
//! ("Port of XFOIL/QFoil `ROUTINE`"); `docs/PORTING_PLAN.md` maps the whole code base.
//!
//! Licensed GPL-2.0-or-later, like XFOIL and QFOIL. See `NOTICE` for provenance.

pub mod error;
pub mod geometry;
pub mod limits;
pub mod paneling;
pub mod spline;

mod fortran;

pub use error::{GeometryError, ParseError, SplineError};
pub use geometry::{Airfoil, AirfoilSource, Geometry, NacaDesignation};
pub use paneling::{Paneling, PanelingMode, PanelingSettings};
