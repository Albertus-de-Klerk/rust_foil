//! Airfoil polar analysis: a Rust port of QFOIL 0.9, an XFOIL 6.99 derivative.
//!
//! The numerical code performs no I/O. Each item names the Fortran routine it ports
//! ("Port of XFOIL/QFoil `ROUTINE`"); `docs/PORTING_PLAN.md` maps the whole code base.
//!
//! Licensed GPL-2.0-or-later, like XFOIL and QFOIL. See `NOTICE` for provenance.

pub mod error;
pub mod geometry;
pub mod inviscid;
pub mod limits;
pub mod linalg;
pub mod operating;
pub mod paneling;
pub mod settings;
pub mod spline;
pub mod wake;

mod fortran;

pub use error::{GeometryError, ParseError, SplineError};
pub use geometry::{Airfoil, AirfoilSource, Geometry, NacaDesignation};
pub use inviscid::InviscidSolution;
pub use operating::OperatingPoint;
pub use paneling::{Paneling, PanelingMode, PanelingSettings};
pub use settings::Settings;
pub use wake::Wake;
