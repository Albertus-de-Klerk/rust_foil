//! QFoil array dimensions (`XFOIL.INC`, QFoil values, PORTING_PLAN D13).
//!
//! The Rust port allocates to the actual problem size. These constants only reproduce
//! behaviour that depends on them (e.g. the NACA point count) and the input limits
//! QFoil enforces.

/// Surface panel nodes + 6.
pub const IQX: usize = 1400;
/// Wake panel nodes: `IQX/8 + 2`.
pub const IWX: usize = IQX / 8 + 2;
