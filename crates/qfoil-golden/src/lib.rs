//! Readers for the QFoil golden reference data in `tests/golden/`.
//!
//! The dump format is written by the instrumentation patch
//! `tools/patches/dump/0001-golden-dumps.patch` and is documented in `tests/golden/README.md`.
//! Every dump file holds named records. A directory of dumps is ordered by the call
//! sequence number in the file name (`NNNNN_<tag>.txt`).
//!
//! Records and dumps are kept in [`BTreeMap`]s, so iteration and error output are
//! deterministic.

use std::collections::BTreeMap;
use std::fmt;
use std::path::{Path, PathBuf};

mod compare;
mod polar;

pub use compare::{AsSlice, Mismatch, compare_slices, ulp_distance};
pub use polar::{GoldenPoint, GoldenPolar};

/// Errors raised while reading golden files.
#[derive(Debug, thiserror::Error)]
pub enum GoldenError {
    /// The file or directory could not be read.
    #[error("cannot read {path}: {source}")]
    Io {
        /// Offending path.
        path: PathBuf,
        /// Underlying I/O error.
        source: std::io::Error,
    },
    /// The file content does not follow the dump format.
    #[error("{path}:{line}: {message}")]
    Format {
        /// Offending file.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
        /// Description.
        message: String,
    },
    /// A requested record or dump is absent.
    #[error("missing {0}")]
    Missing(String),
}

/// One named record of a dump file.
#[derive(Debug, Clone, PartialEq)]
pub enum Record {
    /// `#S`: real scalar.
    Real(f64),
    /// `#I`: integer scalar.
    Int(i64),
    /// `#R`: real vector (Fortran index 1 is element 0).
    Reals(Vec<f64>),
    /// `#V`: integer vector.
    Ints(Vec<i64>),
    /// `#M`: real matrix stored column-major, as Fortran does.
    Matrix {
        /// Number of rows.
        rows: usize,
        /// Number of columns.
        cols: usize,
        /// Column-major data, `rows * cols` values.
        data: Vec<f64>,
    },
}

/// One dump file: a tag and its records by name.
#[derive(Debug, Clone)]
pub struct Dump {
    /// Call sequence number from the file name.
    pub seq: u32,
    /// Dump tag (`pangen`, `ggcalc`, ...).
    pub tag: String,
    /// Records by name.
    pub records: BTreeMap<String, Record>,
    /// Source file, for messages.
    pub path: PathBuf,
}

impl Dump {
    /// Parses one dump file.
    pub fn read(path: impl AsRef<Path>) -> Result<Self, GoldenError> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|source| GoldenError::Io {
            path: path.to_owned(),
            source,
        })?;
        let seq = path
            .file_name()
            .and_then(|n| n.to_str())
            .and_then(|n| n.split('_').next())
            .and_then(|n| n.parse().ok())
            .unwrap_or(0);
        Self::parse(&text, seq, path)
    }

    fn parse(text: &str, seq: u32, path: &Path) -> Result<Self, GoldenError> {
        let err = |line: usize, message: String| GoldenError::Format {
            path: path.to_owned(),
            line,
            message,
        };
        let mut lines = text.lines().enumerate().map(|(i, l)| (i + 1, l.trim()));
        let mut tag = String::new();
        let mut records = BTreeMap::new();

        while let Some((ln, line)) = lines.next() {
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            let (kind, name) = match fields.as_slice() {
                [kind, name, ..] => (*kind, (*name).to_owned()),
                _ => return Err(err(ln, format!("malformed header `{line}`"))),
            };
            let field = |k: usize| {
                fields
                    .get(k)
                    .copied()
                    .ok_or_else(|| err(ln, format!("missing field {k} in `{line}`")))
            };
            let count = |k: usize| -> Result<usize, GoldenError> {
                field(k)?
                    .parse()
                    .map_err(|_| err(ln, format!("bad count in `{line}`")))
            };
            let real =
                |v: &str, l: usize| parse_f64(v).ok_or_else(|| err(l, format!("bad real `{v}`")));
            let int = |v: &str, l: usize| {
                v.parse::<i64>()
                    .map_err(|_| err(l, format!("bad integer `{v}`")))
            };

            let record = match kind {
                "#TAG" => {
                    tag = name;
                    continue;
                }
                "#S" => Record::Real(real(field(2)?, ln)?),
                "#I" => Record::Int(int(field(2)?, ln)?),
                "#R" | "#V" | "#M" => {
                    let rows = count(2)?;
                    let cols = if kind == "#M" { count(3)? } else { 1 };
                    let body: Vec<(usize, &str)> = lines.by_ref().take(rows * cols).collect();
                    if body.len() != rows * cols {
                        return Err(err(
                            ln,
                            format!("`{name}`: expected {} values", rows * cols),
                        ));
                    }
                    match kind {
                        "#V" => Record::Ints(
                            body.iter()
                                .map(|&(l, v)| int(v, l))
                                .collect::<Result<_, _>>()?,
                        ),
                        "#R" => Record::Reals(
                            body.iter()
                                .map(|&(l, v)| real(v, l))
                                .collect::<Result<_, _>>()?,
                        ),
                        _ => Record::Matrix {
                            rows,
                            cols,
                            data: body
                                .iter()
                                .map(|&(l, v)| real(v, l))
                                .collect::<Result<_, _>>()?,
                        },
                    }
                }
                other => return Err(err(ln, format!("unknown record kind `{other}`"))),
            };
            records.insert(name, record);
        }
        Ok(Self {
            seq,
            tag,
            records,
            path: path.to_owned(),
        })
    }

    fn get(&self, name: &str) -> Result<&Record, GoldenError> {
        self.records.get(name).ok_or_else(|| {
            GoldenError::Missing(format!("record `{name}` in {}", self.path.display()))
        })
    }

    /// Real scalar record.
    pub fn real(&self, name: &str) -> f64 {
        match self.get(name) {
            Ok(Record::Real(v)) => *v,
            other => panic!("{name}: expected real scalar, got {other:?}"),
        }
    }

    /// Integer scalar record.
    pub fn int(&self, name: &str) -> i64 {
        match self.get(name) {
            Ok(Record::Int(v)) => *v,
            other => panic!("{name}: expected integer scalar, got {other:?}"),
        }
    }

    /// Real vector record.
    pub fn reals(&self, name: &str) -> &[f64] {
        match self.get(name) {
            Ok(Record::Reals(v)) => v,
            other => panic!("{name}: expected real vector, got {other:?}"),
        }
    }

    /// Integer vector record.
    pub fn ints(&self, name: &str) -> &[i64] {
        match self.get(name) {
            Ok(Record::Ints(v)) => v,
            other => panic!("{name}: expected integer vector, got {other:?}"),
        }
    }

    /// Matrix record as `(rows, cols, column-major data)`.
    pub fn matrix(&self, name: &str) -> (usize, usize, &[f64]) {
        match self.get(name) {
            Ok(Record::Matrix { rows, cols, data }) => (*rows, *cols, data),
            other => panic!("{name}: expected matrix, got {other:?}"),
        }
    }
}

/// Fortran `ES26.17E3` output, which may use a 3-digit exponent without `E` for |e|>99.
fn parse_f64(s: &str) -> Option<f64> {
    s.parse().ok().or_else(|| {
        // e.g. "1.0-100" (Fortran drops the E for 3-digit exponents in some formats)
        let pos = s.rfind(['+', '-']).filter(|&p| p > 0)?;
        format!("{}E{}", &s[..pos], &s[pos..]).parse().ok()
    })
}

/// All dumps of one reference run, ordered by call sequence.
#[derive(Debug, Clone)]
pub struct DumpSet {
    /// Directory the dumps came from.
    pub dir: PathBuf,
    /// Dumps by call sequence number.
    pub dumps: BTreeMap<u32, Dump>,
}

impl DumpSet {
    /// Reads every `NNNNN_<tag>.txt` in `dir`.
    pub fn read(dir: impl AsRef<Path>) -> Result<Self, GoldenError> {
        let dir = dir.as_ref().to_owned();
        let entries = std::fs::read_dir(&dir).map_err(|source| GoldenError::Io {
            path: dir.clone(),
            source,
        })?;
        let mut dumps = BTreeMap::new();
        for entry in entries {
            let path = entry
                .map_err(|source| GoldenError::Io {
                    path: dir.clone(),
                    source,
                })?
                .path();
            let is_dump = path.file_name().and_then(|n| n.to_str()).is_some_and(|n| {
                n.ends_with(".txt") && n.as_bytes().first().is_some_and(u8::is_ascii_digit)
            });
            if is_dump {
                let d = Dump::read(&path)?;
                dumps.insert(d.seq, d);
            }
        }
        Ok(Self { dir, dumps })
    }

    /// All dumps with `tag`, in call order.
    pub fn all<'a>(&'a self, tag: &'a str) -> impl Iterator<Item = &'a Dump> + 'a {
        self.dumps.values().filter(move |d| d.tag == tag)
    }

    /// The first dump with `tag`.
    pub fn first(&self, tag: &str) -> &Dump {
        self.dumps
            .values()
            .find(|d| d.tag == tag)
            .unwrap_or_else(|| panic!("no `{tag}` dump in {}", self.dir.display()))
    }

    /// The last dump with `tag`.
    pub fn last(&self, tag: &str) -> &Dump {
        self.dumps
            .values()
            .rfind(|d| d.tag == tag)
            .unwrap_or_else(|| panic!("no `{tag}` dump in {}", self.dir.display()))
    }

    /// The command script that produced this run (`cmd.txt`), if present.
    pub fn command_script(&self) -> Option<String> {
        std::fs::read_to_string(self.dir.join("cmd.txt")).ok()
    }
}

/// Where the golden data lives: `$QFOIL_GOLDEN` or `<workspace>/tests/golden`.
pub fn golden_root() -> PathBuf {
    std::env::var_os("QFOIL_GOLDEN").map_or_else(
        || Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden"),
        PathBuf::from,
    )
}

/// A committed fixture run from `tests/golden/fixtures/<name>`.
pub fn fixture(name: &str) -> DumpSet {
    let dir = golden_root().join("fixtures").join(name);
    DumpSet::read(&dir).unwrap_or_else(|e| panic!("fixture {name}: {e}"))
}

/// A full (gitignored) dump run from `tests/golden/dumps/<name>`.
///
/// Panics with instructions when the dumps have not been generated.
pub fn full_dump(name: &str) -> DumpSet {
    let dir = golden_root().join("dumps").join(name);
    DumpSet::read(&dir).unwrap_or_else(|e| {
        panic!(
            "{e}\nfull dumps are not committed: run tools/build_reference.sh && tools/gen_golden.sh"
        )
    })
}

/// Path of a golden airfoil coordinate file.
pub fn airfoil_path(name: &str) -> PathBuf {
    golden_root().join("airfoils").join(format!("{name}.dat"))
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Record::Real(v) => write!(f, "{v:e}"),
            Record::Int(v) => write!(f, "{v}"),
            Record::Reals(v) => write!(f, "[{} reals]", v.len()),
            Record::Ints(v) => write!(f, "[{} ints]", v.len()),
            Record::Matrix { rows, cols, .. } => write!(f, "[{rows}x{cols}]"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_all_record_kinds() {
        let text = "#TAG demo\n#I N 2\n#S A 1.5E+000\n#R V 2\n 1.0E+000\n-2.5E-003\n#V IV 2\n 3\n 4\n#M M 2 2\n1\n2\n3\n4\n";
        let d = Dump::parse(text, 7, Path::new("x")).unwrap();
        assert_eq!(d.tag, "demo");
        assert_eq!(d.int("N"), 2);
        assert_eq!(d.real("A"), 1.5);
        assert_eq!(d.reals("V"), &[1.0, -2.5e-3]);
        assert_eq!(d.ints("IV"), &[3, 4]);
        assert_eq!(d.matrix("M"), (2, 2, &[1.0, 2.0, 3.0, 4.0][..]));
    }

    #[test]
    fn fortran_three_digit_exponent() {
        assert_eq!(parse_f64("1.25-100"), Some(1.25e-100));
        assert_eq!(parse_f64("-3.00000000000000000E-003"), Some(-3e-3));
    }
}
