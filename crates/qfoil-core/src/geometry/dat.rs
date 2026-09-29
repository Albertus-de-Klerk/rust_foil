//! Coordinate-file parsing. Port of XFOIL `AREAD` and the number scanning of `GETFLT`.
//!
//! Supported layouts (AREAD `ITYPE`):
//! 1. plain: every line is `x y`;
//! 2. labelled: a name line followed by `x y` lines;
//! 3. MSES single element: name line, a line of four or more domain numbers, then `x y`.
//!
//! Lines starting with `#` or `!` are comments. A line with fewer than two numbers is
//! skipped. Multi-element files (`999.0 999.0` separators) are rejected; QFoil prompts
//! for an element number interactively there.

use super::{Airfoil, AirfoilSource};
use crate::error::ParseError;

/// AREAD reads lines into `CHARACTER*80`, and GETFLT scans at most 128 characters.
const LINE_LEN: usize = 80;

/// Scans up to `max` numbers from a line like GETFLT. Returns `None` if a token that
/// GETFLT would try to read is not a number (list-directed READ error).
fn getflt(line: &str, max: usize) -> Option<Vec<f64>> {
    let line = &line[..line.len().min(LINE_LEN)];
    let line = line.split('!').next().unwrap_or_default();
    line.split([' ', '\t', ','])
        .filter(|t| !t.is_empty())
        .take(max)
        .map(parse_fortran_real)
        .collect()
}

/// A Fortran real literal: also accepts `D`/`d` exponents and a trailing `.`.
fn parse_fortran_real(token: &str) -> Option<f64> {
    token
        .parse()
        .ok()
        .or_else(|| token.replace(['D', 'd'], "E").parse().ok())
}

fn is_comment(line: &str) -> bool {
    line.starts_with(['#', '!'])
}

pub(super) fn parse(text: &str) -> Result<Airfoil, ParseError> {
    let lines: Vec<(usize, &str)> = text
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l.trim_end_matches('\r')))
        .filter(|(_, l)| !is_comment(l))
        .collect();
    let (&(_, line1), &(ln2, line2)) = match lines.as_slice() {
        [a, b, ..] => (a, b),
        _ => return Err(ParseError::UnexpectedEnd),
    };

    let (name, body) = match getflt(line1, 2) {
        // Two valid numbers on the first line: plain file.
        Some(v) if v.len() >= 2 => (String::new(), &lines[..]),
        // Otherwise a name line; the second line decides labelled vs MSES.
        _ => match getflt(line2, 4) {
            Some(v) if v.len() >= 4 => (line1.trim().to_owned(), &lines[2..]),
            Some(v) if v.len() >= 2 => (line1.trim().to_owned(), &lines[1..]),
            _ => return Err(ParseError::UnrecognisedFormat { line: ln2 }),
        },
    };

    let (mut x, mut y) = (Vec::new(), Vec::new());
    for &(ln, line) in body {
        let v = getflt(line, 2).ok_or_else(|| ParseError::BadLine {
            line: ln,
            text: line.to_owned(),
        })?;
        if v.len() < 2 {
            continue;
        }
        if v[0] == 999.0 && v[1] == 999.0 {
            return Err(ParseError::MultiElement);
        }
        x.push(v[0]);
        y.push(v[1]);
    }
    Ok(Airfoil {
        name,
        x,
        y,
        source: AirfoilSource::File,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labelled_file() {
        let a = parse("E387\n 1.0 0.0\n0.5, 0.05 ! comment\n# skip\n0.0 0.0\n").unwrap();
        assert_eq!(a.name, "E387");
        assert_eq!(a.x, vec![1.0, 0.5, 0.0]);
        assert_eq!(a.y, vec![0.0, 0.05, 0.0]);
    }

    #[test]
    fn plain_file_and_fortran_exponent() {
        let a = parse("1.0 0.0\n5.0D-1\t1.0d-2\n0 0\n").unwrap();
        assert_eq!(a.name, "");
        assert_eq!(a.x, vec![1.0, 0.5, 0.0]);
        assert_eq!(a.y[1], 0.01);
    }

    #[test]
    fn mses_header_is_skipped() {
        let a = parse("foil\n -2.0 3.0 -2.5 3.5\n1 0\n0 0\n1 0.001\n").unwrap();
        assert_eq!(a.x.len(), 3);
    }

    #[test]
    fn rejects_bad_formats() {
        assert_eq!(parse("only one line\n"), Err(ParseError::UnexpectedEnd));
        assert!(matches!(
            parse("name\nfoo bar\n1 0\n"),
            Err(ParseError::UnrecognisedFormat { line: 2 })
        ));
        assert_eq!(
            parse("a\n1 0\n999 999\n0 0\n"),
            Err(ParseError::MultiElement)
        );
    }
}
