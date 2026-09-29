//! QFoil/XFOIL polar save-file format (`POLWRIT`, iopol.f), byte-compatible with QFoil's
//! `PACC` output.

use std::fmt::Write as _;

use qfoil_core::{Polar, PolarPoint};

/// Fortran `Fw.d` edit descriptor: right-justified, `d` decimals, the leading zero dropped
/// if the field is too narrow for it, asterisks if the value still does not fit.
pub fn fortran_f(v: f64, w: usize, d: usize) -> String {
    let mut s = format!("{v:.d$}");
    if s.len() > w {
        if let Some(rest) = s.strip_prefix("0.") {
            s = format!(".{rest}");
        } else if let Some(rest) = s.strip_prefix("-0.") {
            s = format!("-.{rest}");
        }
    }
    if s.len() > w {
        "*".repeat(w)
    } else {
        format!("{s:>w$}")
    }
}

/// Column specs of the data rows: (name, width, decimals), from PINDEX.INC CPOLFORM.
const COLUMNS: [(&str, usize, usize); 5] = [
    ("alpha", 7, 3),
    ("CL", 9, 4),
    ("CD", 10, 5),
    ("CDp", 10, 5),
    ("CM", 9, 4),
];
/// Side columns (Xtr, Itr), each written for top and bottom.
const SIDE_COLUMNS: [(&str, usize, usize); 2] = [("Xtr", 9, 4), ("Itr", 9, 4)];

/// Writes the column-label and dash lines exactly as POLWRIT builds them.
fn column_header(out: &mut String) {
    let mut label = vec![b' '; 128];
    let mut dash = vec![b' '; 128];
    // Fortran strings are 1-based: KL, KD start at 1
    let (mut kl, mut kd) = (1usize, 1usize);
    let put = |buf: &mut Vec<u8>, from: usize, text: &[u8]| {
        buf[from - 1..from - 1 + text.len()].copy_from_slice(text);
    };
    for (name, w, _) in COLUMNS {
        let nblank = (w as i64 - name.len() as i64 + 2).max(0) as usize / 2;
        put(&mut label, kl + 1 + nblank, name.as_bytes());
        kl += w;
        put(&mut dash, kd + 2, &b"-".repeat(w - 1));
        kd += w;
    }
    for (name, w, _) in SIDE_COLUMNS {
        let nblank = ((w as i64 - name.len() as i64 - 2) / 2).max(0) as usize;
        for prefix in ["Top_", "Bot_"] {
            put(
                &mut label,
                kl + 1 + nblank,
                format!("{prefix}{name}").as_bytes(),
            );
            kl += w;
            put(&mut dash, kd + 2, &b"-".repeat(w - 1));
            kd += w;
        }
    }
    out.push_str(std::str::from_utf8(&label[..kl]).unwrap());
    out.push('\n');
    out.push_str(std::str::from_utf8(&dash[..kd]).unwrap());
    out.push('\n');
}

/// The header of a polar file (POLWRIT with LHEAD).
pub fn header(polar: &Polar) -> String {
    let s = &polar.settings;
    let flow = &s.settings.flow;
    let (ncrit, xtrip) = s
        .viscous
        .as_ref()
        .map_or(([0.0; 2], [1.0; 2]), |v| (v.ncrit, v.xtrip));
    let mut out = String::new();
    // list-directed WRITE(LU,*) ' ' produces two blanks
    let blank = "  \n";
    out.push_str(blank);
    // 8000 FORMAT(7X,A,9X,'Version',F5.2): QFoil keeps XFOIL's code name and version
    let _ = writeln!(out, "{:7}XFOIL{:9}Version{}", "", "", fortran_f(6.99, 5, 2));
    out.push_str(blank);
    // 9001: NAME is CHARACTER*48
    let _ = writeln!(out, " Calculated polar for: {:<48}", polar.name);
    out.push_str(blank);
    // 9005 FORMAT(1X,I1,I2,2A29)
    let _ = writeln!(
        out,
        " 1 1{:<29}{:<29}",
        " Reynolds number fixed", "   Mach number fixed"
    );
    out.push_str(blank);
    // 9011
    let _ = writeln!(
        out,
        " xtrf = {} (top)    {} (bottom)  ",
        fortran_f(xtrip[0], 7, 3),
        fortran_f(xtrip[1], 9, 3)
    );
    // 9015
    let _ = writeln!(
        out,
        " Mach = {}     Re = {} e 6     Ncrit = {}{}",
        fortran_f(flow.mach, 7, 3),
        fortran_f(flow.reynolds / 1.0e6, 9, 3),
        fortran_f(ncrit[0], 7, 3),
        fortran_f(ncrit[1], 7, 3)
    );
    out.push_str(blank);
    column_header(&mut out);
    out
}

/// One data row.
pub fn row(p: &PolarPoint) -> String {
    let values = [p.alpha, p.cl, p.cd, p.cdp, p.cm];
    let mut s = String::from(" ");
    for (v, (_, w, d)) in values.iter().zip(COLUMNS) {
        s.push_str(&fortran_f(*v, w, d));
    }
    for (vals, (_, w, d)) in [p.xtr, p.itr].iter().zip(SIDE_COLUMNS) {
        for v in vals {
            s.push_str(&fortran_f(*v, w, d));
        }
    }
    s
}

/// A complete polar file. Unconverged points are omitted, as QFoil does.
pub fn write(polar: &Polar) -> String {
    let mut out = header(polar);
    for p in polar.converged() {
        out.push_str(&row(p));
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::fortran_f;

    #[test]
    fn f_edit_descriptor() {
        assert_eq!(fortran_f(-3.7e-8, 9, 4), "  -0.0000");
        assert_eq!(fortran_f(0.00539, 10, 5), "   0.00539");
        assert_eq!(fortran_f(-0.5, 6, 4), "-.5000");
        assert_eq!(fortran_f(12345.0, 5, 2), "*****");
        assert_eq!(fortran_f(6.99, 5, 2), " 6.99");
    }
}
