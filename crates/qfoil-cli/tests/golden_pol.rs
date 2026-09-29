//! The CLI output must equal QFoil's polar files byte for byte.

use std::process::Command;

use qfoil_golden::{airfoil_path, golden_root};

/// Runs `qfoil` for a golden case and returns (our file, QFoil's file).
fn run(case: &str, re: &str, alpha: [&str; 3]) -> (String, String) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_qfoil"));
    match case.strip_prefix("naca") {
        Some(code) => cmd.args(["--naca", code]),
        None => cmd.arg(airfoil_path(case)),
    };
    cmd.args(["--re", re, "--alpha", alpha[0], alpha[1], alpha[2]]);
    let out = cmd.output().expect("run qfoil");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let ours = String::from_utf8(out.stdout).unwrap();
    let theirs = std::fs::read_to_string(
        golden_root()
            .join("polars")
            .join(format!("{case}_re{re}.pol")),
    )
    .unwrap();
    (ours, theirs)
}

#[test]
fn header_and_rows_match_qfoil() {
    let (ours, theirs) = run("naca0012", "1e6", ["0", "2", "1"]);
    let header: Vec<&str> = theirs.lines().take(12).collect();
    let ours_lines: Vec<&str> = ours.lines().collect();
    assert_eq!(ours_lines[..12], header[..], "header");
    for row in &ours_lines[12..] {
        assert!(
            theirs.lines().any(|l| l == *row),
            "row not in QFoil's file: {row}"
        );
    }
    assert_eq!(ours_lines.len(), 15);
}

#[test]
#[ignore = "all 915 points: cargo test --release -p qfoil-cli -- --ignored"]
fn all_polar_files_match_qfoil() {
    let cases: Vec<(&str, &str)> = ["naca0012", "naca0020", "naca4412", "du91w2250", "e387"]
        .into_iter()
        .flat_map(|c| ["1e5", "1e6", "5e6"].map(move |r| (c, r)))
        .collect();
    std::thread::scope(|s| {
        let handles: Vec<_> = cases
            .iter()
            .map(|&(c, r)| s.spawn(move || (c, r, run(c, r, ["-10", "20", "0.5"]))))
            .collect();
        for h in handles {
            let (c, r, (ours, theirs)) = h.join().unwrap();
            assert_eq!(ours, theirs, "{c} Re {r}");
        }
    });
}
