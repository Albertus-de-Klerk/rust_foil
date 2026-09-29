//! `qfoil`: polar sweeps from the command line, writing QFoil's polar file format.

mod polar_file;

use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use qfoil_core::settings::{BlParams, FlowConditions};
use qfoil_core::{
    Airfoil, AlphaSchedule, NacaDesignation, PanelingMode, PanelingSettings, PolarSettings,
    Settings, ViscousSettings, analyse_polar,
};

/// Airfoil polar analysis (Rust port of QFOIL 0.9, an XFOIL 6.99 derivative).
///
/// Every angle of attack is solved from a cold start, as QBlade runs QFoil. Output is
/// QFoil's polar file format (converged points only).
#[derive(Debug, Parser)]
#[command(version, about)]
struct Cli {
    /// Airfoil coordinate file (plain, labelled or MSES single-element format).
    #[arg(required_unless_present = "naca", conflicts_with = "naca")]
    airfoil: Option<PathBuf>,

    /// Built-in NACA 4- or 5-digit section instead of a file (e.g. 0012, 23012).
    #[arg(long)]
    naca: Option<String>,

    /// Reynolds number. Omit together with --inviscid for an inviscid polar.
    #[arg(long, required_unless_present = "inviscid")]
    re: Option<f64>,

    /// Inviscid analysis.
    #[arg(long)]
    inviscid: bool,

    /// Mach number.
    #[arg(long, default_value_t = 0.0)]
    mach: f64,

    /// Angles of attack: START END STEP (degrees).
    #[arg(long, num_args = 3, value_names = ["START", "END", "STEP"], allow_negative_numbers = true,
          default_values_t = [-10.0, 20.0, 0.5])]
    alpha: Vec<f64>,

    /// Critical amplification exponent Ncrit (both sides).
    #[arg(long, default_value_t = 9.0)]
    ncrit: f64,

    /// Forced transition x/c on the top side (1 = free).
    #[arg(long, default_value_t = 1.0)]
    xtr_top: f64,

    /// Forced transition x/c on the bottom side (1 = free).
    #[arg(long, default_value_t = 1.0)]
    xtr_bot: f64,

    /// Newton iteration limit per point.
    #[arg(long, default_value_t = 100)]
    iter: usize,

    /// Re-panel with this many nodes (PANE). Default: files use their own points (as
    /// QFoil's LOAD does), NACA sections get 160 curvature-based panels.
    #[arg(long)]
    panels: Option<usize>,

    /// QFoil wake-drag correction factor GWAKE.
    #[arg(long, default_value_t = 0.40)]
    gwake: f64,

    /// Output file (default: standard output).
    #[arg(short, long)]
    output: Option<PathBuf>,
}

fn run(cli: Cli) -> Result<String, String> {
    let airfoil = match (&cli.naca, &cli.airfoil) {
        (Some(code), _) => {
            let des: NacaDesignation = code.parse().map_err(|e| format!("--naca {code}: {e}"))?;
            Airfoil::naca(des).map_err(|e| e.to_string())?
        }
        (None, Some(path)) => {
            let text =
                std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
            Airfoil::from_dat(&text).map_err(|e| format!("{}: {e}", path.display()))?
        }
        (None, None) => unreachable!("clap requires one of them"),
    };
    let paneling = match cli.panels {
        Some(npan) => PanelingMode::Generate(PanelingSettings {
            npan,
            ..PanelingSettings::default()
        }),
        None => PanelingMode::Auto,
    };
    let settings = Settings {
        flow: FlowConditions {
            mach: cli.mach,
            reynolds: cli.re.unwrap_or(0.0),
            ..FlowConditions::default()
        },
        paneling,
        ..Settings::default()
    };
    let viscous = (!cli.inviscid).then(|| ViscousSettings {
        ncrit: [cli.ncrit; 2],
        xtrip: [cli.xtr_top, cli.xtr_bot],
        max_iterations: cli.iter,
        bl: BlParams {
            gwake: cli.gwake,
            ..BlParams::default()
        },
        ..ViscousSettings::default()
    });
    let alpha = AlphaSchedule::Range {
        start: cli.alpha[0],
        end: cli.alpha[1],
        step: cli.alpha[2],
    };
    let polar = analyse_polar(
        &airfoil,
        &PolarSettings {
            alpha,
            settings,
            viscous,
        },
    )
    .map_err(|e| e.to_string())?;
    Ok(polar_file::write(&polar))
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let output = cli.output.clone();
    match run(cli) {
        Ok(text) => match output {
            Some(path) => match std::fs::write(&path, text) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("qfoil: {}: {e}", path.display());
                    ExitCode::FAILURE
                }
            },
            None => {
                print!("{text}");
                ExitCode::SUCCESS
            }
        },
        Err(e) => {
            eprintln!("qfoil: {e}");
            ExitCode::FAILURE
        }
    }
}
