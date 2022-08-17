use axum::extract::{Json, Path, Query};
use axum::response::{IntoResponse, Response};
use axum::routing::{get, post};
use axum::{Extension, Router};
use http::header::{ACCEPT, ACCEPT_ENCODING, AUTHORIZATION, CONTENT_TYPE, ORIGIN};
use http::StatusCode;
use sage::database::IndexedDatabase;
use sage::ion_series::Kind;
use sage::mass::{Tolerance, H2O, NH3, PROTON};
use sage::mzml::{MzMlReader, Spectrum};
use sage::scoring::{Percolator, Scorer};
use sage::spectrum::{Peak, ProcessedSpectrum, SpectrumProcessor};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::BufReader;
use std::sync::Arc;
use tower_http::{compression::CompressionLayer, cors::CorsLayer};
use tracing::info;

#[derive(Deserialize, Serialize)]
pub struct ScoreRequest {
    precursor_tolerance: Tolerance,
    fragment_tolerance: Tolerance,
    report_psms: usize,
    chimera: bool,
    deisotope: bool,
}

pub struct State {
    db: IndexedDatabase,
    spectra: Vec<Spectrum>,
}

fn find_scan_id(spectra: &[Spectrum], scan: usize) -> Option<&Spectrum> {
    if let Some(spectrum) = spectra.get(scan) {
        if spectrum.scan_id == scan {
            return Some(spectrum);
        }
    }

    let idx = spectra.binary_search_by(|a| a.scan_id.cmp(&scan)).ok()?;
    spectra.get(idx)
}

#[derive(Deserialize)]
pub struct Peptide {
    sequence: String,
    modifications: HashMap<char, f32>,
    nterm: Option<f32>,
    fragment_tol: Tolerance,
    deisotope: Option<bool>,
}

#[derive(Serialize)]
pub struct MatchedPeaks {
    mz: Option<f32>,
    intensity: Option<f32>,
    charge: Option<u8>,
    fragment_mz: f32,
    fragment_loss: f32,
    fragment_kind: Kind,
    fragment_idx: usize,
}

fn score_peptide(
    state: &State,
    mut query: ProcessedSpectrum,
    request: &Peptide,
) -> Result<Vec<MatchedPeaks>, (StatusCode, String)> {
    let digest = sage::fasta::Digest {
        protein: "unknown",
        sequence: &request.sequence,
    };
    let mut peptide = sage::peptide::Peptide::try_from(&digest)
        .map_err(|e| (StatusCode::BAD_REQUEST, format!("Invalid AA: {}", e)))?;

    if let Some(nterm) = request.nterm {
        peptide.set_nterm_mod(nterm);
    }

    for (r, m) in &request.modifications {
        peptide.static_mod(*r, *m);
    }

    let mut matches = Vec::new();
    // let mut mz = query.peaks.iter().map(|peak| peak.mass).collect::<Vec<_>>();
    // mz.sort_unstable_by(|a, b| a.total_cmp(&b));
    query
        .peaks
        .sort_unstable_by(|a, b| a.mass.total_cmp(&b.mass));

    let charge = query.precursors[0].charge.unwrap_or(2);

    for kind in [Kind::B, Kind::Y] {
        for (idx, ion) in sage::ion_series::IonSeries::new(&peptide, kind)
            .enumerate()
            .filter(|(_, ion)| ion.monoisotopic_mass <= 2000.0)
        {
            for charge in 1..charge {
                // for loss in [0.0, H2O, NH3] {
                let loss = 0.0;
                let mass = (ion.monoisotopic_mass - loss) / charge as f32;

                if let Some(peak) =
                    sage::spectrum::select_closest_peak(&query.peaks, mass, request.fragment_tol)
                {
                    matches.push(MatchedPeaks {
                        mz: Some(peak.mass + PROTON),
                        intensity: Some(peak.intensity),
                        charge: Some(charge),
                        fragment_mz: ion.monoisotopic_mass + PROTON,
                        fragment_loss: loss,
                        fragment_kind: ion.kind,
                        fragment_idx: match ion.kind {
                            Kind::Y => peptide.sequence.len() - idx,
                            Kind::B => idx + 1,
                        },
                    });
                }
                // }
            }
        }
        //     } else if kind == Kind::B {
        //         if let Some(peak) = most_intense_peak(
        //             &query.peaks,
        //             ion.monoisotopic_mass - 18.0105,
        //             request.fragment_tol,
        //         ) {
        //             matches.push(MatchedPeaks {
        //                 mz: Some(peak.mass + PROTON - 18.0105),
        //                 intensity: Some(peak.intensity),
        //                 fragment_mz: ion.monoisotopic_mass + PROTON - 18.0105,
        //                 fragment_kind: ion.kind,
        //                 fragment_idx: match ion.kind {
        //                     Kind::Y => peptide.sequence.len() - idx,
        //                     Kind::B => idx + 1,
        //                 },
        //             });
        //         } else {
        //             matches.push(MatchedPeaks {
        //                 mz: None,
        //                 intensity: None,
        //                 fragment_mz: ion.monoisotopic_mass + PROTON,
        //                 fragment_kind: ion.kind,
        //                 fragment_idx: match ion.kind {
        //                     Kind::Y => peptide.sequence.len() - idx,
        //                     Kind::B => idx + 1,
        //                 },
        //             });
        //         }
        //     } else {
        //         matches.push(MatchedPeaks {
        //             mz: None,
        //             intensity: None,
        //             fragment_mz: ion.monoisotopic_mass + PROTON,
        //             fragment_kind: ion.kind,
        //             fragment_idx: match ion.kind {
        //                 Kind::Y => peptide.sequence.len() - idx,
        //                 Kind::B => idx + 1,
        //             },
        //         });
        //     }
        // }
    }
    Ok(matches)
}

async fn score_spectrum(
    Path(scan_id): Path<usize>,
    Json(query): Json<ScoreRequest>,
    Extension(state): Extension<Arc<State>>,
) -> Result<Response<String>, (StatusCode, String)> {
    info!(message = "score_spectrum", scan_id = scan_id);
    let scorer = Scorer::new(
        &state.db,
        query.precursor_tolerance,
        query.fragment_tolerance,
        -1,
        3,
        None,
        query.chimera,
    );
    let spectra = find_scan_id(&state.spectra, scan_id)
        .ok_or_else(|| (StatusCode::BAD_REQUEST, "Cannot find scan id".into()))?;

    let spectra = SpectrumProcessor::new(150, 2000.0, query.deisotope).process(spectra.clone());

    if spectra.level != 2 {
        return Err((StatusCode::BAD_REQUEST, "Not an MS2 scan".into()));
    }

    let scores = scorer.score(&spectra, query.report_psms);

    info!(
        message = "score_spectrum",
        scan_id = scan_id,
        results = scores.len()
    );
    let body = serde_json::to_string_pretty(&scores)
        .map_err(|e| (StatusCode::BAD_REQUEST, e.to_string()))?;

    let builder = Response::builder()
        .header("content-type", "application/json")
        .status(StatusCode::OK);
    Ok(builder.body(body).unwrap())
}

async fn score_spectrum_peptide(
    Path(scan_id): Path<usize>,
    Json(peptide): Json<Peptide>,
    Extension(state): Extension<Arc<State>>,
) -> Result<Json<Vec<MatchedPeaks>>, (StatusCode, String)> {
    info!(
        message = "score_spectrum_peptide",
        scan_id = scan_id,
        sequence = peptide.sequence.as_str()
    );
    let spectra = find_scan_id(&state.spectra, scan_id)
        .ok_or_else(|| (StatusCode::BAD_REQUEST, "Cannot find scan id".into()))?;

    let spectra = SpectrumProcessor::new(100, 2000.0, peptide.deisotope.unwrap_or(false))
        .process(spectra.clone());

    if spectra.level != 2 {
        return Err((StatusCode::BAD_REQUEST, "Not an MS2 scan".into()));
    }

    score_peptide(&state, spectra, &peptide).map(Json)
}

#[derive(Deserialize)]
struct SpectrumQuery {
    deisotope: bool,
    max_peaks: usize,
}

async fn get_spectrum(
    Path(scan_id): Path<usize>,
    Query(deisotope): Query<SpectrumQuery>,
    Extension(state): Extension<Arc<State>>,
) -> Result<Json<ProcessedSpectrum>, (StatusCode, String)> {
    info!(message = "get_spectrum", scan_id = scan_id);
    let spectra = find_scan_id(&state.spectra, scan_id)
        .ok_or_else(|| (StatusCode::BAD_REQUEST, "Cannot find scan id".into()))?;

    let processed = SpectrumProcessor::new(deisotope.max_peaks, 2000.0, deisotope.deisotope)
        .process(spectra.clone());

    Ok(Json(processed))
}

pub type Error = Box<dyn std::error::Error + Send + Sync + 'static>;

#[tokio::main]
async fn main() -> Result<(), Error> {
    tracing_subscriber::fmt()
        .with_ansi(true)
        .with_max_level(tracing::Level::INFO)
        .json()
        .init();

    let db = std::fs::read_to_string("params.json")?;
    let builder: sage::database::Builder = serde_json::from_str(&db)?;
    let db = builder.make_parameters().build().unwrap();
    info!(message = "database created");

    // let file = std::fs::File::open("b1906_293T_proteinID_01A_QE3_122212.mzML")?;
    let path = std::env::args()
        .nth(1)
        .expect("expecting mzML path argument");
    let file = std::fs::File::open(path)?;
    let read = BufReader::new(file);

    let spectra = MzMlReader::default().parse(read)?;

    let state = State { db, spectra };
    let state = Arc::new(state);

    info!(message = "ready to serve!");

    // Set up CORS
    let cors_layer = CorsLayer::new()
        .allow_credentials(true)
        .allow_headers(vec![
            ACCEPT,
            ACCEPT_ENCODING,
            AUTHORIZATION,
            CONTENT_TYPE,
            ORIGIN,
        ])
        .allow_methods(tower_http::cors::Any)
        .allow_origin(tower_http::cors::Any);

    let app = Router::new()
        .route("/spectrum/:scan_id", post(score_spectrum))
        .route("/spectrum/:scan_id/peptide", post(score_spectrum_peptide))
        .route("/spectrum/:scan_id", get(get_spectrum))
        .layer(Extension(state))
        .layer(cors_layer)
        .layer(CompressionLayer::new().gzip(true).deflate(true));

    let addr = std::net::SocketAddr::from(([127, 0, 0, 1], 3000));

    axum::Server::bind(&addr)
        .serve(app.into_make_service())
        .await
        .unwrap();

    println!("Hello, world!");

    Ok(())
}
