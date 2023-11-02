use axum::extract::{Json, State};
use axum::routing::post;
use axum::Router;
use sage_core::database::IndexedDatabase;
use sage_core::mass::Tolerance;
use sage_core::scoring::{Feature, Scorer};
use sage_core::spectrum::{Precursor, RawSpectrum, SpectrumProcessor};
use serde::{Deserialize, Serialize};
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

    precursor_mz: f32,
    precursor_charge: u8,
    mz: Vec<f32>,
    intensity: Vec<f32>,
}

async fn score_v1(
    State(db): State<Arc<IndexedDatabase>>,
    Json(query): Json<ScoreRequest>,
) -> Result<Json<Vec<Feature>>, (axum::http::StatusCode, String)> {
    let scorer = Scorer {
        db: &db,
        precursor_tol: query.precursor_tolerance,
        fragment_tol: query.fragment_tolerance,
        min_matched_peaks: 4,
        min_isotope_err: 0,
        max_isotope_err: 0,
        min_precursor_charge: 1,
        max_precursor_charge: 6,
        max_fragment_charge: Some(1),
        min_fragment_mass: 125.0,
        max_fragment_mass: 2500.0,
        chimera: false,
        report_psms: query.report_psms,
        wide_window: false,
    };

    let spectra = RawSpectrum {
        file_id: 0,
        ms_level: 2,
        id: "real-time".into(),
        precursors: vec![Precursor {
            mz: query.precursor_mz,
            intensity: None,
            charge: Some(query.precursor_charge),
            spectrum_ref: None,
            isolation_window: None,
        }],
        representation: sage_core::spectrum::Representation::Centroid,
        scan_start_time: 0.0,
        ion_injection_time: 0.0,
        total_ion_current: 0.0,
        mz: query.mz,
        intensity: query.intensity,
    };

    let spectra =
        SpectrumProcessor::new(150, 0.0, 2000.0, query.deisotope).process(spectra.clone());

    let scores = scorer.score(&spectra);

    Ok(Json(scores))
}

pub type Error = Box<dyn std::error::Error + Send + Sync + 'static>;

#[tokio::main]
async fn main() -> Result<(), Error> {
    tracing_subscriber::fmt()
        .with_ansi(true)
        .with_max_level(tracing::Level::TRACE)
        .init();

    let parameters: sage_core::database::Builder =
        serde_json::from_str(&tokio::fs::read_to_string("params.json").await.unwrap()).unwrap();

    let parameters = parameters.make_parameters();
    let contents = tokio::fs::read_to_string(&parameters.fasta).await.unwrap();

    let fasta =
        sage_core::fasta::Fasta::parse(contents, &parameters.decoy_tag, parameters.generate_decoys);

    let db = parameters.build(fasta);

    let app = Router::new()
        .route("/v1/score/", post(score_v1))
        .with_state(Arc::new(db))
        .layer(CorsLayer::very_permissive())
        .layer(CompressionLayer::new().gzip(true).deflate(true));

    let addr = std::net::SocketAddr::from(([127, 0, 0, 1], 3000));

    axum::Server::bind(&addr)
        .serve(app.into_make_service())
        .await
        .unwrap();
    Ok(())
}
