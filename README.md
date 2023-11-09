# SageRTS

This crate implements a basic HTTP server that performs real-time database search using the [Sage](https://github.com/lazear/sage) search engine. By leveraging Sage's design as a library, it's possible to import core functionality into other Rust-based projects, like this one - implementing millisecond-scale real-time open-search in ~100 lines of code!

## Building from scratch

**It is very important to build with the `--release` flag!**

- Install the Rust programming language compiler toolchain
- Run the following, which will clone this repository, download the human proteome from UniProt, and start a HTTP server for real-time search
- Configuration is currently performed via a single `params.json` file located in the same directory as the binary

```sh
git clone https://github.com/lazear/sage-rts.git
cd sage-rts

curl https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz | gzip -d > UP000005640_9606.fasta

cat << EOF > params.json
{
  "bucket_size": 4096,
  "fragment_min_mz": 120.0,
  "fragment_max_mz": 2000.0,
  "peptide_min_len": 6,
  "peptide_max_len": 50,
  "peptide_min_mass": 300.0,
  "peptide_max_mass": 5000.0,
  "missed_cleavages": 2,
  "decoy_prefix": "rev_",
  "static_mods": {
    "C": 57.0215
  },
  "variable_mods": {
    "M": 15.9949,
    "^": 304.29,
    "K": 304.29
  },
  "fasta": "UP000005640_9606.fasta"
}
EOF

cargo run --release 
```

# API
At the moment, SageRTS only supports a single, JSON-API based endpoint:

```http
POST http://localhost:3000/v1/score/
```

The following query runs a 1000 Da open search round trip in ~5 ms (against the full human tryptic proteome with TMT, MetOx as variable mods)
```json
{
  "precursor_tolerance": {
    "da": [
      -1000,
      1.25
    ]
  },
  "fragment_tolerance": {
    "ppm": [
      -10,
      10
    ]
  },
  "report_psms": 20,
  "chimera": true,
  "deisotope": true,
  "precursor_mz": 1051.0995,
  "precursor_charge": 3,
  "mz": [
    120.08081055,
    123.84158325,
    132.10189819,
    136.11207581,
  ],
  "intensity": [
    1411.77197266,
    928.66247559,
    5548.22070312,
    810.61688232,
  ]
}
```

Which will return a list of PSMs (up to `report_psms`). Note that discriminant scoring and FDR control is not carried out at this time

```json
[
  {
    "peptide": "HHSSETHEVDSDLSYDSSDDSSPSNK",
    "proteins": "rev_sp|P61129|ZC3H6_HUMAN",
    "feature": {
      "peptide_len": 26,
      "spec_id": "real-time",
      "file_id": 0,
      "rank": 1,
      "label": -1,
      "expmass": 3150.2764,
      "calcmass": 2861.1492,
      "charge": 3,
      "rt": 0.0,
      "aligned_rt": 0.0,
      "predicted_rt": 0.0,
      "delta_rt_model": 0.999,
      "delta_mass": 96192.555,
      "isotope_error": 0.0,
      "average_ppm": 6.4890523,
      "hyperscore": 21.440988204515595,
      "delta_next": 0.9691859962323015,
      "delta_best": 0.0,
      "matched_peaks": 5,
      "longest_b": 2,
      "longest_y": 1,
      "longest_y_pct": 0.03846154,
      "missed_cleavages": 0,
      "matched_intensity_pct": 5.2684937,
      "scored_candidates": 31913,
      "poisson": -2.7720455939201067,
      "discriminant_score": 0.0,
      "posterior_error": 1.0,
      "spectrum_q": 1.0,
      "peptide_q": 1.0,
      "protein_q": 1.0,
      "ms2_intensity": 13831.666,
      "ms1_intensity": 0.0
    }
  }
]
```