from dataclasses import dataclass
from typing import Dict, List, Optional
import requests
import matplotlib.pyplot as plt


@dataclass
class MatchedPeak:
    mz: float
    intensity: float
    fragment_mz: float
    fragment_kind: str
    fragment_idx: float


@dataclass
class Score:
    peptide: str
    peptide_len: int
    proteins: str
    scannr: int
    label: int
    expmass: float
    calcmass: float
    charge: int
    rt: float
    delta_mass: float
    hyperscore: float
    delta_hyperscore: float
    matched_peaks: int
    longest_b: int
    longest_y: int
    matched_intensity_pct: float
    scored_candidates: int
    poisson: float
    q_value: float
    specid: int


@dataclass
class Spectrum:
    scan: int
    monoisotopic_mass: float
    charge: int
    rt: float
    mz: List[float]
    intensity: List[float]


class Api:
    def __init__(self):
        self.url = "http://localhost:3000/spectrum"

    def get_spectrum(self, scan_id: int) -> Optional[Spectrum]:
        resp = requests.get(
            f"{self.url}/{scan_id}?deisotope=true", headers={"accept-encoding": "gzip"}
        )
        if resp.ok:
            return Spectrum(**resp.json())
        else:
            print(resp.text)

    def get_matches(self, scan_id: int, data) -> List[Score]:
        resp = requests.post(
            f"{self.url}/{scan_id}", json=data, headers={"accept-encoding": "gzip"}
        )
        if resp.ok:
            return [Score(**x) for x in resp.json()]
        else:
            print(resp.text)

    def get_matched_peaks(
        self,
        scan_id: int,
        peptide: str,
        modifications: Dict[str, float],
        ppm: float = 10,
        deisotope: bool = False,
    ) -> List[MatchedPeak]:
        data = {
            "sequence": ''.join(c for c in peptide if c.isalpha()),
            "modifications": modifications,
            "fragment_tol": {"ppm": [-ppm, ppm]},
            "deisotope": deisotope,
        }

        resp = requests.post(
            f"{self.url}/{scan_id}/peptide",
            json=data,
            headers={"accept-encoding": "gzip"},
        )
        if resp.ok:
            return [MatchedPeak(**x) for x in resp.json()]
        else:
            print(resp.text)


class SpectrumGraph:
    def __init__(self, scan_id: int):
        self.api = Api()
        self.scan_id = scan_id
        self.spectrum = self.api.get_spectrum(scan_id)
        self.anns = []
        self.matches = self.api.get_matches(
            scan_id,
            {
                "precursor_tolerance": {"da": [-3.6, 1.2]},
                "fragment_tolerance": {"ppm": [-10.0, 10.0]},
                "report_psms": 10,
                "chimera": True,
                "deisotope": True,
            },
        )
    
    def base_plot(self):
        self.fig, self.ax = plt.subplots()

        for mz, intensity in zip(self.spectrum.mz, self.spectrum.intensity):
            self.ax.vlines(mz, 0, intensity, colors="grey", linewidths=1, alpha=0.8)

    def plot_match(self, match: Score):
        peaks = self.api.get_matched_peaks(self.scan_id, match.peptide, dict(C=57.0215), deisotope=True)
        for peak in peaks:
            color = '#1976D2' if peak.fragment_kind == 'B' else '#D32F2F'
            self.ax.vlines(peak.mz, 0, peak.intensity, colors=color, linewidths=1)
            self.anns += [self.ax.annotate(f"{peak.fragment_kind.lower()}{peak.fragment_idx}+", xy=(peak.mz, peak.intensity+20), rotation=90, color=color)]
    
    def plot(self):
        plt.show()

# peaks = Api().get_matched_peaks(30069, matches[1].peptide, dict(C=57.0215))
# for peak in peaks:
#     ax.vlines(peak.mz, 0, peak.intensity, colors='blue', linewidths=1)
#     ax.annotate(f"{peak.fragment_kind.lower()}{peak.fragment_idx}+", xy=(peak.mz, peak.intensity), rotation=90)


s = SpectrumGraph(30091)
s.base_plot()
for m in s.matches[:2]:
    s.plot_match(m)
    print(m)
    # plt.suptitle(m.peptide)
plt.show()
plt.close()