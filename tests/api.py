import json
import requests
import matplotlib.pyplot as plt
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional

PROTON = 1.0072764

@dataclass
class MatchedPeak:
    mz: Optional[float]
    intensity: Optional[float]
    charge: Optional[int]
    fragment_mz: float
    fragment_loss: float
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
class Precursor:
    mz: float
    intensity: Optional[float]
    charge: Optional[int]
    scan: Optional[int]

@dataclass
class Peak:
    mass: float
    intensity: float

@dataclass
class Spectrum:
    level: int
    scan: int
    scan_start_time: float
    ion_injection_time: float
    precursors: List[Precursor]
    peaks: List[Peak]


class Api:
    def __init__(self):
        self.url = "http://localhost:3000/spectrum"

    def get_spectrum(
        self, scan_id: int, deisotope=True, max_peaks=150
    ) -> Optional[Spectrum]:
        resp = requests.get(
            f"{self.url}/{scan_id}?deisotope={str(deisotope).lower()}&max_peaks={max_peaks}",
            headers={"accept-encoding": "gzip"},
        )
        if resp.ok:
            sp = Spectrum(**resp.json())
            sp.peaks = [Peak(**p) for p in sp.peaks]
            sp.precursors = [Precursor(**p) for p in sp.precursors]
            return sp
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
        tolerance: Dict[str, List[float]],
        deisotope: bool = False,
        nterm: Optional[float] = None,
    ) -> List[MatchedPeak]:
        data = {
            "sequence": "".join(c for c in peptide if c.isalpha()),
            "modifications": modifications,
            "fragment_tol": tolerance,
            "deisotope": deisotope,
            "nterm": nterm,
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
    def __init__(
        self,
        scan_id: int,
        tolerance: Dict[str, List[float]],
        modifications: Dict[str, float],
        nterm: Optional[float] = None,
        deisotope=False,
        chimera=False,
        max_peaks=150,
    ):
        self.api = Api()
        self.scan_id = scan_id
        self.spectrum = self.api.get_spectrum(
            scan_id, deisotope=deisotope, max_peaks=max_peaks
        )
        self.deisotope = deisotope
        self.tolerance = tolerance
        self.modifications = modifications
        self.nterm = nterm
        self.fig = None
        self.anns = []
        self.matches = self.api.get_matches(
            scan_id,
            {
                "precursor_tolerance": {"da": [-3.6, 1.2]},
                "fragment_tolerance": tolerance,
                "report_psms": 10,
                "chimera": chimera,
                "deisotope": deisotope,
            },
        )

    def base_plot(self, sign=1):
        if self.fig is None:
            self.fig, self.ax = plt.subplots()
        for peak in self.spectrum.peaks:
            self.ax.vlines(peak.mass + PROTON, 0, sign * peak.intensity, colors="grey", linewidths=1, alpha=0.8)

    def plot_match(self, peptide: str, sign=1):
        peaks = self.api.get_matched_peaks(
            self.scan_id,
            peptide,
            modifications=self.modifications,
            nterm=self.nterm,
            tolerance=self.tolerance,
            deisotope=self.deisotope,
        )
        print(len(peaks))
        for peak in peaks:
            color = "#1976D2" if peak.fragment_kind == "B" else "#D32F2F"
            y = -50
            x = peak.fragment_mz
            if peak.mz and peak.intensity:
                self.ax.vlines(peak.mz, 0, sign * peak.intensity, colors=color, linewidths=1)
                x = peak.mz
                y = peak.intensity * sign
            else:
                self.ax.vlines(
                    peak.fragment_mz, 0, y, colors="purple", linewidths=1
                )
                # color = "purple"

            pluses = ''.join('+' for _ in range(peak.charge))
            loss = ''
            if peak.fragment_loss >= 18:
                loss = '-H2O'
            elif peak.fragment_loss >= 17:
                loss = '-NH3'
            
            self.anns += [
                self.ax.annotate(
                    f"{peak.fragment_kind.lower()}{peak.fragment_idx}{pluses}{loss}",
                    xy=(x, y),
                    rotation=90,
                    color=color,
                )
            ]

    def plot(self):
        plt.show()


# peaks = Api().get_matched_peaks(30069, matches[1].peptide, dict(C=57.0215))
# for peak in peaks:
#     ax.vlines(peak.mz, 0, peak.intensity, colors='blue', linewidths=1)
#     ax.annotate(f"{peak.fragment_kind.lower()}{peak.fragment_idx}+", xy=(peak.mz, peak.intensity), rotation=90)


s = SpectrumGraph(
    40324,
    tolerance=dict(da=[-0.3, 0.3]),
    deisotope=False,
    chimera=False,
    max_peaks=100,
    modifications=dict(C=57.0215, K=304.20715),
    nterm=304.20715
)
# ann = 'GVLGYTEDAVVSSDFLGDSNSSIFDAAAGIQLSPK'
# ann = 'QLTFEEIAK'
ann = 'CAVVSSAGSLK'
s.base_plot()
s.base_plot(sign=-1)
s.ax.hlines(0, 0, 1200, colors='black', linewidths=2)
s.plot_match(ann)
s.plot_match(s.matches[0].peptide, sign=-1)
plt.suptitle(f"{ann}, {s.matches[0].peptide}")
plt.tight_layout()
plt.show()
plt.close()
for m in s.matches[:3]:
    print(json.dumps(asdict(m), indent=2))
