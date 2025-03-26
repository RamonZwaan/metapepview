from pathlib import Path

from backend.spectral_ref_builder import build_ref_data

build_ref_data(
    Path("O:/SIHUMIx"),
    Path("C:/Users/ramonvanderzwa/Repositories/Metapepview/data/share/qc_refs"),
    "sihumix_refs.json",
    "Peaks 10",
    "Peaks 10")