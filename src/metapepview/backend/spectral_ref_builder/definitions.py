from pathlib import Path
from typing import Dict, List
from metapepview.backend import DbSearchSource, DeNovoSource
from dataclasses import dataclass


db_search_file_name: Dict[DbSearchSource, str] = {
    'Peaks 11': "db.psms.csv",
    'Peaks 10': "DB search psm.csv",
    'MaxQuant': "evidence.txt",
    # 'ProteomeDiscoverer': ...,
    'Sage': "*.sage.tsv",
}

de_novo_file_name: Dict[DeNovoSource, str] = {
    'Peaks 11': "*.denovo.csv",
    'Peaks 10': "de novo peptides.csv",
    'Novor': "*.novor.csv",
    'Casanovo': "*.mztab"
}


@dataclass
class RefBuilderOptions():
    root_dir: Path
    db_search_format: DbSearchSource
    db_search_file_pattern: str
    de_novo_format: DeNovoSource
    de_novo_file_pattern: str
    db_search_thresholds: List[float]
    de_novo_thresholds: List[float]
    intensity_percentiles: List[int | float]
    transmission_loss_percentiles: List[int | float]