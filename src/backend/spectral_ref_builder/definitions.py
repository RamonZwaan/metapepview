from pathlib import Path
from typing import Dict, List
from backend import DbSearchSource, DeNovoSource
from dataclasses import dataclass


db_search_file_name: Dict[DbSearchSource, str] = {
    'Peaks': "db.psms.csv",
    'Peaks10': "DB search psm.csv",
    'MaxQuant': "evidence.txt",
    # 'ProteomeDiscoverer': ...,
    'Sage': "*.sage.tsv",
}

de_novo_file_name: Dict[DeNovoSource, str] = {
    'Peaks': "*.denovo.csv",
    'Peaks10': "de novo peptides.csv",
    'Novor': "*.novor.csv"
}


@dataclass
class RefBuilderOptions():
    root_dir: Path
    db_search_format: DbSearchSource
    db_search_file_pattern: str
    de_novo_format: DeNovoSource
    de_novo_file_pattern: str
    db_search_thresholds: List[int | float]
    de_novo_thresholds: List[int | float]
    intensity_percentiles: List[int | float]
    transmission_loss_percentiles: List[int | float]