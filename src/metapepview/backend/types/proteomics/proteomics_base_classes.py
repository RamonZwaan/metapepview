"""
This module defines base classes that implement behavior of proteomics
classes.
"""

from typing import Protocol
from abc import abstractmethod
from pathlib import Path

import pandas as pd

from metapepview.backend.types.definitions import *
from metapepview.backend.types.base_classes import *
from metapepview.backend.types.metapep_table import MetaPepDbSearch, MetaPepDeNovo


class DbSearchMethods(DataIO, DataValidator, Protocol):
    
    ACCESSION_DELIMITER: str
    DATA_FORMAT: DbSearchSource
    CONFIDENCE_FORMAT: DbSearchConfFormat
    
    @classmethod
    @abstractmethod
    def get_source_file(cls, file_name: Path | str) -> str:
        ...
        
    @abstractmethod
    def get_source_files(self) -> List[str]:
        ...

    @abstractmethod
    def to_metapep_db_search(self,
                             sample_name: str | None = None,
                             crap_dataset: pd.Series | None = None) -> MetaPepDbSearch:
        ...
    
    @abstractmethod
    def data(self) -> pd.DataFrame:
       ...
    

class DeNovoMethods(DataIO, DataValidator, Protocol):
    
    DATA_FORMAT: DeNovoSource
    CONFIDENCE_FORMAT: DeNovoConfFormat
    
    @classmethod
    @abstractmethod
    def get_source_file(cls, file_name: Path | str) -> str:
        ...
    
    @abstractmethod
    def to_metapep_de_novo(self, 
                           sample_name: str | None = None,
                           crap_dataset: pd.Series | None = None) -> MetaPepDeNovo:
        ...
    
    @abstractmethod
    def data(self) -> pd.DataFrame:
       ...
    