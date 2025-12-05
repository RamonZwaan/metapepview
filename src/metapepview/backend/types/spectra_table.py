"""
This module defines types that contain spectral data.
Spectral data can be derived from mzML, mzXML, or raw formats.

Unused data is filtered out during import.
"""


from metapepview.backend.types.base_classes import DataValidator, DataIO

import pandas as pd
import numpy as np
from typing import IO, Tuple, Self, Dict

from metapepview.backend.io import mzml_to_df



class PeaksData:
    """Custom table of MS spectral data. Rows represent separate spectra.
    
    
    Fields:
        peaks                   # Binary encoded peaks list (Only mzXML)
        peaks mz                # binary encoded peak mz list (only mzML)
        intensities             # Binary encoded peak intensities list (only mzML)
        compression type        (Optional)
        byte order              (Optional, not in mzML format)
        precision               (Optional)
    """
    def __init__(self,
                 peaks_data: Dict[str, str]):
        # contains file_format, mzml or mzxml, compression properties and peaks data
        self._peaks_data = peaks_data


    @classmethod
    def from_json(cls) -> Self:
        ...

    
    @classmethod
    def from_mzxml(cls) -> Self:
        ...


    @classmethod
    def from_mzml(cls) -> Self:
        ...
    

    def fetch_peak_data(self) -> Tuple[np.ndarray, np.ndarray]:
        ...


    def store_json(self) -> str:
        ...


class SpectralData(DataValidator, DataIO):
    """Custom table of MS spectral data. Rows represent separate spectra.
    
    Fields:
        scan number
        MS level
        peaks count
        retention time
        total ion current
        injection time          (Not in mzXML format)
        precursor intensity
        precursor mz
    """    

    # expected columns within input data
    REQUIRED_FIELDS = [
        "scan number",
        "MS level",
        "peaks count",
        "retention time",
        "ion injection time"
        "total ion current", 
        "precursor intensity",
        "precursor mz"
        "precursor charge"
        "precursor scan number"
        # "m/z array",
        # "intensity array"
    ]
    
    # columns expected to be numeric
    NUMERIC_FIELDS = REQUIRED_FIELDS
    
    
    def __init__(self, data: pd.DataFrame):
        # 'injection time' is not present within mzxml format, then add empty column
        if "ion injection time" not in data.columns:
            data["ion injection time"] = np.nan

        success, msg = self.validate_input(data)

        if success is False:
            raise ValueError(msg)
        
        self._data = data
    

    @property
    def data(self) -> pd.DataFrame:
       return self._data


    @classmethod
    def from_df(cls,
                df: pd.DataFrame) -> Self:
        return cls(df)


    @classmethod
    def from_mzml(cls,
                  mzml_file: IO[bytes]) -> Self:
        data, metadata = mzml_to_df(mzml_file)
        return cls(data)


class SpectralMetaData:
    def __init__(self, *args, **kwargs):
        ...
