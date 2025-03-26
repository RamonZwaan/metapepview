"""
This module defines types that contain spectral data.
Spectral data can be derived from mzML, mzXML, or raw formats.

Unused data is filtered out during import.
"""


from .base_classes import DataValidator, DataIO

import pandas as pd
import numpy as np
from typing import IO, Tuple, Self, Dict

from backend.io import mzxml_to_df, mzml_to_df



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
        
    @classmethod
    def from_mzxml_file(cls,
                        mzxml_file: IO[bytes]) -> Tuple[Self, PeaksData]:
        # fields used in data import
        mzxml_fields = [
            'num',
            'msLevel',
            'peaksCount',
            'retentionTime',
            'totIonCurrent',
            'precursorIntensity',
            'precursorMz',
            'peaks',
            'compressionType',
            'byteOrder',
            'precision'
        ]
        
        data = mzxml_to_df(mzxml_file, mzxml_fields)
        
        numeric_fields = [
            'num',
            'msLevel',
            'peaksCount',
            'retentionTime',
            'totIonCurrent',
            'precursorIntensity',
            'precursorMz'
        ]
        peaks_fields = [
            'peaks',
            'compressionType',
            'byteOrder',
            'precision'
        ]
        
        data[numeric_fields] = data[numeric_fields].astype(float)
        
        # split peaks data from spectral data
        peaks_data = data[peaks_fields]
        peaks_data.rename(columns={'compressionType': 'compression type', 
                                   'byteOrder': 'byte order'})
        

        peaks_obj = PeaksData(peaks_data.to_dict('index'))

        return cls(cls.__rename_mzxml_cols(data)), peaks_obj
    
    
    @staticmethod
    def __rename_mzxml_cols(mzxml_df: pd.DataFrame) -> pd.DataFrame:
        """Rename official mzXML field names into custom SpectralData column
        names.

        Args:
            mzxml_df (pd.DataFrame): Wrangled mzXML dataset.

        Returns:
            pd.DataFrame: Dataset with renamed columns
        """
        name_map = {
            "num": "scan number",
            "msLevel":"MS level",
            "peaksCount":"peaks count",
            "retentionTime":"retention time",
            "totIonCurrent":"total ion current",
            "precursorIntensity":"precursor intensity",
            "precursorMz":"precursor mz",
        }
        return mzxml_df.rename(columns=name_map)



class SpectralMetaData:
    def __init__(self, *args, **kwargs):
        ...


def import_mzxml_data(mzxml_file: IO[bytes]) -> Tuple[SpectralData,
                                                      PeaksData,
                                                      SpectralMetaData]:
    """Fully import mzxml file into python data structures. To improve import
    speed and maintain expected behavior of tabular data (identical to pandas),
    the data is separated into three object types. "SpectralData" stores scan
    specific data (excluding peaks and intensities), "PeaksData" stores peaks
    and intensities for each scan, "SpectralMetaData" stores data at MS
    experiment level.
    
    Note:
        File can be imported directly from class constructors. However, this is
        slower as the file will be repeatedly parsed.

    Args:
        mzxml_file (IO[bytes]): File-buffer object of mzxml file (object
            returned when executing "open()")
    Returns:
        Tuple[SpectralData, PeaksData, SpectralMetaData]: _description_
    """
    # fields used in data import
    mzxml_fields = [
        "num",
        "msLevel",
        "peaksCount",
        "retentionTime",
        "totIonCurrent",
        "precursorIntensity",
        "precursorMz",
        "peaks",
        "compressionType",
        "byteOrder",
        "precision",
    ]
    
    data = mzxml_to_df(mzxml_file, mzxml_fields)
    
    numeric_fields = [
        'num',
        'msLevel',
        'peaksCount',
        'retentionTime',
        'totIonCurrent',
        'precursorIntensity',
        'precursorMz'
    ]
    
    data[numeric_fields] = data[numeric_fields].astype(float)
    
    # store peaks inside separate dataset
    peaks_data = data["peaks"]
    peaks_data = peaks_data.to_json()
    data = data.drop(labels="peaks", axis=1)
    
    # create metadata dict
    metadata = dict()
    
    # add peaks compression parameters
    const_fields = ["compressionType", "byteOrder", "precision"]
    for const_field in const_fields:
        metadata[const_field] = data.loc[0, const_field]
    data.drop(const_fields, axis=1, inplace=True)
    
    metadata["scanCount"] = data.shape[0]
    metadata["totalRetentionTime"] = data.iloc[-1]['retentionTime']
    
    # rename columns to consistent SpectralData format
    data = SpectralData.__rename_mzxml_cols(data)
    
    return (SpectralData(data),
            PeaksData(peaks_data),
            SpectralMetaData(metadata))
