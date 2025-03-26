
from typing import Tuple
import numpy as np
import pandas as pd

from backend.utils import *



def wrangle_psm_dataset(psm_data: pd.DataFrame,
                        peptide_col: str='Peptide',
                        area_col: str='Area'
                        ) -> pd.DataFrame:
    """Process data to consistent and expected format.
    """
    # Wrangle db search psm dataset
    psm_data.loc[:, 'Sequence'] = psm_data[peptide_col].transform(wrangle_peptides)
    psm_data[area_col] = psm_data[area_col].fillna(0.0)
    return psm_data


def validate_psm_peaks_11(psm_data: pd.DataFrame) -> Tuple[bool, str | None]:
    """Check if input data follows expected format.
    The following criteria are tested:
        - File contains columns as expexted from format
        - correct columns are numeric
        
    Args:
        psm_data (pd.DataFrame): peaks 11 format psm file

    Returns:
        Tuple[bool, str | None]: Whether input data satisfies criteria, and 
            error message if issue encountered
    """
    # columns expected from psm input data
    expected_cols = ["Peptide", "-10LgP", "Mass", "Length", "ppm", "m/z",
                     "z", "RT", "Scan", "Area", "Feature Id", "Source File",
                     "Accession", "PTM", "AScore"]
    
    # columns expected to be numeric
    numeric_cols = ["-10LgP", "Mass", "Length", "ppm", "m/z",
                     "z", "RT", "Scan", "Area", "Feature Id"]
    
    # check that all expected columns are present within the file
    if not all(i in psm_data.columns for i in expected_cols):
        err_msg = "Invalid format, missing expected columns in dataset"
        return (False, err_msg)
    
    # check if expected columns are numeric
    psm_dtypes = psm_data.dtypes
    if not all(pd.api.types.is_numeric_dtype(psm_data.loc[:, i]) 
               for i in numeric_cols):
        err_msg = "Invalid format, non-numeric values encountered in numeric columns"
        return (False, err_msg)
    
    return (True, None)
        

def validate_de_novo_peaks_11(de_novo_data: pd.DataFrame) -> Tuple[bool, str | None]:
    """Check if de novo input data follows expected format.
    The following criteria are tested:
        - File contains columns as expexted from format
        - correct columns are numeric
        
    Args:
        de_novo_data (pd.DataFrame): peaks 11 format de novo file

    Returns:
        Tuple[bool, str | None]: Whether input data satisfies criteria, and 
            error message if issue encountered
    """
    # columns expected from psm input data
    expected_cols = ["Source File", "Scan", "Peptide", "Tag length",
                     "ALC (%)", "Length", "m/z", "z", "RT", "Area", "Mass",
                     "ppm", "PTM", "local confidence (%)", "mode",
                     "tag(>=0.0%)", "Feature Id"]
    
    # columns expected to be numeric
    numeric_cols = ["Scan", "Tag length", "ALC (%)", "Length", "m/z", "z",
                    "RT", "Area", "Mass", "ppm", "Feature Id"]
    
    # check that all expected columns are present within the file
    if not all(i in de_novo_data.columns for i in expected_cols):
        err_msg = "Invalid format, missing expected columns in dataset"
        return (False, err_msg)
    
    # check if expected columns are numeric
    psm_dtypes = de_novo_data.dtypes
    if not all(pd.api.types.is_numeric_dtype(de_novo_data.loc[:, i]) 
               for i in numeric_cols):
        err_msg = "Invalid format, non-numeric values encountered in numeric columns"
        return (False, err_msg)
    
    return (True, None)
