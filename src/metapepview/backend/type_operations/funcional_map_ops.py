import base64
import io
from typing import List, Tuple

import numpy as np
import pandas as pd

from metapepview.backend.utils import *
from metapepview.backend.types.object_mappings import functional_db_importers
from metapepview.backend.types.definitions import FuncAnnotFormat
from metapepview.backend.types import FunctionDbMapper


def import_func_map(upload_contents: str,
                    db_format: FuncAnnotFormat,
                    accession_pattern: str | None = None,
                    max_evalue = 1e-6,
                    archive_format: str | None = None) -> FunctionDbMapper:
    """Import functional annotation mapper object based on format settings
    supplied.
    """
    importer = functional_db_importers.get(db_format)
    
    if importer is None:
        raise ValueError("Invalid file format supplied")
    else:
        # extract file if delivered in archive    
        file_buffer = memory_to_stringio(upload_contents, archive_format)
        
        func_mapper_obj = importer.read_file_buffer(file_buffer, 
                                                    max_evalue=max_evalue,
                                                    acc_regex=accession_pattern)
        
        return func_mapper_obj


def validate_func_map(file_buffer: str | IO[str],
                      file_format: FuncAnnotFormat,
                      archive_format: str | None = None) -> Tuple[bool, str | None]:
    """Check if imported functional annotation dataset adheres to expected formats

    Args:
        file_buffer (str | IO[str]): String contents or string buffer of dataset.
        file_format (FuncAnnotFormat): Functional annotation format.
        archive_format (str | None, optional): Archive format of data if compressed.
            Defaults to None.

    Returns:
        Tuple[bool, str | None]: Success status with error message if failed.
    """
    importer = functional_db_importers.get(file_format)
    
    if importer is None:
        return False, "Invalid file format supplied"
    else:
        try:
            if isinstance(file_buffer, str):
                file_buffer = memory_to_stringio(file_buffer, archive_format)
            importer.read_file_buffer(file_buffer)
            return True, None
        except Exception as e:
            return False, repr(e)
