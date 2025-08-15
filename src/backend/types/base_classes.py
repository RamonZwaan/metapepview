"""
This module defines base classes that implement shared behavior for the
type classes defined elsewhere.
"""

from typing import List, Tuple, TypeVar, Type, IO, Self, Protocol
from abc import ABC, abstractmethod
from pathlib import Path
import io
from functools import wraps

import pandas as pd
import numpy as np

from backend.utils import memory_to_stringio
from .definitions import *


class DataValidator(Protocol):
    """Implements validator function that checks fields within input data
    for consistency with field names described in class definitions
    """
    
    REQUIRED_FIELDS: List[str]
    NUMERIC_FIELDS: List[str]
    
    @classmethod
    def validate_input(cls,
                       input_df: pd.DataFrame,
                       optional_fields: List[str] | None = None) -> Tuple[bool, str | None]:
        """Check if input data follows expected format.
        The following criteria are tested:
            - File contains columns as expexted from format
            - correct columns are numeric
        
        Args:
            input_df (pd.DataFrame): Input data.
            optional_fields (List[str] | None, optional): Fields part of the
                input_df but not required to be present. Defaults to None.
        
        Returns:
            Tuple[bool, str | None]: Whether input data satisfies criteria, and 
            error message if issue encountered
        """
        # check that all expected columns are present within the file
        field_array = np.array(cls.REQUIRED_FIELDS)
        missing_cols = field_array[~np.isin(field_array, input_df.columns)]
        if len(missing_cols) > 0:
            err_msg = "Invalid format, missing expected columns in dataset: {}".format(missing_cols)
            return (False, err_msg)
        
        # get columns that are expected numeric but object
        not_numeric = []
        for col in cls.NUMERIC_FIELDS:
            if not pd.api.types.is_numeric_dtype(input_df.loc[:, col]):
                not_numeric.append(col)
        
        # check if unexpected non-numeric columns encountered
        if len(not_numeric) > 0:
            # Convert ambiguous representation of missing values within numeric columns, such as '-' or ''
            input_df.loc[:, not_numeric] = input_df.loc[:, not_numeric]\
                .replace(["", "-"],
                        np.nan)
            
            # try to convert columns to floats, this may fail if other ambiguous strings are present
            try:
                input_df.loc[:, not_numeric] = input_df.loc[:, not_numeric].astype(float)
            except ValueError:
                err_msg = "Invalid format, non-numeric values encountered in numeric columns"
                return (False, err_msg)
        
        return (True, None)


T = TypeVar('T', bound='DataIO')

class DataIO(Protocol):
    """Provide wrapper over pandas IO functions that return dataframe in
    correct class objects.
    """
    
    
    @classmethod
    def read_json(cls: Type[T], path: str | Path) -> T:
        """Read json file and return data as instance of class object.

        Args:
            cls (Type[T]): Class type
            path (str | Path): Location of json file

        Returns:
            T: Instance of class object
        """
        file_name = Path(path).stem
        df = pd.read_json(path)
        
        cls_obj = cls(df, file_name) # type: ignore

        return cls_obj
        
    @classmethod
    def import_pandas(cls: Type[T], *args, **kwargs) -> T:
        """Import data into class object through pandas DataFrame
        constructor interface.

        Args:
            cls (Type[T]): Class type.

        Returns:
            T: Instance of class object
        """
        df = pd.DataFrame(*args, **kwargs)
        return cls(df) # type: ignore
    

    @classmethod
    def _read_csv(cls: Type[T], path: str | Path, delim: str | None=None) -> T:
        """Read csv file and return data as instance of class object.

        Args:
            cls (Type[T]): Class type
            path (str | Path): Location of csv file
            delim (str | None, Optional): Custom delimiter for import of tsv
                or other types of table files. If not specified, then regular
                csv delimiter used Defaults to None

        Returns:
            T: Instance of class object
        """
        file_name = Path(path).stem
        df = pd.read_csv(path, delimiter=delim, low_memory=False)
        
        cls_obj = cls(df, file_name) # type: ignore

        return cls_obj
    

    @classmethod
    def _read_csv_buffer(cls: Type[T],
                         file_buffer: str | IO[str],
                         file_name: str | None = None,
                         delim: str | None = None) -> T:
        """Read csv file buffer and return data instance of class object.
        Use this when data is imported into memory as str string buffer.
        
        >>> Class.__read_csv("example.csv")
        Is equal to
        >>> Class.__read_csv_buffer(open("example.csv"))

        Args:
            cls (Type[T]): Class type
            file_buffer (str | IO[str]): csv file data or buffer
            file_name (str | None, optional): Name of csv file. 
                Defaults to None.
            delim (str | None, Optional): Custom delimiter for import of tsv
                or other types of table files. If not specified, then regular
                csv delimiter used Defaults to None

        Returns:
            T: Instance of class object
        """
        # if data has been decoded, directly return csv DataFrame
        if isinstance(file_buffer, io.TextIOBase):
            df = pd.read_csv(file_buffer, delimiter=delim, low_memory=False)
        elif isinstance(file_buffer, str):
            df = pd.read_csv(memory_to_stringio(file_buffer), 
                             delimiter=delim,
                             low_memory=False)
        else:
            raise TypeError("invalid content type supplied...")
        
        cls_obj = cls(df, file_name) # type: ignore
        
        return cls_obj
    

    @classmethod
    def read_file(cls: Type[T], path: str | Path) -> T:
        """Read input file and return data as instance of class object.
        Overwrite method with function that adheres to input format.

        Args:
            cls (Type[T]): Class type
            path (str | Path): Location of file

        Returns:
            T: Instance of class object
        """
        # default behavior reads csv format
        return cls._read_csv(path)


    @classmethod
    def read_file_buffer(cls: Type[T],
                         file_buffer: str | IO[str],
                         file_name: str | None = None) -> T:
        """Read input file buffer and return data instance of class object.
        Use this when data is imported into memory as string buffer. Overwrite
        method with function that adheres to input format.
        
        >>> Class.read_file("example.csv")
        Is equal to
        >>> Class.read_file_buffer(open("example.csv"))

        Args:
            cls (Type[T]): Class type
            file_buffer (str | IO[str]): input file data or buffer
            file_name (str | None, optional): Name of input file. 
                Defaults to None.

        Returns:
            T: Instance of class object
        """
        # default behavior reads csv buffer
        return cls._read_csv_buffer(file_buffer, file_name)


