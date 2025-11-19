import re
import io
import numpy as np
import pandas as pd

from pathlib import Path
from typing import Type, Self, Sequence, IO, List

from metapepview.backend.types.proteomics.proteomics_base_classes import DbSearchMethods, DeNovoMethods
from metapepview.backend.types.metapep_table import MetaPepDbSearch, MetaPepDeNovo
from metapepview.backend.utils import filter_crap, memory_to_stringio, wrangle_peptides


class NovorDeNovo(DeNovoMethods):
    """Novor de novo table format, as supplied by SearchGUI, not to be confused
    with the latest Novor tool accessed online. Below is a summary of expected
    fields for this format.
    
    Note:
        While scan number is a parameter within the format, data is not
        correctly displayed (all values are zero). For peptide-to-spectrum
        matching, retention time matching can be used as a fallback strategy as
        long as the spectral file contains serial spectra with unique retention
        times.
    
    Fields:
        id
        scanNum
        RT
        mz(data)
        z
        pepMass(denovo)
        err(data-denovo)
        ppm(1e6*err/(mz*z))
        score
        peptide
        aaScore
    """
    
    # expected columns within input data
    REQUIRED_FIELDS = ["id", "scanNum", "RT", "mz(data)", "z", "pepMass(denovo)",
                       "err(data-denovo)", "ppm(1e6*err/(mz*z))", "score",
                       "peptide", "aaScore", "Source File"]
    
    # columns expected to be numeric
    NUMERIC_FIELDS = ["id", "scanNum", "RT", "mz(data)", "z", "pepMass(denovo)",
                      "err(data-denovo)", "ppm(1e6*err/(mz*z))", "score"]
    
    DATA_FORMAT = 'Novor'
    CONFIDENCE_FORMAT = 'Score'
    
    # match everything after '# input file = '
    SOURCE_FILE_PATTERN = re.compile(r"(?<=^# input file = ).+$")
    TABLE_HEADER_PATTERN = re.compile(r"^# id, ")

    def __init__(self,
                 data: pd.DataFrame,
                 file_name: str | None = None):
        success, msg = self.validate_input(data)

        if success is False:
            raise ValueError(msg)
        
        data = data[self.REQUIRED_FIELDS]

        self._data = data
        self._file_name = file_name
   
    @property
    def data(self) -> pd.DataFrame:
       return self._data
    
    @property
    def file_name(self) -> str | None:
       return self._file_name


    @classmethod
    def read_file(cls: Type[Self], path: str | Path) -> Self:
        """Read input file and return data as instance of class object.
        
        Note:
            Source file is not included in data table, but as metadata row on
            top of the input file. During import, the source file is added as
            separate column in the data table.

        Args:
            cls (Type[T]): Class type
            path (str | Path): Location of file

        Returns:
            T: Instance of class object
        """
        file_name = Path(path).stem
        source_file = cls.get_source_file(path)
        
        df = pd.read_csv(path, skiprows=20)
        df['Source File'] = source_file
        
        cls_obj = cls(df, file_name) 

        return cls_obj


    @classmethod
    def read_file_buffer(cls: Type[Self],
                         file_buffer: str | IO[str],
                         file_name: str | None = None) -> Self:
        """Read input file buffer and return data instance of class object.
        Use this when data is imported into memory as string buffer.
        
        >>> Class.read_file("example.csv")
        Is equal to
        >>> Class.read_file_buffer(open("example.csv"))

        Note:
            Source file is not included in data table, but as metadata row on
            top of the input file. During import, the source file is added as
            separate column in the data table.

        Args:
            cls (Type[T]): Class type
            file_buffer (str | IO[str]): input file data or buffer
            file_name (str | None, optional): Name of input file. 
                Defaults to None.

        Returns:
            T: Instance of class object
        """
        # get source file
        source_file = cls.get_source_file_buffer(file_buffer)
        
        # determine start of data table
        skip_rows = cls.get_table_start_row(file_buffer)
        
        # if data has been decoded, directly return csv DataFrame
        if isinstance(file_buffer, io.TextIOBase):
            df = pd.read_csv(file_buffer, 
                             skiprows=skip_rows, 
                             delimiter=', ')
        elif isinstance(file_buffer, str):
            df = pd.read_csv(memory_to_stringio(file_buffer),
                             skiprows=skip_rows,
                             delimiter=', ')
        else:
            raise TypeError("invalid content type supplied...")
        
        df['Source File'] = source_file
        
        # remove comment prefix and traling comma for header
        df = df.rename(columns={"# id": "id",
                                "aaScore,": "aaScore"})
        
        cls_obj = cls(df, file_name) # type: ignore
        
        return cls_obj
    

    @classmethod
    def get_source_file(cls, file_name: Path | str) -> str:
        """Return spectral file name from novor path.

        Args:
            file_name (Path | str): Location of novor file.

        Returns:
            str: Name of spectral file
        """
        with open(file_name) as contents:
            return cls.get_source_file_buffer(contents)
    

    @classmethod
    def get_source_file_buffer(cls, file_buffer: IO[str] | str) -> str:
        """Return spectral file name from novor file stream buffer.

        Args:
            file_buffer (Path | str): novor file stream buffer.

        Returns:
            str: Name of spectral file
        """
        if isinstance(file_buffer, str):
            file_buffer = memory_to_stringio(file_buffer)
        
        source = None
        # fetch metadata content, iterate until source file fetched
        c = 0
        while source == None:
            row = file_buffer.readline()
            source_match = re.search(cls.SOURCE_FILE_PATTERN, row)
        
            if source_match is not None:
                source = source_match.group(0).rstrip()
            
            # only parse first 20 rows 
            c += 1
            if c > 20:
                break
            
        if source is None:
            raise ValueError("unable to fetch source file from Novor data")
        
        # reset buffer to start of stream to allow subsequent functions to process same stream object
        file_buffer.seek(0)
        
        # return source
        return Path(source).stem
    
    
    @classmethod
    def get_table_start_row(cls, file_buffer: IO[str] | str) -> int:
        """Return line number where data table starts.

        Args:
            file_buffer (Path | str): novor file stream buffer.

        Returns:
            int: row number of header names
        """
        if isinstance(file_buffer, str):
            file_buffer = memory_to_stringio(file_buffer)
        
        row_num = None
        # fetch metadata content, iterate until table headers located
        idx = 0
        while row_num == None:
            row = file_buffer.readline()
            row_match = re.search(cls.TABLE_HEADER_PATTERN, row)
        
            if row_match is not None:
                row_num = idx
                break
            
            # only parse first 100 rows 
            idx += 1
            if idx > 100:
                break
            
        if row_num is None:
            raise ValueError("unable to locate data table")
        
        # reset buffer to start of stream to allow subsequent functions to process same stream object
        file_buffer.seek(0)
        
        # return source
        return row_num
    

    def get_source_files(self) -> Sequence[str]:
        """Return all raw spectral file names from dataset, excluding file type
        suffix.

        Returns:
            Sequence[str]: All raw spectral file names in dataset.
        """
        source_file_col = self.data['Source File']
        source_files: List[str] =  source_file_col\
            .dropna()\
            .drop_duplicates()\
            .apply(lambda x: Path(x).stem)\
            .tolist()

        return source_files


    def to_metapep_de_novo(self,
                           sample_name: str | None = None,
                           crap_dataset: pd.Series | None = None) -> MetaPepDeNovo:
        """Convert data table to MetaPepDeNovo format.

        Args:
            sample_name (str | None, optional): Name of sample to pass to
                MetaPepDeNovo object, overwrites file name specified in object.
                Defaults to None.
            crap_dataset (pd.Series | None, optional): Pandas series of peptide
                sequences part of the cRAP dataset. These are filtered out of
                the de novo peptide data. defaults to None.

        Returns:
            MetaPepDeNovo: De novo data table in Metapep de novo format.
        """
        # drop unrelevant columns and rename columns to metapep format
        df = self.data[["scanNum", "RT", "mz(data)", "z", "pepMass(denovo)",
                        "ppm(1e6*err/(mz*z))", "score", "peptide", "Source File"]]\
            .rename(columns={'peptide': 'Peptide',
                             'scanNum': 'Scan',
                             'mz(data)': 'm/z',
                             'z': 'Charge',
                             'ppm(1e6*err/(mz*z))': 'ppm',
                             'score': 'Confidence',
                             'pepMass(denovo)': 'Mass'})
        
        # Wrangle sequence into consistent format (remove PTM, equalte L, I)
        df.loc[:, 'Sequence'] = df.loc[:, 'Peptide'].apply(wrangle_peptides)
        df.loc[:, 'Length'] = df.loc[:, 'Sequence'].apply(len)
        
        # remove file type suffix from Source File column
        df.loc[:, 'Source File'] = df.loc[:, 'Source File'].apply(
            lambda x: Path(x).stem
        )

        # Signal intensities are not provided by Casanovo
        if df.shape[0] > 0:
            df.loc[:, "Area"] = np.nan
            df.loc[:, "PTM"] = np.nan
        else:
            df.loc[:, "Area"] = []
            df.loc[:, "PTM"] = []

        # filter out peptides from cRAP dataset
        if crap_dataset is not None:
            df = filter_crap(df, 'Sequence', crap_dataset)

        # fetch all spectral file names present in dataset
        source_files: List[str] = df.loc[:, 'Source File'].dropna().unique().tolist()
        
        # configure file name, if given the function argument, else the class file name
        if sample_name is None:
            sample_name = self.file_name
        
        # Put dataset into MetaPepDbSearch class object
        return MetaPepDeNovo(df,
                             self.DATA_FORMAT,
                             self.CONFIDENCE_FORMAT,
                             source_files,
                             sample_name)
        