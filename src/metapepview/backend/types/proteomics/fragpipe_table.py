from __future__ import annotations

import pandas as pd
import numpy as np

from pathlib import Path
from typing import List, Sequence, Type, Self, IO

from metapepview.backend.types.proteomics.proteomics_base_classes import DbSearchMethods
from metapepview.constants import *
from metapepview.backend.types.metapep_table import MetaPepDbSearch
from metapepview.backend.utils import filter_crap, wrangle_peptides, mz_diff_to_ppm


class FragPipeDbSearch(DbSearchMethods):
    """FragPipe database search psm table.
    Below is a summary of table columns when default parameters are specified
    for DB search. Some columns may have been supplemented from post processing
    modules within FragPipe.

    Hyperscore is taken as confidence parameter.
    
    Fields:
    Spectrum	
    Spectrum File	
    Peptide	
    Modified Peptide	
    Extended Peptide	
    Prev AA	
    Next AA	
    Peptide Length	
    Charge	
    Retention	
    Observed Mass	
    Calibrated Observed Mass	
    Observed M/Z	
    Calibrated Observed M/Z	
    Calculated Peptide Mass	
    Calculated M/Z	
    Delta Mass	
    SpectralSim	
    RTScore	
    Expectation	
    Hyperscore	
    Nextscore	
    Probability	
    Qvalue	
    Number of Enzymatic Termini	
    Number of Missed Cleavages	
    Protein Start	
    Protein End	
    Intensity	
    Assigned Modifications	
    Observed Modifications	
    Purity	
    Is Decoy	
    Is Contaminant	
    Is Unique	
    Protein	
    Protein ID	
    Entry Name	
    Gene	
    Protein Description	
    Mapped Genes	
    Mapped Proteins
    """
    
    # expected columns within input data
    REQUIRED_FIELDS = ["Peptide", "Modified Peptide", "Intensity",
                       "Hyperscore", "Observed Mass", "Observed M/Z",
                       "Calculated M/Z", "Peptide Length", "Delta Mass", 
                       "Charge", "Retention", "Spectrum", "Spectrum File", 
                       "Protein", "Is Decoy", "Is Contaminant", "Mapped Proteins"]
    
    # columns expected to be numeric
    NUMERIC_FIELDS = ["Hyperscore", "Intensity", "Observed Mass", 
                      "Observed M/Z", "Calculated M/Z", "Peptide Length", 
                      "Delta Mass", "Charge", "Retention"]
    
    ACCESSION_DELIMITER = ', '

    DATA_FORMAT = 'FragPipe'
    CONFIDENCE_FORMAT = 'Hyperscore'

    ONLY_RAZOR_PROTEIN = False
    
    
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

        Args:
            cls (Type[T]): Class type
            path (str | Path): Location of file

        Returns:
            T: Instance of class object
        """
        # Execute csv reader with tab delimiter from DbSearchMethods
        return cls._read_csv(path, delim='\t')


    @classmethod
    def read_file_buffer(cls: Type[Self],
                         file_buffer: str | IO[str],
                         file_name: str | None = None) -> Self:
        """Read input file buffer and return data instance of class object.
        Use this when data is imported into memory as string buffer.
        
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
        return cls._read_csv_buffer(file_buffer, file_name, delim='\t')


    @classmethod
    def get_source_file(cls,
                        file_name: Path | str,
                        delimiter: str = '\t') -> str:
        """Return raw spectral file name from Sage db search file.

        Args:
            file_name (Path | str): Location of FragPipe db search psm file.
            delimiter (str, Optional): value delimiter used in tabular file.
                Defaults to '\t'.

        Returns:
            str: Name of raw spectral file
        """
        with open(file_name) as contents:
            header = contents.readline().replace('\n', '')
            if not all(x in header.split(delimiter) for x in cls.REQUIRED_FIELDS):
                raise TypeError("Invalid format supplied, expects Sage db search format.")
            # get first row
            line = contents.readline()
            line_cells = line.split(delimiter)

            # get source column
            source_idx = header.index["Spectrum File"]
            return Path(line_cells[source_idx]).stem


    def get_source_files(self) -> Sequence[str]:
        """Return all raw spectral file names from dataset, excluding file type
        suffix.

        Returns:
            Sequence[str]: All raw spectral file names in dataset.
        """
        source_file_col = self.data['Spectrum File']
        source_files: List[str] =  source_file_col\
            .dropna()\
            .drop_duplicates()\
            .apply(lambda x: Path(x).stem)\
            .tolist()

        return source_files

    
    def to_metapep_db_search(self, sample_name: str | None = None,
                             crap_dataset: pd.Series | None = None) -> MetaPepDbSearch:
        """Convert data table to MetaPepDbSearch format.
        
        Args:
            sample_name (str | None, optional): Name of sample to pass to 
                MetaPepDbSearch object, overwrites file name stored in self.
                Defaults to None.
            crap_dataset (pd.Series | None, optional): Set of crap peptides to 
                filter out of input peptide dataset. Defaults to None.

        Returns:
            MetaPepDbSearch: Db search data table in Metapep db search format.
        """
        df = self.data.copy(deep=True)

        # FragPipe reports the protein ID with highest evidence as determined by 
        # ProteinProphet. Other proteins are stored in 'Mapped Proteins'
        if self.ONLY_RAZOR_PROTEIN is True:
            df["Protein"] = (df["Protein"] + 
                            ", " + 
                            df["Mapped Proteins"].fillna("")).str.removesuffix(", ")
            df.drop('Mapped Proteins', axis=1)

        # drop unrelevant columns and rename columns to metapep format  
        df = self.data[self.REQUIRED_FIELDS]\
            .rename(columns={
                'Retention': 'RT',
                'Peptide Length': 'Length',
                'Hyperscore': 'Confidence', 
                'Intensity': 'Area',
                'Observed M/Z': 'm/z',
                'Observed Mass': 'Mass',
                'Protein': 'Accession',
                'Modified Peptide': 'PTM',
                'Spectrum File': 'Source File'
                })

        # filter out decoy or contamination sequences
        df = df[(df["Is Decoy"] == False) & (df["Is Contaminant"] == False)]
        
        # Wrangle sequence into consistent format (remove PTM, equalte L, I)
        df.loc[:, 'Sequence'] = df.loc[:, 'Peptide'].apply(wrangle_peptides)

        # remove file type suffix from Source File column
        df.loc[:, 'Source File'] = df.loc[:, 'Source File'].apply(
            lambda x: Path(x).stem
        )

        df.loc[:, "ppm"] = df[['Calculated M/Z', 'm/z']].apply(
            lambda x: mz_diff_to_ppm(x['Calculated M/Z'],
                                     x['m/z']),
            axis=1
        )
        df.drop('Calculated M/Z', axis=1)
        
        # extract scan number from 'scannr' column
        df['Scan'] = df['Spectrum'].str.extract(r"((?<=\.)\d+(?=\.))", expand=False)
        df.drop('Spectrum', axis=1)
        
        # filter out peptides from cRAP dataset
        if crap_dataset is not None:
            df = filter_crap(df, 'Sequence', crap_dataset)

        # replace accession delimiter to metapep format
        if self.ACCESSION_DELIMITER != MetaPepDbSearch.ACCESSION_DELIMITER:
            df.loc[:, 'Accession'] = df.loc[:, 'Accession'].str.replace(
                self.ACCESSION_DELIMITER, MetaPepDbSearch.ACCESSION_DELIMITER)
        
        # fetch all spectral file names present in dataset
        source_files: List[str] = self.get_source_files()
        
        # configure file name, if given the function argument, else the class file name
        if sample_name is None:
            sample_name = self.file_name
        
        # Put dataset into MetaPepDbSearch class object
        return MetaPepDbSearch(df,
                               self.DATA_FORMAT,
                               self.CONFIDENCE_FORMAT,
                               source_files,
                               sample_name)
