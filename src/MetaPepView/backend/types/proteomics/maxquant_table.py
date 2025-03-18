from pathlib import Path
import pandas as pd
from typing import Type, List, Sequence, Self, IO

from .proteomics_base_classes import DbSearchMethods
from ..metapep_table import MetaPepDbSearch
from backend.utils import filter_crap, wrangle_peptides


class MaxQuantDbSearch(DbSearchMethods):
    """MaxQuant evidence datatable from 'evidence.txt'.
    Below is a summary of table columns within the format
    
    Note: Table format depends on selected PTM's during database search, as
    columns are generated for each PTM. This is denoted below with {X}, where X
    is the PTM. An example for X is 'Deamidation (NQ)'. Currently, these columns
    are kept in the dataset but not considered in processing.
    
    Fields:
        Sequence
        Length
        Modifications
        Modified sequence
        {X Probabilities}   # any PTM selected in maxquant gets its own column
        {X Score Diffs}     # any PTM selected in maxquant gets its own column
        Acetyl (Protein N-term)
        {X}                 # any PTM selected in maxquant gets its own column
        Missed cleavages
        Proteins
        Leading proteins
        Leading razor protein
        Type
        Raw file
        Experiment
        MS/MS m/z
        Charge
        m/z
        Mass
        Uncalibrated - Calibrated m/z [ppm]
        Uncalibrated - Calibrated m/z [Da]
        Mass error [ppm]
        Mass error [Da]
        Uncalibrated mass error [ppm]
        Uncalibrated mass error [Da]
        Max intensity m/z 0
        Retention time
        Retention length
        Calibrated retention time
        Calibrated retention time start
        Calibrated retention time finish
        Retention time calibration
        Match time difference
        Match m/z difference
        Match q-value
        Match score
        Number of data points
        Number of scans
        Number of isotopic peaks
        PIF
        Fraction of total spectrum
        Base peak fraction
        PEP
        MS/MS count
        MS/MS scan number
        MS/MS scan numbers
        MS3 scan numbers
        Score
        Delta score
        Combinatorics
        Intensity
        Reverse
        Potential contaminant
        id
        Protein group IDs
        Peptide ID
        Mod. peptide ID
        MS/MS IDs
        Best MS/MS
        Deamidation (NQ) site IDs
        Oxidation (M) site IDs
        Taxonomy IDs
        Taxonomy names
        Mass deficit

    """
    
    # expected columns within input data
    REQUIRED_FIELDS = ['Sequence', 'Length', 'Modifications', 'Modified sequence', 
              'Acetyl (Protein N-term)', 'Missed cleavages', 'Proteins',
              'Leading proteins', 'Leading razor protein', 'Type', 'Raw file',
              'Experiment', 'MS/MS m/z', 'Charge', 'm/z', 'Mass',
              'Uncalibrated - Calibrated m/z [ppm]',
              'Uncalibrated - Calibrated m/z [Da]', 'Mass error [ppm]',
              'Mass error [Da]', 'Uncalibrated mass error [ppm]',
              'Uncalibrated mass error [Da]', 'Max intensity m/z 0',
              'Retention time', 'Retention length', 'Calibrated retention time',
              'Calibrated retention time start', 'Calibrated retention time finish',
              'Retention time calibration', 'Match time difference',
              'Match m/z difference', 'Match q-value', 'Match score',
              'Number of data points', 'Number of scans',
              'Number of isotopic peaks', 'PIF', 'Fraction of total spectrum',
              'Base peak fraction', 'PEP', 'MS/MS count', 'MS/MS scan number',
              'MS/MS scan numbers', 'MS3 scan numbers', 'Score', 'Delta score',
              'Combinatorics', 'Intensity', 'Reverse', 'Potential contaminant',
              'id', 'Protein group IDs', 'Peptide ID', 'Mod. peptide ID',
              'MS/MS IDs', 'Best MS/MS', 'Deamidation (NQ) site IDs',
              'Oxidation (M) site IDs', 'Taxonomy IDs', 'Taxonomy names',
              'Mass deficit']

    # columns expected to be numeric. Note, only relevant columns considered
    NUMERIC_FIELDS = ['Length', 'Missed cleavages', 'Charge', 'm/z', 'Mass',
                      'Retention time', 'MS/MS count', 'MS/MS scan number',
                      'Score', 'Delta score', 'Intensity']
    
    ACCESSION_DELIMITER = ';'
    
    DATA_FORMAT = 'MaxQuant'
    CONFIDENCE_FORMAT = '-10lgp'
    
    
    def __init__(self,
                 data: pd.DataFrame,
                 file_name: str | None = None):
        success, msg = MaxQuantDbSearch.validate_input(data)

        if success is False:
            raise ValueError(msg)
        
        # keep all fields, including specific PTM fields
        # data = data[self.REQUIRED_FIELDS]
        
        self._data = data
        self._file_name = file_name
   
    @property
    def data(self) -> pd.DataFrame:
       return self._data

    @property
    def file_name(self) -> str | None:
       return self._file_name

    @classmethod
    def get_source_file(cls, file_name: Path | str) -> str:
        """Return raw spectral file name from peaks 10 db search path.

        Args:
            file_name (Path | str): Location of peaks 10 db search psm file.

        Returns:
            str: Name of raw spectral file
        """
        with open(file_name) as contents:
            header = contents.readline().replace('"', '')
            if header.split(",") != cls.REQUIRED_FIELDS:
                raise TypeError("Invalid format supplied, expects Peaks Studio 10 db search format.")
            # get first row
            line = contents.readline()
            line_cells = line.split(",")

            # get source column and match against dict
            return Path(line_cells[10]).stem


    # overwrite read_file method to change delimiter to tab
    @classmethod
    def read_file(cls: Type[Self], path: str | Path) -> Self:
        """Read input file and return data as instance of class object.

        Args:
            cls (Type[MaxQuantDbSearch]): Class type
            path (str | Path): Location of file

        Returns:
            MaxQuantDbSearch: Instance of class object
        """
        return cls._read_csv(path, '\t')
    
    
    # overwrite read_file method to change delimiter to tab
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
            cls (Type[MaxQuantDbSearch]): Class type
            file_buffer (str | IO[str]): input file data or buffer
            file_name (str | None, optional): Name of input file. 
                Defaults to None.

        Returns:
            MaxQuantDbSearch: Instance of class object
        """
        return super()._read_csv_buffer(file_buffer,
                                        file_name, 
                                        '\t')

    def get_source_files(self) -> Sequence[str]:
        """Return all raw spectral file names from dataset, excluding file type
        suffix.

        Returns:
            Sequence[str]: All raw spectral file names in dataset.
        """
        source_file_col = self.data['Raw file']
        source_files: List[str] =  source_file_col\
            .dropna()\
            .unique()\
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
        # drop unrelevant columns and rename columns to metapep format
        df = self.data[['Modified sequence', 'Sequence', 'Retention time',
                        'MS/MS scan number', 'm/z', 'Charge', 'Mass error [ppm]',
                        'Length', 'Score', 'Intensity', 'Mass', 'Proteins',
                        'Modifications', 'Raw file']]\
            .rename(columns={
                'Modified sequence': 'Peptide',
                'Retention time': 'RT',
                'MS/MS scan number': 'Scan',
                'Mass error [ppm]': 'ppm',
                'Score': 'Confidence',
                'Intensity': 'Area',
                'Proteins': 'Accession',
                'Modifications': 'PTM',
                'Raw file': 'Source File'})
        
        # Wrangle Sequence into consistent format (equalte L, I)
        df.loc[:, 'Sequence'] = df.loc[:, 'Sequence'].apply(wrangle_peptides)
        
        # remove file type suffix from Source File column
        df.loc[:, 'Source File'] = df.loc[:, 'Source File'].apply(
            lambda x: Path(x).stem
        )

        # filter out peptides from cRAP dataset
        if crap_dataset is not None:
            df = filter_crap(df, 'Sequence', crap_dataset)

        # replace accession delimiter to metapep format
        if self.ACCESSION_DELIMITER != MetaPepDbSearch.ACCESSION_DELIMITER:
            df.loc[:, 'Accession'] = df.loc[:, 'Accession'].str.replace(
                self.ACCESSION_DELIMITER, MetaPepDbSearch.ACCESSION_DELIMITER)
        
        # fetch all spectral file names present in dataset, ignore nan
        source_files: List[str] = df.loc[:, 'Source File'].dropna().unique().tolist()
        
        # configure file name, if given the function argument, else the class file name
        if sample_name is None:
            sample_name = self.file_name
        
        # Put dataset into MetaPepDbSearch class object
        return MetaPepDbSearch(df,
                               self.DATA_FORMAT,
                               self.CONFIDENCE_FORMAT,
                               source_files,
                               sample_name)
        