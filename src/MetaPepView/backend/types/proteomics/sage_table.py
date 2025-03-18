import pandas as pd
import numpy as np

from pathlib import Path
from typing import List, Sequence, Type, Self, IO

from .proteomics_base_classes import DbSearchMethods
from constants import *
from ..metapep_table import MetaPepDbSearch
from backend.utils import filter_crap, wrangle_peptides


class SageDbSearch(DbSearchMethods):
    """Sage database search psm table.
    Below is a summary of table columns within the format.
    Hyperscore is taken as confidence parameter.
    
    Fields:
        peptide
        proteins (; delimited)
        num_proteins
        filename
        scannr
        rank
        label
        expmass
        calcmass
        charge
        peptide_len
        missed_cleavages
        isotope_error
        precursor_ppm
        fragment_ppm
        hyperscore
        delta_next
        delta_best
        rt
        aligned_rt
        predicted_rt
        delta_rt_model
        matched_peaks
        longest_b
        longest_y
        longest_y_pct
        matched_intensity_pct
        scored_candidates
        poisson
        sage_discriminant_score
        posterior_error
        spectrum_q
        peptide_q
        protein_q
        ms2_intensity
    """
    
    # expected columns within input data
    REQUIRED_FIELDS = ["peptide", "hyperscore", "expmass", "peptide_len",
                       "rank", "label", "precursor_ppm", "charge", "rt", 
                       "scannr", "filename", "proteins"]
    
    # columns expected to be numeric
    NUMERIC_FIELDS = ["hyperscore", "expmass", "peptide_len", "precursor_ppm", 
                      "charge", "rt"]
    
    ACCESSION_DELIMITER = ';'

    DATA_FORMAT = 'Sage'
    CONFIDENCE_FORMAT = 'Hyperscore'
    
    
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
            file_name (Path | str): Location of peaks 11 db search psm file.
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
            return Path(line_cells[3]).stem


    def get_source_files(self) -> Sequence[str]:
        """Return all raw spectral file names from dataset, excluding file type
        suffix.

        Returns:
            Sequence[str]: All raw spectral file names in dataset.
        """
        source_file_col = self.data['filename']
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
        df = self.data[["peptide", "hyperscore", "expmass", "peptide_len", 
                        "precursor_ppm", "charge", "rt", "scannr", "rank", 
                        "label", "filename", "proteins"]]\
            .rename(columns={
                'peptide': 'Peptide',
                'rt': 'RT',
                'charge': 'Charge',
                'precursor_ppm': 'ppm',
                'peptide_len': 'Length',
                'hyperscore': 'Confidence', 
                'expmass': 'Mass',
                'proteins': 'Accession',
                'filename': 'Source File'
                })
        # Keep only best match per scan and filter out decoy sequences
        df = df[df['rank'] == 1]
        df = df[df['label'] == 1]
        df.drop(['rank', 'label'], axis=1)
        
        # Wrangle sequence into consistent format (remove PTM, equalte L, I)
        df.loc[:, 'Sequence'] = df.loc[:, 'Peptide'].apply(wrangle_peptides)

        # remove file type suffix from Source File column
        df.loc[:, 'Source File'] = df.loc[:, 'Source File'].apply(
            lambda x: Path(x).stem
        )
        
        # extract scan number from 'scannr' column
        df['Scan'] = df['scannr'].str.extract(r"((?<=scan=)\d+$)")
        df.drop('scannr', axis=1)
        
        # compute m/z from measured mass + protons divided by charge
        h_mass = PhysicalConstants.proton_mass
        df['m/z'] = (df['Mass'] + df['Charge'] * h_mass) / df['Charge']
        
        # TODO: Extract PTM data from peptide column
        df['PTM'] = np.nan

        # precursor peak intensity is not supplied in Sage output
        df['Area'] = np.nan
        
        # filter out peptides from cRAP dataset
        if crap_dataset is not None:
            df = filter_crap(df, 'Sequence', crap_dataset)
        
        # replace accession delimiter to metapep format
        if self.ACCESSION_DELIMITER != MetaPepDbSearch.ACCESSION_DELIMITER:
            df.loc[:, 'Accession'] = df.loc[:, 'Accession'].str.replace(
                self.ACCESSION_DELIMITER, MetaPepDbSearch.ACCESSION_DELIMITER)
        
        # fetch all spectral file names present in dataset
        source_files: List[str] = df.loc[:, 'Source File'].unique().tolist()
        
        # configure file name, if given the function argument, else the class file name
        if sample_name is None:
            sample_name = self.file_name
        
        # Put dataset into MetaPepDbSearch class object
        return MetaPepDbSearch(df,
                               self.DATA_FORMAT,
                               self.CONFIDENCE_FORMAT,
                               source_files,
                               sample_name)
