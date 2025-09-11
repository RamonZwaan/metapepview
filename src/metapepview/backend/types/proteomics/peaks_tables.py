import pandas as pd

from pathlib import Path
from typing import Type, Tuple, Sequence, List

from metapepview.backend.types.proteomics.proteomics_base_classes import DbSearchMethods, DeNovoMethods
from metapepview.backend.types.metapep_table import MetaPepDbSearch, MetaPepDeNovo
from metapepview.backend.utils import filter_crap, wrangle_peptides, peaks_10_drop_fraction_id


class PeaksDbSearchPsm11(DbSearchMethods):
    """Peaks database search psm table in Peaks Studio 11 format.
    Below is a summary of table columns within the format.
    
    Fields:
        "Peptide
        -10LgP
        Mass
        Length
        ppm
        m/z
        z
        RT
        Scan
        Area
        Feature Id
        Source File
        Accession
        PTM
        AScore
    """
    
    # expected columns within input data
    REQUIRED_FIELDS = ["Peptide", "-10LgP", "Mass", "Length", "ppm", "m/z", "z", "RT",
              "Scan", "Area", "Feature Id", "Source File", "Accession", "PTM",
              "AScore"]
    
    # columns expected to be numeric
    NUMERIC_FIELDS = ["-10LgP", "Mass", "Length", "ppm", "m/z",
                      "z", "RT", "Scan", "Area", "Feature Id"]
    
    ACCESSION_DELIMITER = ';'

    DATA_FORMAT = 'Peaks 11'
    CONFIDENCE_FORMAT = '-10lgp'
    
    
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
    def get_source_file(cls, file_name: Path | str) -> str:
        """Return raw spectral file name from peaks 11 db search path.

        Args:
            file_name (Path | str): Location of peaks 11 db search psm file.

        Returns:
            str: Name of raw spectral file
        """
        with open(file_name) as contents:
            header = contents.readline().replace('"', '').replace('\n', '')
            if not all(x in header.split(",") for x in cls.REQUIRED_FIELDS):
                raise TypeError("Invalid format supplied, expects Peaks Studio 11 db search format.")
            # get first row
            line = contents.readline()
            line_cells = line.split(",")

            # get source field and return without file type suffix
            return Path(line_cells[11]).stem


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
        df = self.data[["Peptide", "-10LgP", "Mass", "Length", "ppm", "m/z",
                        "z", "RT", "Scan", "Area", "Source File", "Accession",
                        "PTM"]]\
            .rename(columns={'-10LgP': 'Confidence', 'z': 'Charge'})
        
        # Wrangle sequence into consistent format (remove PTM, equalte L, I)
        df.loc[:, 'Sequence'] = df.loc[:, 'Peptide'].apply(wrangle_peptides)
        
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
                


class PeaksDbSearchPsm10(DbSearchMethods):
    """Peaks database search psm table in Peaks Studio 10 format.
    Below is a summary of table columns within the format.
    
    Fields:
        Peptide
        -10lgP
        Mass
        Length
        ppm
        m/z
        Z
        RT
        Area
        Fraction
        Id
        Scanfrom Chimera
        Source File
        Accession
        PTM
        AScore
        (Found By)      # Present in some files, but omitted here
    """
    
    # expected columns within input data
    REQUIRED_FIELDS = ["Peptide", "-10lgP", "Mass", "Length", "ppm", "m/z", "Z", "RT",
              "Area", "Fraction", "Id", "Scan", "from Chimera", "Source File",
              "Accession", "PTM", "AScore"]
    
    # columns expected to be numeric
    NUMERIC_FIELDS = ["-10lgP", "Mass", "Length", "ppm", "m/z", "Z", "RT",
                      "Area", "Fraction"]
    
    ACCESSION_DELIMITER = ':'
    
    DATA_FORMAT = 'Peaks 10'
    CONFIDENCE_FORMAT = '-10lgp'
    
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
    def get_source_file(cls, file_name: Path | str) -> str:
        """Return raw spectral file name from peaks 10 db search path.

        Args:
            file_name (Path | str): Location of peaks 10 db search psm file.

        Returns:
            str: Name of raw spectral file
        """
        with open(file_name) as contents:
            header = contents.readline().replace('\n', '')
            if not all(x in header.split(",") for x in cls.REQUIRED_FIELDS):
                raise TypeError("Invalid format supplied, expects Peaks Studio 10 db search format.")
            # get first row
            line = contents.readline()
            line_cells = line.split(",")

            # get source column and match against dict
            return Path(line_cells[13]).stem
     
        
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
        
    
    def to_peaks_db_psm_11(self) -> PeaksDbSearchPsm11:
        """Convert data table to Peaks studio 11 db search psm format.

        Returns:
            MetaPepDeNovo: Db search psm data table in Peaks studio 11 format.
        """
        # dict to map peaks 10 col names (keys) with peaks 11 names (values)
        col_name_changes = {"-10lgP": "-10LgP",
                            "Z": "z",
                            "Id": "Feature Id"}
        # columns to remove from data
        drop_cols = ["Fraction", "from Chimera"]
        # columns with fraction information added
        frac_id_cols = ["Scan"]
        
        # apply changes of peaks 10 dataset
        peaks_11_dataset = self.data.copy()
        peaks_11_dataset[frac_id_cols] = peaks_11_dataset[frac_id_cols].map(peaks_10_drop_fraction_id)
        peaks_11_dataset.rename(columns=col_name_changes, inplace=True)
        peaks_11_dataset.drop(drop_cols, axis=1, inplace=True)
        peaks_11_dataset["Accession"].apply(lambda x: x.replace(":", ";") if isinstance(x, str) else x)
        
        return PeaksDbSearchPsm11(peaks_11_dataset)
    
    
    def to_metapep_db_search(self,
                             sample_name: str | None = None,
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
        df = self.data[["Peptide", "-10lgP", "Mass", "Length", "ppm", "m/z",
                        "Z", "RT", "Area", "Scan", "Source File", "Accession",
                        "PTM"]]\
            .rename(columns={'-10lgP': 'Confidence', 'Z': 'Charge'})
        
        # Wrangle sequence into consistent format (remove PTM, equalte L, I)
        df.loc[:, 'Sequence'] = df.loc[:, 'Peptide'].apply(wrangle_peptides)
        
        # scan column includes fraction if multiple present, need to be removed
        df.loc[:, 'Scan'] = fraction_col_to_numeric(df, 'Scan')
        
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
        
        # fetch all spectral file names present in dataset
        source_files: List[str] = df.loc[:, 'Source File'].unique().tolist()
        
        # configure file name, if given the function argument, else the class file name
        if sample_name is None:
            sample_name = self.file_name
        
        return MetaPepDbSearch(df,
                               self.DATA_FORMAT,
                               self.CONFIDENCE_FORMAT,
                               source_files,
                               sample_name)
     
        
class PeaksDeNovo11(DeNovoMethods):
    """Peaks de novo table format in Peaks Studio 11 format.
    Below is a summary of expected fields for this format
    
    Fields:
        Source File
        Scan
        Peptide
        Tag length
        ALC (%)
        Length
        m/z
        z
        RT
        Area
        Mass
        ppm
        PTM
        local confidence (%)
        mode
        tag(>=0.0%)
        Feature Id
    """
    
    # expected columns within input data
    REQUIRED_FIELDS = ["Source File", "Scan", "Peptide", "Tag length", "ALC (%)",
              "Length", "m/z", "z", "RT", "Area", "Mass", "ppm", "PTM",
              "local confidence (%)", "mode", "tag(>=0.0%)", "Feature Id"]
    
    # columns expected to be numeric
    NUMERIC_FIELDS = ["Scan", "Tag length", "ALC (%)", "Length", "m/z", "z",
                      "RT", "Area", "Mass", "ppm", "Feature Id"]
    
    DATA_FORMAT = 'Peaks 11'
    CONFIDENCE_FORMAT = 'ALC'

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
    def get_source_file(cls, file_name: Path | str) -> str:
        """Return raw spectral file name from peaks 11 de novo path.

        Args:
            file_name (Path | str): Location of peaks 11 de novo file.

        Returns:
            str: Name of raw spectral file
        """
        with open(file_name) as contents:
            header = contents.readline().replace('"', '').replace('\n', '')
            if not all(x in header.split(",") for x in cls.REQUIRED_FIELDS):
                raise TypeError("Invalid format supplied, expects Peaks Studio 10 db search format.")
            # get first row
            line = contents.readline()
            line_cells = line.split(",")

            # get source column and match against dict
            return Path(line_cells[0]).stem
        

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
        df = self.data[["Source File", "Scan", "Peptide", "ALC (%)", "Length",
                        "m/z", "z", "RT", "Area", "Mass", "ppm", "PTM"]]\
            .rename(columns={'ALC (%)': 'Confidence', 'z': 'Charge'})
        
        # Wrangle sequence into consistent format (remove PTM, equalte L, I)
        df.loc[:, 'Sequence'] = df.loc[:, 'Peptide'].apply(wrangle_peptides)
        
        # remove file type suffix from Source File column
        df.loc[:, 'Source File'] = df.loc[:, 'Source File'].apply(
            lambda x: Path(x).stem
        )
        

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
                

class PeaksDeNovo10(DeNovoMethods):
    """Peaks de novo table format in Peaks Studio 10 format.
    Below is a summary of expected fields for this format
    
    Fields:
        Fraction
        Source File
        Feature
        Peptide
        Scan
        Tag Length
        Denovo Score
        ALC (%)
        length
        m/z
        z
        RT
        Predict RT
        Area
        Mass
        ppm
        PTM
        local confidence (%)
        tag (>=0%)
        mode
    """
    
    # expected columns within input data
    REQUIRED_FIELDS = ["Fraction", "Source File", "Feature", "Peptide", "Scan",
              "Tag Length", "ALC (%)", "length", "m/z", "z", "RT", "Area",
              "Mass", "ppm", "PTM", "local confidence (%)", "tag (>=0%)", "mode"]
    
    # columns expected to be numeric
    NUMERIC_FIELDS = ["Fraction", "Tag Length","ALC (%)", "length", "m/z", "z",
                      "RT", "Area", "Mass", "ppm"]

    DATA_FORMAT = 'Peaks 10'
    CONFIDENCE_FORMAT = "ALC"
    
    def __init__(self,
                 data: pd.DataFrame,
                 file_name: str | None):
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
    def get_source_file(cls, file_name: Path | str) -> str:
        """Return raw spectral file name from peaks 10 de novo path.

        Args:
            file_name (Path | str): Location of peaks 10 de novo file.

        Returns:
            str: Name of raw spectral file
        """
        with open(file_name) as contents:
            header = contents.readline().replace('"', '').replace('\n', '')
            if not all(x in header.split(",") for x in cls.REQUIRED_FIELDS):
                raise TypeError("Invalid format supplied, expects Peaks Studio 10 de novo format.")
            # get first row
            line = contents.readline()
            line_cells = line.split(",")

            # get source column and match against dict
            return Path(line_cells[1]).stem


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
    
        
    def to_peaks_de_novo_11(self) -> PeaksDeNovo11:
        """Convert data table to Peaks studio 11 de novo format.

        Returns:
            PeaksDeNovo11: De novo data table in peaks studio 11 format
        """
        # dict to map peaks 10 col names (keys) with peaks 11 names (values)
        col_name_changes = {"Tag Length": "Tag length",
                            "length": "Length",
                            "Feature": "Feature Id",
                            "tag (>=0%)": "tag(>=0.0%)"}
        
        drop_cols = ("Fraction")

        frac_id_cols = ["Feature", "Scan"]
        
        # apply header changes
        peaks_11_dataset = self.data.copy()
        peaks_11_dataset[frac_id_cols] = peaks_11_dataset[frac_id_cols].map(peaks_10_drop_fraction_id)
        peaks_11_dataset.rename(columns=col_name_changes, inplace=True)
        peaks_11_dataset.drop(drop_cols, axis=1, inplace=True)
        
        # accession only given in denovo only data
        if "Accession" in peaks_11_dataset.columns:
            peaks_11_dataset["Accession"].apply(lambda x: x.replace(":", ";") if isinstance(x, str) else x)
        
        return PeaksDeNovo11(peaks_11_dataset)
    
    
    def to_metapep_de_novo(self,
                           sample_name: str | None = None,
                           crap_dataset: pd.Series | None = None) -> MetaPepDeNovo:
        """Convert data table to MetaPepDeNovo format.

        Returns:
            MetaPepDeNovo: De novo data table in Metapep de novo format.
        """
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
        df = self.data[["Source File", "Peptide", "Scan", "ALC (%)", "length",
                        "m/z", "z", "RT", "Area", "Mass", "ppm", "PTM"]]\
            .rename(columns={'ALC (%)': 'Confidence', 'z': 'Charge', 'length': 'Length'})
        
        # for multiple fractions, 'Scan' column contains fraction identity,
        # this needs to be removed in order to match spectral files to peptide
        df.loc[:, 'Scan'] = fraction_col_to_numeric(df, 'Scan')
        
        # Wrangle sequence into consistent format (remove PTM, equalte L, I)
        df.loc[:, 'Sequence'] = df.loc[:, 'Peptide'].apply(wrangle_peptides)
        
        # remove file type suffix from Source File column
        df.loc[:, 'Source File'] = df.loc[:, 'Source File'].apply(
            lambda x: Path(x).stem
        )
        
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


################################################################################
# utils
################################################################################

def fraction_col_to_numeric(df: pd.DataFrame, col: str) -> pd.Series:
    """Convert string column to numerical by removing prefix data from numbers.
    Some numeric columns within Peaks (like 'Scan') get altered with prefix
    data in some configurations. This hinders matching values down the processing
    pipeline. 
    

    Args:
        df (pd.DataFrame): Input dataset.
        col (str): Column to convert.

    Raises:
        ValueError: Conversion failed.
    """
    # scan column includes fraction if multiple present, need to be removed
    scan_val = df.loc[df.index[0], 'Scan']
    if type(scan_val) == str and len(scan_val.split(':')) > 1:
        col_series = df.loc[:, 'Scan'].apply(lambda x: x.split(':')[-1] if not isinstance(x, float) else x)
    else:
        col_series = df.loc[:, 'Scan']

    try:
        col_series = pd.to_numeric(col_series)
    except Exception as e:
        raise ValueError("Scan column conversion failed")

    return col_series