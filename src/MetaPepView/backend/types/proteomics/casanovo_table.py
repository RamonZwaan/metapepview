import io
import re
from pathlib import Path
import pandas as pd
import numpy as np
from typing import List, IO, Self, Dict

from metapepview.backend.types.proteomics.proteomics_base_classes import DeNovoMethods
from metapepview.backend.types.metapep_table import MetaPepDeNovo
from metapepview.backend.utils import memory_to_stringio,\
    wrangle_peptides,\
    spectrum_id_to_scan_number,\
    mz_diff_to_ppm,\
    filter_crap

from metapepview.constants import GlobalConstants, PhysicalConstants

class CasanovoDeNovo(DeNovoMethods):

    REQUIRED_FIELDS = [
        "sequence",
        "PSM_ID",
        "accession",
        "unique",
        "database",
        "database_version",
        "search_engine",
        "search_engine_score[1]",
        "modifications",
        "retention_time",
        "charge",
        "exp_mass_to_charge",
        "calc_mass_to_charge",
        "spectra_ref",
        "opt_ms_run[1]_aa_scores"]
    NUMERIC_FIELDS = ["search_engine_score[1]", 
                      "retention_time", 
                      "charge", 
                      "exp_mass_to_charge",
                      "calc_mass_to_charge"]

    DATA_FORMAT = 'Casanovo'
    # Confidence score between [-1, 1], 
    # where score <0 means mismatch between predicted and measured mass
    CONFIDENCE_FORMAT = 'Score'

    # custom regex patterns to parse dataset
    MS_RUN_PATTERN = re.compile(r"ms_run\[[0-9+]\]")
    PTM_PATTERN = re.compile(r"\[.*\]")

    def __init__(self,
                 data: pd.DataFrame,
                 file_name: str | None = None,
                 source_files: Dict[str, str] | None = None):
        success, msg = self.validate_input(data)

        if success is False:
            raise ValueError(msg)
        
        data = data[self.REQUIRED_FIELDS]
        
        self._data = data
        self._file_name = file_name

        # mapping run name as present in Casanovo dataset to source file
        self._source_files = source_files

    @property
    def data(self) -> pd.DataFrame:
       return self._data
    
    @property
    def file_name(self) -> str | None:
       return self._file_name
    
    @classmethod
    def get_source_file(cls, file_name: Path | str) -> str:
        """Return raw spectral file name from Casanovo mztab path.

        Note:
            Casanovo stores the name of the input file as source, this is a
            'mzML' or 'MGF' format spectral file. This may differ from the name
            of the vendor specific raw spectral file.

        Args:
            file_name (Path | str): Location of Casanovo mztab file.

        Returns:
            str: Name of raw spectral file
        """
        failed_msg = "Invalid format supplied, cannot extract raw file name"

        with open(file_name) as contents:
            for line in contents.readlines():
                if line.startswith("MTD\tms_run[1]"):
                    line_fields = line.split("\t")
                    if len(line_fields) > 3:
                        raise TypeError(failed_msg)
                    # return file name without suffix
                    else:
                        return Path(line_fields[2]).stem
            raise TypeError(failed_msg)


    def get_source_files(self):
        """Return all raw spectral file names from dataset, excluding file type
        suffix.

        Note:
            Casanovo stores the name of the input file as source, this is a
            'mzML' or 'MGF' format spectral file. This may differ from the name
            of the vendor specific raw spectral file.

        Returns:
            Sequence[str]: All raw spectral file names in dataset.
        """
        return self._source_files
       
    @classmethod
    def __fetch_all_source_files(cls, contents: IO[str]) -> Dict[str, str]:
        failed_msg = "Invalid format supplied, cannot extract raw file names"

        source_files = dict()

        # get first row
        for line in contents.readlines():
            if line.startswith("MTD\tms_run"):
                line_fields = line.split("\t")
                if len(line_fields) > 3:
                    continue
                # extract internal name of MS run without suffix
                run_name = re.search(cls.MS_RUN_PATTERN, line_fields[1])[0]
                if run_name is None:
                    continue

                # extract file name
                file_name = Path(line_fields[2]).stem

                source_files[run_name] = file_name

            elif not line.startswith("MTD"):
                break
        
        # only return if at least one source file was found
        if len(source_files) > 0:
            return source_files
        else:
            raise TypeError(failed_msg)


    @classmethod
    def read_file(cls, path: str | Path) -> Self:
        """Read input file and return data as instance of class object.
        Overwrite method with function that adheres to input format.

        Args:
            path (str | Path): Location of file

        Returns:
            T: Instance of class object
        """
        file_name = Path(path).stem

        # open file to parse metadata fields
        with open(path) as file_contents:
            cls_obj = cls.read_file_buffer(file_contents, file_name)

        return cls_obj


    @classmethod
    def read_file_buffer(cls,
                         file_buffer: str | IO[str],
                         file_name: str | None = None) -> Self:
        """Read input file buffer and return data instance of class object.
        Use this when data is imported into memory as string buffer. Overwrite
        method with function that adheres to input format.
        
        >>> CasanovoDeNovo.read_file("example.mztab")
        Is equal to
        >>> CasanovoDeNovo.read_file_buffer(open("example.mztab"))

        Args:
            file_buffer (str | IO[str]): input file data or buffer
            file_name (str | None, optional): Name of input file. 
                Defaults to None.

        Returns:
            T: Instance of class object
        """
        if isinstance(file_buffer, str):
            file_buffer = memory_to_stringio(file_buffer)

        try:
            source_files = cls.__fetch_all_source_files(file_buffer)
        except TypeError:
            return TypeError("failed to parse Casanovo metadata: No source file found.")
        
        # go back to start of buffer
        file_buffer.seek(0) 

        # get to header row and count number of metadata rows
        for c, line in enumerate(file_buffer.readlines()):
            if line.startswith("PSH"):
                metadata_rows = c
                break
            elif not line.startswith("MTD"):
                raise TypeError("Invalid mztab format: header row not found after metadata.")
        file_buffer.seek(0) 
        
        df = pd.read_csv(file_buffer, 
                         delimiter="\t",
                         skiprows=metadata_rows,
                         low_memory=False).iloc[:, 1:]
    
        cls_obj = cls(df, file_name, source_files) # type: ignore
        return cls_obj


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
        df = self.data[[
            "sequence", 
            "search_engine_score[1]", 
            "modifications", 
            "retention_time",
            "charge",
            "exp_mass_to_charge",
            "calc_mass_to_charge",
            "spectra_ref"]]\
            .rename(columns={
                'sequence': 'Peptide',
                'retention_time': 'RT',
                'calc_mass_to_charge': 'm/z',
                'search_engine_score[1]': 'Confidence',
                'charge': 'Charge',
                'modifications': 'PTM'
                })
        
        # Wrangle sequence into consistent format: equalte L, I.
        # Apply custom PTM filter to remove text inside squate brackets
        df.loc[:, 'Sequence'] = df.loc[:, 'Peptide']\
            .apply(lambda x: re.sub(self.PTM_PATTERN, "", x))\
            .apply(wrangle_peptides)
        

        # use spectral id column to get source file name for each row.
        def spectra_ref_to_run(spectra_ref: str) -> str | None:
            run_match = re.search(self.MS_RUN_PATTERN, spectra_ref)
            if run_match is not None:
                return run_match[0]
            return None

        df.loc[:, 'Source File'] = df.loc[:, 'spectra_ref'].apply(
            lambda x: self._source_files.get(spectra_ref_to_run(x),
                                             np.nan)
        )

        # attempt to fetch scan number from spectra_ref column
        df.loc[:, 'Scan'] = df['spectra_ref'].apply(
            spectrum_id_to_scan_number
        )

        df.loc[:, 'ppm'] = df[['exp_mass_to_charge', 'm/z']].apply(
            lambda x: mz_diff_to_ppm(x['exp_mass_to_charge'],
                                     x['m/z']),
            axis=1
        )

        df.loc[:, "Length"] = df['Sequence'].apply(
            lambda x: len(x) if isinstance(x, str) else np.nan
        )

        # Signal intensities are not provided by Casanovo
        df.loc[:, "Area"] = np.nan

        h_mass = PhysicalConstants.proton_mass

        # calculate theoretical mass of peptides as [M]
        df.loc[:, "Mass"] = df['exp_mass_to_charge'] * df['Charge'] - (df['Charge'] * h_mass)

        # drop columns not part of MetaPepDeNovo dataset
        df = df.drop(['exp_mass_to_charge', 'spectra_ref'], axis=1)
        
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
    