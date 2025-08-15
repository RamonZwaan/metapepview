"""
This module contains custom Metapepview data formats for processing within the
dashboard. All proteomics input formats are converted into Metapepview format
prior to taxonomic/functional/QC processing to ensure common behavior between
any input type.
"""

import io
import pandas as pd
import numpy as np
from typing import Tuple, Sequence, Set, Self, TypeVar, Dict, Callable, Any, List
from pathlib import Path
from itertools import chain
import json

from .base_classes import DataValidator
from backend.utils import mode_func, to_json, decompress_string, convert_deprecated_metapeptable_naming
from constants import GlobalConstants
from .definitions import *




class MetaPepTable(DataValidator):
    """Table of Annotated peptides used for taxonomy and functional analysis
    within the MetaPepView interface. Below is a summary of the table format
    expected within the tool.

    Peptide info:
        Peptide
        Sequence
    
    DB search:
        PSM Count
        RT
        Scan
        m/z
        Charge
        ppm
        Length
        Feature Id
        Confidence
        Area
        Mass
        Accession
        PTM
        Source File
    
    De novo (Optional):
        De Novo Confidence
        De Novo Area
        De Novo Match Count
        De Novo Scan
        De Novo Source File
    
    Taxonomy annotation 
    (lineage = [superkingdom, phylum, class, order, family, genus, species]):
        Taxonomy Id
        Taxonomy Name
        {lineage} Id
        {lineage} Name
        
    Global taxonomy annotation (Optional):
        Global Taxonomy Id
        Global Taxonomy Name
        {lineage} Id (global search)
        {lineage} Name (global search)

    Functional annotation:
        KEGG KO
        CAZy
        GO
        
    Metapep Metadata fields:
        Taxonomy Annotation
        Taxonomy DB Name
        Functional Annotation
        Functional Annotation DB Name
        DB Search Imported
        De Novo Imported
        Global Taxonomy Annotation
        Sample Name
    
    Global Metapep metadata:
        Taxonomy DB Format
        Functional Annotation DB Format
        DB Search Format
        DB Search Confidence Format
        De Novo Format
        De Novo Confidence Format
        

    """
    
    ACCESSION_DELIMITER: str = ';'
    
    # Describe different column groups within the complete metapep table
    DB_SEARCH_FIELDS = ['Peptide', 'Sequence', 'PSM Count', 'RT', 'Scan', 'm/z',
                        'Charge', 'ppm', 'Length', 'Confidence', 'Area', 'Mass', 
                        'Accession', 'PTM', 'Source File']
    DE_NOVO_FIELDS = ['De Novo Confidence', 'De Novo Area',
                      'De Novo Match Count', 'De Novo Scan',
                      'De Novo Source File']
    LINEAGE_IDS = [i + ' Id' for i in GlobalConstants.standard_lineage_ranks]
    LINEAGE_NAMES = [i + ' Name' for i in GlobalConstants.standard_lineage_ranks]
    ANNOTATION_FIELDS = ['Taxonomy Id', 'Taxonomy Name'] + LINEAGE_IDS + LINEAGE_NAMES

    # annotation of peptide sequences against global reference database (OPTIONAL)
    GLOBAL_LINEAGE_IDS = [i + f' Id{GlobalConstants.global_annot_suffix}' for i in GlobalConstants.standard_lineage_ranks]
    GLOBAL_LINEAGE_NAMES = [i + f' Name{GlobalConstants.global_annot_suffix}' for i in GlobalConstants.standard_lineage_ranks]
    GLOBAL_ANNOTATION_FIELDS = ['Global Taxonomy Id', 'Global Taxonomy Name'] + GLOBAL_LINEAGE_IDS + GLOBAL_LINEAGE_NAMES

    METADATA_FIELDS = ['Taxonomy DB Name', 
                       'Functional Annotation DB Name',
                       'Global Taxonomy Annotation', 
                       'De Novo Imported',
                       'Sample Name']
    
    OBJ_ARGS = ['Taxonomy DB Format', 
                'Functional DB Format', 
                'DB Search Format', 
                'DB Search Confidence Format', 
                'De Novo Format', 
                'De Novo Confidence Format', 
                'Experiment Name']
    
    # TODO: Describe fields required for functional annotation of peptides
    # TODO: Manage both mandatory fields and optional fields
    REQUIRED_FIELDS = DB_SEARCH_FIELDS + DE_NOVO_FIELDS + ANNOTATION_FIELDS + METADATA_FIELDS
    NUMERIC_FIELDS = ['PSM Count', 'RT', 'Scan', 'm/z', 'Charge', 'ppm',
                      'Length', 'Confidence', 'Area', 'Mass',
                      'De Novo Confidence', 'De Novo Area',
                      'De Novo Match Count', 'De Novo Scan']


    def __init__(self,
                 data: pd.DataFrame,
                 taxonomy_db_format: TaxonomyFormat | None,
                 functional_db_format: FuncAnnotFormat | None,
                 db_search_format: DbSearchSource | None,
                 db_search_confidence_format: DbSearchConfFormat | None,
                 de_novo_format: DeNovoSource | None,
                 de_novo_confidence_format: DeNovoConfFormat | None,
                 experiment_name: str | None = None):
        success, msg = MetaPepTable.validate_input(data)

        if success is False:
            raise ValueError(msg)

        # extra columns in input data are not caught in validation, but dropped
        # data = data[self.REQUIRED_FIELDS]
        
        # add data
        self._data = data
        
        # global metadata, these have to be consistent across complete dataset
        self._taxonomy_db_format = taxonomy_db_format
        self._functional_db_format = functional_db_format
        self._db_search_format = db_search_format
        self._db_search_confidence_format = db_search_confidence_format
        self._de_novo_format = de_novo_format
        self._de_novo_confidence_format = de_novo_confidence_format
        self._experiment_name = experiment_name
        
        # grab source files. Used for comparing metapep tables
        self._source_files = self._data['Source File'].dropna().unique().tolist()
        self._sample_names = self._data['Sample Name'].dropna().unique().tolist()
    
    @property
    def data(self) -> pd.DataFrame:
        return self._data
    
    @data.setter
    def data(self, new_data: pd.DataFrame):
        # test that converted data follows input format
        valid, msg = self.validate_input(new_data)
        
        if not valid:
            raise ValueError(msg)
        
        self._data = new_data
    
    @property
    def functional_annotation_present(self) -> bool:
        return self._functional_db_format != None
     
    @property
    def taxonomy_db_format(self) -> TaxonomyFormat | None:
        return self._taxonomy_db_format # type: ignore
                
    @property
    def functional_db_format(self) -> FuncAnnotFormat | None:
        return self._functional_db_format # type: ignore
    
    @property
    def db_search_format(self) -> DbSearchSource | None:
        return self._db_search_format # type: ignore
    
    @property
    def db_search_confidence_format(self) -> DbSearchConfFormat | None:
        return self._db_search_confidence_format # type: ignore
    
    @property
    def de_novo_format(self) -> DeNovoSource | None:
        return self._de_novo_format # type: ignore
    
    @property
    def de_novo_confidence_format(self) -> DeNovoConfFormat | None:
        return self._de_novo_confidence_format # type: ignore
    
    @property
    def experiment_name(self) -> str | None:
        return self._experiment_name
    
    @property
    def source_files(self) -> List[str]:
        return self._source_files

    @property
    def get_metadata(self) -> Dict[str, Any]:
        """Return metadata from class object as dictionary. Metadata comprise
        all instance attributes outside of the data table itself.

        Returns:
            Dict[str, Any]: Metadata as dictionary
        """
        return {
            'Taxonomy DB Format': self.taxonomy_db_format,
            'Functional DB Format': self.functional_db_format,
            'DB Search Format': self.db_search_format,
            'DB Search Confidence Format': self.db_search_confidence_format,
            'De Novo Format': self.de_novo_format,
            'De Novo Confidence Format': self.de_novo_confidence_format,
            'Experiment Name': self.experiment_name
        }

    @property
    def __get_format_data(self) -> Dict[str, Any]:
        """Return specific metadata fields for comparing compatibility
        of dataset against different MetaPepTable dataset.

        Returns:
            Dict[str, Any]: Metadata as dictionary
        """
        return {
            'Taxonomy DB Format': self.taxonomy_db_format,
            'Functional DB Format': self.functional_db_format,
            'DB Search Format': self.db_search_format,
            'DB Search Confidence Format': self.db_search_confidence_format,
            'De Novo Format': self.de_novo_format,
            'De Novo Confidence Format': self.de_novo_confidence_format,
        }
        
    @property
    def sample_names(self) -> List[str]:
        return self._sample_names
    
    @classmethod
    def validate_json(cls, json_str: str) -> Tuple[bool, str]:
        """Check if json data represents a valid MetaPepTable object.

        Args:
            json_str (str): json data.

        Returns:
            Tuple[bool, str]: True if data is valid
        """
        # Deserialize the JSON string
        try:
            json_dict = json.loads(json_str)
        except Exception as err:
            return (False, 
                    "Failed to parse data as json.")

        # Extract the DataFrame and custom variables
        csv_str = decompress_string(json_dict['dataframe'])
        try:
            df = pd.read_csv(io.StringIO(csv_str), index_col=0, low_memory=False)
        except Exception as err:
            return (False,
                    "Failed to read dataset as DataFrame.")
        
        df = convert_deprecated_metapeptable_naming(df)
        valid, msg = cls.validate_input(df)
        
        if not valid:
            return (False, f"Invalid data format: {msg}")
        
        metadata = json_dict['metadata']
        valid_metadata = all([x in cls.OBJ_ARGS for x in metadata.keys()])

        return (valid_metadata, None)


    @classmethod
    def read_json(cls, json_str: str) -> Self:
        # Deserialize the JSON string
        json_dict = json.loads(json_str)

        # Extract the DataFrame and custom variables
        csv_str = decompress_string(json_dict['dataframe'])
        df = pd.read_csv(io.StringIO(csv_str), index_col=0, low_memory=False)
        df = convert_deprecated_metapeptable_naming(df)
        
        taxonomy_db_format = json_dict['metadata']['Taxonomy DB Format']
        functional_db_format = json_dict['metadata']['Functional DB Format']
        db_search_format = json_dict['metadata']['DB Search Format']
        db_search_confidence_format = json_dict['metadata']['DB Search Confidence Format']
        de_novo_format = json_dict['metadata']['De Novo Format']
        de_novo_confidence_format = json_dict['metadata']['De Novo Confidence Format']
        experiment_name = json_dict['metadata']['Experiment Name']
        
        return cls(df, taxonomy_db_format, functional_db_format, db_search_format,
                   db_search_confidence_format, de_novo_format,
                   de_novo_confidence_format, experiment_name)


    def to_json(self) -> str:
        """Write Metapep object to json format and store at file location
        """
        return to_json(self.data, self.get_metadata)


    @classmethod
    def concat_tables(cls, metapep_list: Sequence[Self]) -> Self:
        """Combine list of metapep objects into single object. It combines tables
        derived from multiple db search output files.
        
        While multiple input files can be freely combined, the resulting 
        MetaPepTable should adhere to the same metadata formats as all the input
        objects.

        Args:
            metapep_list (Sequence[metapep_obj]): List of MetaPepTable objects.

        Returns:
            MetaPepTable: Combined data object.
        """
        # check that input is not empty
        if len(metapep_list) == 0:
            raise ValueError("empty input sequence. Cannot generate MetaPepTable")
        
        
            
        # Set initial table for comparison of formats
        fval = metapep_list[0]
        fval_format_data = list(fval.__get_format_data.values())
        
        # Fetch format data from tables, look for conflicts, ingoring undefined formats
        consensus_format = fval.__get_format_data
        for table in metapep_list:
            table_formats = table.__get_format_data
            for format in table_formats.keys():
                if table_formats[format] is None:
                    continue
                elif consensus_format[format] is None:
                    consensus_format[format] = table_formats[format]
                if consensus_format[format] != table_formats[format]:
                    raise ValueError(f"Conflicting format during metapep table concatenation: {format}")
        
        # check that all spectral files are unique across objects
        validate_spectral_redundancy(metapep_list)

        # concatinate rows of all separate tables in one
        comb_df = pd.concat([i.data for i in metapep_list], axis=0).reset_index(drop=True)
        
        return cls(comb_df,
                   fval.taxonomy_db_format,
                   fval.functional_db_format,
                   fval.db_search_format,
                   fval.db_search_confidence_format,
                   fval.de_novo_format,
                   fval.de_novo_confidence_format,
                   fval.experiment_name)


    def remove_samples(self, samples: str | List[str]) -> Self:
        """Remove samples from metapep table dataset. This will filter all
        data from those samples from the data table and update format
        information by removing format values for which no more data is present
        in the dataset.

        Args:
            samples (str | List[str]): Sample names to remove from dataset

        Returns:
            Self: MetaPepTable object with sample data filtered out
        """
        
        if isinstance(samples, str):
            samples = [samples]
        
        # remove samples from data
        peptides_df = self.data
        peptides_df = peptides_df[~peptides_df["Sample Name"].isin(samples)]
        
        # update data tables
        self.data = peptides_df
        self._source_files = self.data['Source File'].dropna().unique().tolist()
        self._sample_names = self.data['Sample Name'].dropna().unique().tolist()
        
        def empty_or_absent(df: pd.DataFrame,
                            col: str) -> bool:
            if col not in df.columns:
                return True
            return df[col].dropna().empty
        
        # remove format information if no more data in dataset
        if (self.data["De Novo Imported"] == False).all():
            self._de_novo_format = None
            self._de_novo_confidence_format = None
        if (self.data["DB Search Imported"] == False).all():
            self._db_search_format = None
            self._db_search_confidence_format = None
        if empty_or_absent(self.data, "Taxonomy Id") and \
            empty_or_absent(self.data, "Global Taxonomy Id"):
            self._taxonomy_db_format = None
        if empty_or_absent(self.data, "KEGG_ko"):
            self._functional_db_format = None
        
        return self
        

db_search_T = TypeVar('db_search_T', bound='MetaPepDbSearch')

class MetaPepDbSearch(DataValidator):
    """Custom table format for DB search PSM data. Each spectrum match
    corresponds towards a single row in the dataset. Proteomics files from
    different sources can be converted into this format to be processed
    consistently within the dashboard.
    
    Multiple db search files can be concatenated. However, they need to be
    from the same input source and belong to a single sample.
    
    Global metadata:
        DATA_SOURCE     {'Peaks 11', 'Peaks 10', 'MaxQuant', 'ProteomeDiscoverer'}
        CONFIDENCE_FORMAT   {'-10lgp', 'p'}
        SAMPLE_NAME
    
    Fields:
        Peptide
        Sequence
        RT
        Scan
        m/z
        Charge
        ppm
        Length
        Confidence
        Area
        Mass
        Accession       # ';' delimited
        PTM
        Source File
    """
    
    ACCESSION_DELIMITER: str = ';'
    
    PEPTIDE_GROUP_NAME = 'PSM Count'
    
    REQUIRED_FIELDS = ['Peptide', 'Sequence', 'RT', 'Scan', 'm/z', 'Charge', 'ppm',
              'Length', 'Confidence', 'Area', 'Mass', 'Accession', 'PTM', 
              'Source File']
    NUMERIC_FIELDS = ['RT', 'Scan', 'm/z', 'Charge', 'ppm', 'Length', 
                      'Confidence', 'Area', 'Mass']
    
    # aggregation settings for grouping peptides
    AGGS_METHODS: Dict[str, str | Callable] = {
        'Peptide': mode_func,
        # 'm/z': 'mean',
        'Area': 'sum',
        'Confidence': 'max',    # TODO: Confidence merge depends on confidence format
        # 'Mass': 'mean',
        # 'ppm': 'mean'
        }
    
    # define parameters whose value should be taken based on the max value of another parameter
    MATCH_IDXMAX: Dict[str, str | List[str]] = {'Confidence': ['RT', 'm/z', 'Mass', 'ppm', 'Scan', 'Source File']}
    

    def __init__(self,
                 data: pd.DataFrame,
                 data_source: DbSearchSource,
                 confidence_format: DbSearchConfFormat,
                 source_files: Sequence[str],
                 sample_name: str | None):
        success, msg = MetaPepDbSearch.validate_input(data)

        if success is False:
            raise ValueError(msg)

        # extra columns in input data are not caught in validation, but dropped
        data = data[self.REQUIRED_FIELDS]
        
        # add data
        self._data = data
        self._data_source = data_source
        self._confidence_format = confidence_format
        self._source_files = source_files
        self._sample_name = sample_name
    
    @property
    def data(self) -> pd.DataFrame:
        return self._data
    
    @data.setter
    def data(self, new_data: pd.DataFrame):
        # test that converted data follows input format
        valid, msg = self.validate_input(new_data)
        
        if not valid:
            raise ValueError(msg)
        
        self._data = new_data
            
        
    @property
    def data_source(self) -> DbSearchSource:
        return self._data_source # type: ignore

    @property
    def confidence_format(self) -> DbSearchConfFormat:
        return self._confidence_format # type: ignore
        
    @property
    def sample_name(self) -> str | None:
        return self._sample_name

    @property
    def source_files(self) -> Sequence[str]:
        return self._source_files
    
    @property
    def source_files_set(self) -> Set[str]:
        return set(self._source_files)


    @property
    def get_metadata(self) -> Dict[str, Any]:
        """Return metadata from class object as dictionary. Metadata comprise
        all instance attributes outside of the data table itself.

        Returns:
            Dict[str, Any]: Metadata as dictionary
        """
        return {
            'Data Source': self.data_source,
            'Confidence Format': self.confidence_format,
            'source_files': self.source_files,
            'Sample Name': self.sample_name
        }


    @classmethod
    def read_json(cls, json_str: str) -> Self:
        # Deserialize the JSON string
        json_dict = json.loads(json_str)

        # Extract the DataFrame and custom variables
        # df = pd.DataFrame(json_dict['dataframe'])
        csv_str = decompress_string(json_dict['dataframe'])
        df = pd.read_csv(io.StringIO(csv_str))
        
        data_source = json_dict['metadata']['Data Source']
        confidence_format = json_dict['metadata']['Confidence Format']
        source_files = json_dict['metadata']['source_files']
        sample_name = json_dict['metadata']['Sample Name']
        
        return cls(df, data_source, confidence_format, source_files, sample_name)


    def to_json(self) -> str:
        """Write Metapep object to json format and store at file location
        """
        return to_json(self.data, self.get_metadata)
    

    def filter_spectral_name(self, spectral_name: str) -> Self | None:
        """Create a new MetaPepDbSearch instance where the dataset is filtered to
        only contain data specified towards a single source file.

        Args:
            spectral_name (str): Name of Source File, or raw file.
        """
        if spectral_name not in self.source_files_set:
            return None
        else:
            filtered_data = self.data[self.data['Source File'] == spectral_name]
            
            return self.__class__(
                filtered_data,
                self.data_source,
                self.confidence_format,
                [spectral_name],
                self.sample_name
            )
            
            
    @staticmethod
    def concat_tables(metapep_list: Sequence[db_search_T]) -> db_search_T:
        return _concat_tables_sample(metapep_list)


de_novo_T = TypeVar('de_novo_T', bound='MetaPepDeNovo')

class MetaPepDeNovo(DataValidator):
    """Custom table format for de novo identification data. Each peptide 
    identification corresponds towards a single row in the dataset.
    Proteomics files from different sources can be converted into this format
    to be processed consistently within the dashboard.
    
    Fields:
        Peptide
        Sequence
        RT
        Scan
        m/z
        Charge
        ppm
        Length
        Confidence
        Area
        Mass
        PTM
        Source File
    """
    PEPTIDE_GROUP_NAME = 'Match Count'
    
    REQUIRED_FIELDS = ['Peptide', 'Sequence', 'RT', 'Scan', 'm/z', 'Charge', 'ppm',
              'Length', 'Confidence', 'Area', 'Mass', 'PTM', 
              'Source File']
    NUMERIC_FIELDS = ['RT', 'Scan', 'm/z', 'Charge', 'ppm', 'Length', 
                      'Confidence', 'Area', 'Mass']
    
    # aggregation settings for grouping peptides
    AGGS_METHODS: Dict[str, str | Callable] = {
        'Peptide': mode_func,
        # 'm/z': 'mean',
        'Area': 'sum',
        'Confidence': 'max',    # TODO: Confidence merge depends on confidence format
        # 'Mass': 'mean',
        # 'ppm': 'mean'
        }
    MATCH_IDXMAX: Dict[str, str | List[str]] = {'Confidence': ['RT', 'm/z', 'Mass', 'ppm', 'Scan', 'Source File']}
    
    
    def __init__(self,
                 data: pd.DataFrame,
                 data_source: DeNovoSource,
                 confidence_format: DeNovoConfFormat,
                 source_files: Sequence[str],
                 sample_name: str | None):
        success, msg = MetaPepDeNovo.validate_input(data)

        if success is False:
            raise ValueError(msg)

        # extra columns in input data are not caught in validation, but dropped
        data = data[self.REQUIRED_FIELDS]
        
        # add data
        self._data = data
        self._data_source = data_source
        self._confidence_format = confidence_format
        self._source_files = source_files
        self._sample_name = sample_name


    @property
    def data(self) -> pd.DataFrame:
        return self._data
    
    @data.setter
    def data(self, new_data: pd.DataFrame):
        # test that converted data follows input format
        valid, msg = self.validate_input(new_data)
        
        if not valid:
            raise ValueError(msg)
        
        self._data = new_data
        
    @property
    def data_source(self) -> DeNovoSource:
        return self._data_source # type: ignore

    @property
    def confidence_format(self) -> DeNovoConfFormat:
        return self._confidence_format # type: ignore
        
    @property
    def sample_name(self) -> str | None:
        return self._sample_name

    @property
    def source_files(self) -> Sequence[str]:
        return self._source_files
    
    @property
    def source_files_set(self) -> Set[str]:
        return set(self._source_files)

    @property
    def get_metadata(self) -> Dict[str, Any]:
        """Return metadata from class object as dictionary. Metadata comprise
        all instance attributes outside of the data table itself.

        Returns:
            Dict[str, Any]: Metadata as dictionary
        """
        return {
            'Data Source': self.data_source,
            'Confidence Format': self.confidence_format,
            'source_files': self.source_files,
            'Sample Name': self.sample_name
        }    


    @classmethod
    def read_json(cls, json_str: str) -> Self:
        # Deserialize the JSON string
        json_dict = json.loads(json_str)

        # Extract the DataFrame and custom variables
        # df = pd.DataFrame(json_dict['dataframe'])
        csv_str = decompress_string(json_dict['dataframe'])
        df = pd.read_csv(io.StringIO(csv_str))
        
        data_source = json_dict['metadata']['Data Source']
        confidence_format = json_dict['metadata']['Confidence Format']
        source_files = json_dict['metadata']['source_files']
        sample_name = json_dict['metadata']['Sample Name']
        
        return cls(df, data_source, confidence_format, source_files, sample_name)


    def to_json(self) -> str:
        """Write Metapep object to json format and store at file location
        """
        return to_json(self.data, self.get_metadata)

    
    def filter_spectral_name(self, spectral_name: str) -> Self | None:
        """Create a new metapepnovo instance where the dataset is filtered to
        only contain data specified towards a single source file.

        Args:
            spectral_name (str): Name of Source File, or raw file.
        """
        if spectral_name not in self.source_files_set:
            return None
        else:
            filtered_data = self.data[self.data['Source File'] == spectral_name]
            
            return self.__class__(
                filtered_data,
                self.data_source,
                self.confidence_format,
                [spectral_name],
                self.sample_name
            )
    
    
    def filter_de_novo_only(self,
                            db_search_object: MetaPepDbSearch,
                            ignore_ptm: bool = True) -> Self:
        """Filter out any peptide in the dataset that is present in the
        db search object supplied.

        Args:
            db_search_object (MetaPepDbSearch): MetaPepDbSearch object.
            ignore_ptm (bool, optional): If true filter out any peptide whose
                amino acid sequence is identified in the db search file.
                Otherwise, filter peptides (including ptm's) that match  to the
                db search file. Defaults to True.

        Returns:
            Self: Filtered MetaPepDeNovo object.
        """
        match_col = 'Sequence' if ignore_ptm is True else 'Peptide'
        filter_peptides = db_search_object.data[match_col].to_list()
        
        filtered_de_novo = self.data[~self.data[match_col].isin(filter_peptides)]

        return self.__class__(
            filtered_de_novo,
            self.data_source,
            self.confidence_format,
            self.source_files,
            self.sample_name
        )
    
    
    @staticmethod
    def concat_tables(metapep_list: Sequence[de_novo_T]) -> de_novo_T:
        return _concat_tables_sample(metapep_list)
        

        
metapep_obj = TypeVar('metapep_obj', MetaPepDbSearch, MetaPepDeNovo)


def _concat_tables_sample(metapep_list: Sequence[metapep_obj]) -> metapep_obj:
    """Combine list of metapep objects into single object. It combines tables
    derived from multiple db search output files.
    
    While multiple input files can be freely combined, the resulting 
    MetaPepDbSearch or MetaPepDeNovo should represent a single 'experiment'
    (e.g. analysis of a single sample). For successful concatenation, metadata
    must be compatible for all objects. Data format, confidence format and
    sample name must be identical. In addition, redundancy of identifications
    should be avoided. Therefore, all raw data files between the input data
    objects should have no duplicates.

    Args:
        metapep_list (Sequence[metapep_obj]): List of MetaPep objects,
            belonging to the same experiment.

    Returns:
        metapep_obj: Combined data object.
    """
    # check that input is not empty
    if len(metapep_list) == 0:
        raise ValueError("empty input sequence. Cannot generate metapep table")
    
    # check that formats and samples are identical
    fval = metapep_list[0]
    if not all(
        (i.data_source, i.confidence_format, i.sample_name) == 
        (fval.data_source, fval.confidence_format, fval.sample_name) 
        for i in metapep_list
        ):
        raise ValueError("metapep metadata not consistent. Cannot concatenate.")
    
    # check that all spectral files are unique
    concat_spectra = validate_spectral_redundancy(metapep_list)
    
    # concatinate rows of all separate tables in one, reset index to ensure uniqueness
    comb_df = pd.concat([i.data for i in metapep_list], axis=0).reset_index(drop=True)
    
    return metapep_list[0].__class__(comb_df,
                                     fval.data_source, # type: ignore
                                     fval.confidence_format, # type: ignore
                                     concat_spectra,
                                     fval.sample_name)


def validate_spectral_redundancy(metapep_list: Sequence[MetaPepDeNovo | MetaPepDbSearch | MetaPepTable]):
    """Check from a list of MetaPep objects that the list of spectral source files
    are exlusively unique from each other. This ensures that during concatenation
    of objects, no redundancy in measurements is encountered.

    Args:
        metapep_list (Sequence[MetaPepDeNovo  |  MetaPepDbSearch  |  MetaPepTable]):
            List of metapep objects to check spectral redundancy

    Raises:
        ValueError: Redundant spectral sources found across MetaPep objects.
    """
    # check that all spectral files are unique across objects
    union_list = chain.from_iterable([i.source_files for i in metapep_list])
    union_list = list(union_list)
    
    if len(set(union_list)) < len(union_list):
        file_format = "DB search" if isinstance(metapep_list[0], MetaPepDbSearch) else "de novo"
        raise ValueError(f"Redundancy found in {file_format} files, cannot concatenate.")
    
    return union_list