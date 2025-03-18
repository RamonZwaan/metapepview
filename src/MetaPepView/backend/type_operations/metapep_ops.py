"""This module contains functions that process MetaPep custom objects.
It is meant for functions that are in dependency hierarchy above the 'types'
objects. This separates it from 'utils's functions, which comprise only basic
functions at the bottom of the hierarchy."""

from typing import IO, Dict, Callable, Type, Union

from ..utils import custom_groupby
from ..types import *
from .object_mappings import db_search_importers, de_novo_importers


def validate_db_search(file_buffer: str | IO[str],
                       file_format: DbSearchSource) -> Tuple[bool, str | None]:
    """Check that user supplied file is parsed correctly following the specified
    format.

    Args:
        file_buffer (str | IO[str]): Input file buffer.
        file_format (DbSearchSource): Format of file contents

    Returns:
        Tuple[bool, str | None]: Parse successful and if not an error message.
    """
    importer = db_search_importers.get(file_format)
    
    if importer is None:
        return False, "Invalid file format supplied"
    else:
        try:
            importer.read_file_buffer(file_buffer)
            return True, None
        except ValueError as e:
            return False, repr(e)


def validate_de_novo(file_buffer: str | IO[str],
                     file_format: DeNovoSource) -> Tuple[bool, str | None]:
    """Check that user supplied file is parsed correctly following the specified
    format.

    Args:
        file_buffer (str | IO[str]): Input file buffer.
        file_format (DbSearchSource): Format of file contents

    Returns:
        Tuple[bool, str | None]: Parse successful and if not an error message.
    """
    importer = de_novo_importers.get(file_format)
    
    if importer is None:
        return False, "Invalid file format supplied"
    else:
        try:
            importer.read_file_buffer(file_buffer)
            return True, None
        except ValueError as e:
            return False, repr(e)


def load_metapep_db_search(file_buffer: str | IO[str],
                           sample_name: str,
                           file_format: DbSearchSource,
                           crap_dataset: pd.Series | None = None) -> MetaPepDbSearch:
    """Load buffer of db search file into specified type object and convert to
    MetaPepDbSearch object.

    Args:
        file_buffer (IO[str]): File buffer object
        sample_name (str): Name of sample.
        file_format (DbSearchSource): Format of db search psm output.

    Returns:
        MetaPepDbSearch: Db search psm data in MetaPep table format
    """
    importer = db_search_importers.get(file_format)
    
    if importer is None:
        raise ValueError("Invalid file format supplied")
    else:
        db_search_obj = importer.read_file_buffer(file_buffer, sample_name)
    
    # convert db search data to metapep format. If multiple source files in dataset,
    # omit all data that are not of the currently processed source file.
    return db_search_obj.to_metapep_db_search(sample_name, crap_dataset)\
        .filter_spectral_name(sample_name)
    
    
def load_metapep_de_novo(file_buffer: str | IO[str],
                         sample_name: str | None,
                         file_format: DeNovoSource,
                         crap_dataset: pd.Series | None = None) -> MetaPepDeNovo:
    """Load buffer of de novo file into specified type object and convert to
    MetaPepDeNovo object.

    Args:
        file_buffer (IO[str]): File buffer object
        sample_name (str): Name of sample.
        file_format (DeNovoSource): Format of db search psm output.
        crap_dataset (pd.Series | None, optional): Series of peptide sequences
            part of the cRAP dataset. These are filtered out of the de novo
            peptide data. defaults to None.

    Returns:
        MetaPepDeNovo: De novo data in MetaPep table format
    """

    importer = de_novo_importers.get(file_format)
    
    if importer is None:
        raise ValueError("Invalid file format supplied")
    else:
        de_novo_obj = importer.read_file_buffer(file_buffer, sample_name)
    
    # convert de novo data to metapep format. If multiple source files in dataset,
    # omit all data that are not of the currently processed source file.
    return de_novo_obj.to_metapep_de_novo(sample_name, crap_dataset)\
        .filter_spectral_name(sample_name)


def metapep_table_to_peptides(metapep_table: MetaPepDbSearch | MetaPepDeNovo,
                              aggs_methods: Dict[str, str | Callable] | None=None,
                              match_idxmax: Dict[str, str | List[str]]| None=None) -> pd.DataFrame:
    """Group spectral scans into peptide sequence rows. Column values are aggregated
    by specific functions that best describe the peptide.
    
    Note:
        Custom aggregation methods described will overwrite all defaults.

    Args:
        metapep_table (pd.DataFrame): metapep table object, db search or de novo.
        custom_aggs (Dict[str, Callable] | None, optional): Dictionary of column names and custom aggregation
            methods to overwrite the default methods. Defaults to None.
        match_idxmax (Dict[str | List[str]] | None, optional): Specify columns who'se value 
            selection is bound to the maximum group value of a different column.
            Defaults to None.
        group_size_name (str, optional): Name of column to store peptide group
            sizes. Defaults to 'Count'

    Returns:
        pd.DataFrame: Aggregated dataset.
    """
    groupby_col = "Sequence"

    # Setup defaults if no customs defined
    if match_idxmax is None:
        match_idxmax = metapep_table.MATCH_IDXMAX    
    if aggs_methods is None:
        aggs_methods = metapep_table.AGGS_METHODS
        
    return custom_groupby(metapep_table.data,
                          groupby_col,
                          aggs_methods,
                          match_idxmax,
                          group_size_name=metapep_table.PEPTIDE_GROUP_NAME)


def metapep_de_novo_to_peptides(metapep_db_search: MetaPepDeNovo,
                                aggs_methods: Dict[str, str | Callable] | None=None,
                                match_idxmax: Dict[str, str | List[str]]| None=None,
                                group_size_name: str = 'PSM Count') -> pd.DataFrame:
    """Group spectral scans into peptide sequence rows. Column values are aggregated
    by specific functions that best describe the peptide.
    
    Note:
        Custom aggregation methods described will overwrite all defaults.

    Args:
        db_search_psm (pd.DataFrame): db search psm dataset.
        custom_aggs (Dict[str, Callable] | None, optional): Dictionary of column names and custom aggregation
            methods to overwrite the default methods. Defaults to None.
        match_idxmax (Dict[str | List[str]] | None, optional): Specify columns who'se value 
            selection is bound to the maximum group value of a different column.
            Defaults to None.
        group_size_name (str, optional): Name of column to store peptide group
            sizes. Defaults to 'PSM Count

    Returns:
        pd.DataFrame: Aggregated dataset.
    """
    groupby_col = "Sequence"

    # Setup defaults if no customs defined
    if match_idxmax is None:
        match_idxmax = {'Area': ['RT', 'Scan']}
    
    if aggs_methods is None:
        mode_func = lambda x: pd.Series.mode(x).iat[0]
        aggs_methods = {'Peptide': mode_func,
                        'm/z': 'mean',
                        'Area': 'sum',
                        'Confidence': 'max',    # TODO: Confidence merge depends on confidence format
                        'Mass': 'mean',
                        'ppm': 'mean'}
    
    
    # perform groupby process
    return custom_groupby(metapep_db_search.data,
                          groupby_col,
                          aggs_methods,
                          match_idxmax,
                          group_size_name=group_size_name)
