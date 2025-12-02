"""This module contains functions that process MetaPep custom objects.
It is meant for functions that are in dependency hierarchy above the 'types'
objects. This separates it from 'utils's functions, which comprise only basic
functions at the bottom of the hierarchy."""

from typing import IO, Dict, Callable, Type, Union
from copy import deepcopy

from metapepview.backend.utils import custom_groupby, re_find_list
from metapepview.backend.types import *
from metapepview.backend.types.object_mappings import db_search_importers, de_novo_importers


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
        except Exception as e:
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
        except Exception as e:
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
    try:
        importer = db_search_importers.get(file_format)
        
        if importer is None:
            raise ValueError("Invalid file format supplied")
        else:
            db_search_obj = importer.read_file_buffer(file_buffer, sample_name)
        
        # convert db search data to metapep format. If multiple source files in dataset,
        # omit all data that are not of the currently processed source file.
        return db_search_obj.to_metapep_db_search(sample_name, 
                                                  crap_dataset)

    except Exception as err:
        raise ValueError(f"Failed to load db search sample '{sample_name}': {err}")
    
    
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
    try:
        if importer is None:
            raise ValueError("Invalid file format supplied")
        else:
            de_novo_obj = importer.read_file_buffer(file_buffer, sample_name)
        
        # convert de novo data to metapep format. If multiple source files in dataset,
        # omit all data that are not of the currently processed source file.
        return de_novo_obj.to_metapep_de_novo(sample_name, crap_dataset)
    except Exception as err:
        raise ValueError(f"Failed to load de novo sample '{sample_name}': {err}")


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


def metapep_de_novo_to_peptides(metapep_de_novo: MetaPepDeNovo,
                                aggs_methods: Dict[str, str | Callable] | None=None,
                                match_idxmax: Dict[str, str | List[str]]| None=None,
                                group_size_name: str = 'PSM Count') -> pd.DataFrame:
    """Group spectral scans into peptide sequence rows. Column values are aggregated
    by specific functions that best describe the peptide.
    
    Note:
        Custom aggregation methods described will overwrite all defaults.

    Args:
        metapep_de_novo (MetaPepDeNovo): de novo dataset.
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
    return custom_groupby(metapep_de_novo.data,
                          groupby_col,
                          aggs_methods,
                          match_idxmax,
                          group_size_name=group_size_name)


def metapep_db_search_accession_regex(
        metapep_db_search: MetaPepDbSearch,
        acc_column: str,
        acc_regex: str | None) -> MetaPepDbSearch:
    """Perform regex operation on protein id's present in MetaPepDbSearch
    dataset. It matches a pattern to each protein id and convert the full id
    to the first capture group value from the regex pattern.

    Args:
        metapep_db_search (MetaPepDbSearch): db search object.
        acc_regex (str | None): regex pattern.

    Returns:
        MetaPepDbSearch: db search object with converted protein id's.
    """
    if acc_regex != "" and acc_regex is not None:
        acc_delim = GlobalConstants.peptides_accession_delimiter
        peptides = metapep_db_search.data
        
        valid_acc = ~peptides[acc_column].isna()
        # split accessions per row, apply regex to each element, and join again.
        peptides.loc[valid_acc, acc_column] = peptides.loc[valid_acc, acc_column].str.split(acc_delim)\
            .apply(re_find_list, pattern = acc_regex)\
            .apply(lambda x: acc_delim.join(x))
        
        metapep_db_search.data = peptides
    
    return metapep_db_search


T = TypeVar("T", bound=MetaPepDbSearch | MetaPepDeNovo)
def rt_from_spectral_data(metapep_table: T, 
                          spectral_df: pd.DataFrame) -> T:
    """In case when retention time is not reported in MetaPep db search or 
    de novo dataset. Attempt to obtain retention time information by mapping
    PSM's to spectral information from the scan number value.

    Note:
        This function will overwrite the complete retention time column with 
        spectral retention time information.

        Depending on the proteomics data source, numbers reported in proteomics 
        datasets might not correspond directly to the scan numbers defined in 
        the mzML dataset. Validation can be performed by matching precursor m/z
        as reported in the mzML dataset to the measure m/z value reported in 
        proteomics data.

    Args:
        metapep_table (MetaPepDbSearch | MetaPepDeNovo): Proteomics dataset.
        spectral_df (pd.DataFrame): mzML dataset

    Returns:
        MetaPepDbSearch | MetaPepDeNovo: Proteomics dataset supplemented with 
            retention time data.
    """
    df = deepcopy(metapep_table.data)

    # check if scan numbers are reported in metapep dataset
    if df["Scan"].isnull().all():
        return metapep_table
    
    # map scan number to retention time in mzml file
    scan_to_rt = spectral_df.set_index("scan number")["retention time"].to_dict()

    df.loc[:, "RT"] = df["Scan"].apply(lambda x: scan_to_rt.get(x, np.nan))

    metapep_table.data = df

    return metapep_table
