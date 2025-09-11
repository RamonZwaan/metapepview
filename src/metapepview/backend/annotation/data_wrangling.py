import io

from typing import List, Dict, Tuple, Any
import pandas as pd

from metapepview.backend.utils import *
from metapepview.backend.type_operations import *
from metapepview.backend.types import *
from metapepview.constants import *


def match_source_denovo(denovo_list: Sequence[IO[str] | str],
                        denovo_format: DeNovoSource,
                        crap_dataset: pd.Series | None = None) -> Dict[str, MetaPepDeNovo]:
    """Extract source files from denovo datasets and match them to the datasets
    by mapping index to source files. If multiple denovo files point to same source
    file, then only the first analyzed file is taken. During the matching, data
    is imported into table files.

    Args:
        denovo_list (Sequence[IO[str]] | List[str]): list of denovo datasets.
        denovo_format (DeNovoSource): Data format of input datasets.
        crap_dataset (pd.Series | None, optional): Series of peptide sequences
            part of the cRAP dataset. These are filtered out of the de novo 
            peptide data. defaults to None.

    Returns:
            Dict[str, MetaPepDeNovo]: Dict of de novo source and MetaPepDeNovo data.
    """
    # Store source files as keys into dictionary with index in denovo_list as value
    denovo_sources = {}
    for denovo_contents in denovo_list:
        de_novo_obj = load_metapep_de_novo(denovo_contents,
                                           None,
                                           denovo_format,
                                           crap_dataset)
        raw_files = de_novo_obj.source_files
        for spec_file in raw_files:
            # skip if raw file already present in the dictionary
            if spec_file in denovo_sources.keys():
                continue
            # if multiple raw source files in de novo data, keep only data from selected source
            denovo_sources[spec_file] = de_novo_obj.filter_spectral_name(spec_file)
    
    return denovo_sources


def match_source_denovo_buffer(denovo_list: Sequence[IO[str] | str],
                               source_file_col: int=0) -> Dict[str, IO[str]]:
    """Extract source files from denovo datasets and match them to the datasets
    by mapping index to source files. If multiple denovo files point to same source
    file, then only the first analyzed file is taken.

    Args:
        denovo_list (Sequence[IO[str]] | Sequence[str]): list of denovo datasets.
        source_file_col (int, optional): column index that contains source file
            name. Defaults to 0.
            
    Returns:
            Dict[str, IO[str]]: Dict of de novo source and file data.
    """
    # Store source files as keys into dictionary with index in denovo_list as value
    denovo_sources = {}
    for i, denovo_contents in enumerate(denovo_list):
        
        # convert data to string file-like buffer
        if isinstance(denovo_contents, str):
            denovo_file = memory_to_stringio(denovo_contents)
        else:
            denovo_file = denovo_contents
        
        # get header row, if no data, print issue and continue to next
        header_line = next(denovo_file, None)
        if header_line is None:
            print("Failed to import denovo data / no data inside file")
            break
        
        header_line = header_line.replace('"', '')
        # before getting data, assert that file is correctly formatted
        if header_line.split(",")[source_file_col] != "Source File":
            print("Format error in denovo data. No 'Source File' column")
        
        # get first line of values, if no data, that means file is empty
        line = denovo_file.readline()
        if line == "" or line == "\n":
            print("Failed to import denovo data / no data inside file")
            break
        
        # return cursor in filebuffer to start of file
        denovo_file.seek(0)
        
        # get source column and match against dict
        source = line.split(",")[source_file_col].replace('"', '')
        if source in denovo_sources.keys():
            print("duplicate sources encountered, discard duplicate")
            continue
        denovo_sources[source] = denovo_file
    
    return denovo_sources


def deduplicate_input_lists(current_sample_names: List[str],
                            sample_name: str,
                            psm_data: Sequence[Any],
                            psm_names: List[str],
                            filter_duplicates: bool = False) -> Tuple[str, Sequence[Any], List[str]]:
    """Preprocessing of input data lists prior to annotation and concatenation into
    peptides dataset. This function does not modify db search data itself, it deduplicates
    input file names or remove duplicates to obtain a input dataset that can be
    added onto the peptides dataset.

    Args:
        current_sample_names (List[str]): List of names in current peptides dataset.
        sample_name (str): Name of new input data.
        psm_data (Sequence[Any]): input psm data.
        psm_names (List[str]): filenames of input data.
        filter_duplicates (bool): Specify if duplicate files should be filtered
            out of the input data. If False, filenames will be deduplicated.
            Defaults to False.
    
    Returns:
        str: Deduplicated sample name.
        Sequence[Any]: Filtered psm input data.
        List[str]: Filtered/Deduplicated psm input filenames.
    """
    # deduplicate sample name from current samples
    sample_name = deduplicate_strings(sample_name, current_sample_names)

    # If specified, filter out files present in the current peptide df
    # (based on psm file name)
    if filter_duplicates is True:
        psm_data, psm_names = filter_sample_duplicates(current_sample_names, 
                                                       psm_data,
                                                       psm_names)

    # if no filtering selected, locate collision files and deduplicate them
    else:
        dupl_present, col_names, col_idx = check_sample_duplicates(current_sample_names, 
                                                                   psm_names)
        
        # deduplicate psm filenames, by adding '(n)' suffix
        # copy sample names to new variable to append only for this operation
        psm_sample_names = current_sample_names
        
        for psm_name, idx in zip(col_names, col_idx):
            unique_name = deduplicate_strings(psm_name, psm_sample_names)
            psm_names[idx] = unique_name
            psm_sample_names.append(unique_name)

    return sample_name, psm_data, psm_names


def filter_sample_duplicates(current_sample_names: List[str],
                             psm_data: Sequence[Any],
                             psm_names:List[str]) -> Tuple[Sequence[Any], List[str]]:
    """Filter psm input data that contain identical names to sample names in
    existing peptides dataset.

    Args:
        current_sample_names (List[str]): sample names from current annotated peptides dataset
        psm_data (Sequence[Any]): List of data for input db search PSM files
        psm_names (List[str]): List of filenames for input db search PSM files
    
    Returns:
        Sequence[Any]: Filtered list of psm sample data.
        List[str]: Filtered list of psm names.
    """
    # fetch indices of samples that are new
    new_samples = []
    for i, sample in enumerate(psm_names):
        if sample not in current_sample_names:
            new_samples.append(i)
            
    # filter psm data and psm names with the unique samples
    psm_data = [psm_data[i] for i in new_samples]
    psm_names = [psm_names[i] for i in new_samples]

    return psm_data, psm_names


def check_sample_duplicates(current_sample_names: List[str],
                            psm_names: List[str]) -> Tuple[bool, List[str], List[int]]:
    """Check if any name collission is present between input PSM data
    and existing sample names in the peptides dataset. Collisions include
    duplicates within input names, but these are only documented from the
    second occurrence onward.

    Args:
        current_sample_names (List[str]): Sample names from current annotated peptides dataset
        psm_names (List[str]): List of file names for PSM data

    Returns:
        bool: True if name collision encountered
        List[str]: list of filenames with collision
        List[int]: list of indices of collision files in psm_names list
    """
    # add name and indices of collisions to list
    dupl_names, dupl_idx = [], []
    for idx, name in enumerate(psm_names):
        if name in current_sample_names:
            dupl_names.append(name)
            dupl_idx.append(idx)
        # add current name to include potential collisions within new names
        else:
            current_sample_names.append(name)
    
    return (dupl_idx != [], dupl_names, dupl_idx)


# @lru_cache(maxsize=128)
def unique_taxa_from_rank(metapep_table: MetaPepTable,
                          rank: str,
                          filter_clade: str | None=None,
                          clade_rank: str | None=None,
                          show_names: bool=True) -> List[str | int]:
    """Return list of unique taxa present in the peptides dataset for a given
    taxonomy rank.

    Args:
        peptide_json (str): peptide dataset
        rank (str): taxonomy rank.
        filter_clade (str | None): taxonomy group to filter dataset from.
            Defaults to None.
        clade_rank (str | None): rank of filter clade. Defaults to None.
        show_names (bool): Return list of names as opposed to list of tax id's.
            Defaults to True.

    Returns:
        List[str | int]: list of unique taxa
    """
    dataset = metapep_table.data
    
    if show_names is True:
        col_suffix = ' Name'
    else:
        col_suffix = ' Id'
    
    # filter dataset to only contain rows belonging to one clade
    if filter_clade and clade_rank and clade_rank != 'Root':
        dataset = dataset[dataset[clade_rank + " Name"] == filter_clade]
    
    return dataset[rank + col_suffix]\
        .dropna()\
        .drop_duplicates()\
        .sort_values()\
        .to_list()