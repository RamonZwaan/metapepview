from typing import List, Any, Sequence, Pattern
import re
from itertools import chain
import numpy as np

from metapepview.constants import GlobalConstants


def peaks_10_drop_fraction_id(input_cell: str | float | None,
                              custom_regex: str | None = None) -> str | float | None:
    """Drop fraction identifier from value columns, such as scan id from Peaks
    Studio 10 db search psm file.
    
    Peaks incudes the fraction number to certain columns in experiments
    that process multiple MS runs. However, since the fraction can be inferred
    directly from the 'source file' column, it is redundant to include the
    fraction information in the other columns. In addition, it allows for more
    consistent data processing as column values stripped from fraction info can
    be directly used to match other datasets.

    Args:
        input_cell (str | float | None): input value.

    Returns:
        str: Value stripped from fraction data
    """
    # if value is numeric or empty, return value directly
    if not isinstance(input_cell, str):
        return input_cell
    
    # filter out fraction information
    if custom_regex is None:
        return input_cell.split(":")[-1]
    else:
        pattern = re.compile(custom_regex)
        value_match = pattern.search(input_cell)
        if value_match is None:
            return None
        else:
            return value_match.group(0)


def truncate_end(string: str | float | None, end: int) -> str | float | None:
    """Truncate end of string if it exeeds user defined character length.

    Args:
        string (str): Input string.
        end (int): Output string.

    Returns:
        str: _description_
    """
    if not isinstance(string, str):
        return string
    else:
        return string[:end] + "..." if len(string) > end else string


def deduplicate_strings(input_string: str,
                        string_set: Sequence[str]) -> str:
    """Add suffix to strings that have a duplicate in a given
    set of strings, iterate through suffixes until unique value
    encountered

    Args:
        input_string (str): _description_
        string_set (Sequence[str]): _description_

    Returns:
        _type_: _description_
    """
    suffix, counter = "", 0
    
    # iterate suffix count until unique value
    while input_string + suffix in string_set:
        counter += 1
        suffix = f"({counter})"
    
    return input_string + suffix



def wrangle_peptides(sequence: str | float,
                     ptm_filter: bool=True,
                     li_swap: bool=True,
                     ptm_regex: Pattern[str] = GlobalConstants.sequence_regex) -> str:
    """Process protein sequences by removing post-translational
    modifications and/or equating Leucin and Isoleucin amino acids.

    Args:
        sequence (str): Protein sequence string
        ptm_filter (bool, optional): Remove PTM,s from sequence.
            Defaults to True.
        li_swap (bool, optional): Equate leucin and isoleucin.
            Defaults to True.

    Returns:
        str: Processed sequence string.
    """
    # handle potential missing values
    if sequence != sequence or sequence is None:
        return sequence

    if ptm_filter is True:
        sequence = "".join(re.findall(ptm_regex, sequence))
    if li_swap is True:
        sequence = sequence.replace("L", "I")
    return sequence


def digest_proteins(sequence: str,
                    custom_cleave_rule: re.Pattern[str] | None = None,
                    min_length: int=7,
                    max_length: int=35,
                    mis_cleave: int=3,
                    m_cleave: bool=True,
                    li_swap: bool=True,
                    reverse: bool=False) -> List[str]:
    """Apply in-silico digestion to protein string. The function splits
    a string consisting of amino acids into a list of peptide strings.
    The argument settings allow for changes in the digestion rules.

    Args:
        protein (str): String of a protein sequence.
        custom_cleave_rule (Union[str, re.Pattern, None], optional): 
        min_length (int, optional): Minimal peptide length to return. Fragments
            below this length are discarded. Defaults to 7.
        max_length (int, optional): Maximal peptide length to return. Fragments
            above this length are discarded. Defaults to 35.
        mis_cleave (int, optional): Include peptides with n splits missed.
            Defaults to 0.
        m_cleave (bool, optional): Remove first Methionine from protein string.
            Defaults to True.
        li_swap (bool, optional): Replace Leucine with Isoleucine. Defaults to True.
        cleave_rule (tuple, optional): Recognition site and cleave side.
            Defaults to (r'[RK]', 'c').
        reverse (bool, optional): Reverse protein string. Defaults to False.

    Returns:
        list: List of strings representing peptide sequences.
    """
    # reverse sequence
    if reverse:
        sequence = sequence[::-1]
    
    # remove leading Met
    if m_cleave and sequence[0] == 'M':
        sequence = sequence[1:]

    # set cleave rule. If none specified, get trypsin rule.
    if custom_cleave_rule is None:
        cleave_rule = GlobalConstants.trypsin_cleave_rule
    else:
        cleave_rule = custom_cleave_rule

    #split sequence
    cut_pept = re.split(cleave_rule, sequence)
    
    # add miscleavages to list of peptide sequences
    join_elements = lambda seqs, n: [''.join(seqs[i:i+2+n]) for i in range(len(seqs) - (n+1))]
    mis_cleave_pept = [join_elements(cut_pept, i) for i in range(mis_cleave)]
    cut_pept = list(chain(cut_pept, *mis_cleave_pept))

    # replace L to I and filter length
    if li_swap:        
        cut_pept = [i.replace("L", "I") for i in cut_pept if min_length <= len(i) <= max_length]
    else:
        cut_pept = [i for i in cut_pept if min_length <= len(i) <= max_length]

    return cut_pept


def re_find_list(input_list: List[str] | float | None,
                 pattern: re.Pattern) -> List[str]:
    """Perform pattern extraction to each string element from a list of inputs.
    If no match is found, the original string is returned.

    Args:
        input_list (List[str]): Input strings.
        pattern (re.Pattern): Pattern to extract.

    Returns:
        List[str]: List of substrings extracted from input
    """
    if input_list != input_list or input_list is None:
        return input_list

    for idx, str_item in enumerate(input_list):
        res = re.search(pattern, str_item)

        if res is None:
            input_list[idx] = str_item
        else:
            input_list[idx] = res.group()
    
    return input_list
    


def spectrum_id_to_scan_number(spectrum_id: str) -> int | float:
    """Extract scan number from a spectrum id. This is useful when dealing with
    sources that only provide a spectrum id in the output data. The function 
    will take into account multiple vendor-specific formats.

    If no valid scan field is encountered, NaN is returned.

    Note:
        The difference in definitions between "spectrum" and "scan" may differ 
        based on the vendor format and used acquisition mode. Notably, in 
        timsTOF ion mobility spectral data, scans have a different meaning from
        regular acquisition methods.

    Args:
        spectrum_id (str): Spectrum id as provided from the spectral dataset.

    Returns:
        int | : Scan number or nan if no scan number found.
    """
    for pattern in (
        GlobalConstants.frame_pattern, 
        GlobalConstants.scan_pattern, 
        GlobalConstants.cycle_pattern, 
        GlobalConstants.spectrum_pattern):
        match_obj = re.search(pattern, spectrum_id)
        if match_obj is not None:
            return int(re.search(GlobalConstants.int_pattern, match_obj[0])[0])
    return np.nan