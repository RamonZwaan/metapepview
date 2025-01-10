from typing import List, Tuple, Any, Sequence
import numpy as np
import pandas as pd

from backend.types import MetaPepDbSearch, MetaPepDeNovo
from backend.spectral_ref_builder.utils import count_threshold_values
from backend.utils import wrangle_peptides, filter_denovo_only


def reference_score_distribution(stat_dict: dict,
                                 formats: List[str],
                                 label_dict: dict[str, List[str]],
                                 normalize_psm: bool=False,
                                 normalize_rt: bool=False,
                                 normalize_fill: bool=False) -> pd.DataFrame:
    """Extract peptide identification confidence thresholds from reference
    dataset and return it as a dataframe. Confidence distributions from
    different identification methods can be extracted. In addition, confidence
    counts can be normalized sample wise towards total MS2 scans, or towards
    retention time.

    Args:
        stat_dict (dict): Reference dataset dictionary.
        formats (List[str]): Identification method to extract data of.
        normalize_psm (bool, optional): Divide counts by total MS2 scans for
            sample. Defaults to False.
        normalize_rt (bool, optional): Divide counts by total run time.
            Defaults to False.
        normalize_fill (bool, optional): For each sample, divide (potentially
            psm or rt normalized) counts such that the largest category counts
            to 100. This way, bars between samples are of same height.
            Defaults to False.

    Raises:
        ValueError: Invalid identification method given in 'formats'.

    Returns:
        pd.DataFrame: Threshold confidence distribution data from ref dataset.
    """
    # output variables
    dist_data = {"ident method": [],
                 "x axis": [],
                 "value": [],
                 "sample": []}
    
    format_to_key = {"db search": "db search confidence dist",
                     "de novo": "de novo confidence dist",
                     "de novo only": "de novo only confidence dist"}
    
    if any([i not in ["db search", "de novo", "de novo only"] for i in formats]):
        raise ValueError("invalid format supplied")

    # parse dictionary and add data
    for sample, data in stat_dict["samples"].items():
        # initialize data lists and add format values
        ident_method = []
        sample_name = []
        threshold_names = []
        counts = []
        for format in formats:
            format_key = format_to_key[format]
            
            # normalize db psm by the total number of scans
            if normalize_psm is True:
                if data["ms2 count"] != data["ms2 count"]:
                    continue
                psm_num = data["ms2 count"]
                category_counts = (np.array(data[format_key]["counts"]) / psm_num).tolist()
            # all sections can be normalized by retention time
            elif normalize_rt is True:
                if data['total rt'] != data['total rt']:
                    continue
                sample_rt = data['total rt']
                category_counts = (np.array(data[format_key]["counts"]) / sample_rt).tolist()
            else:
                category_counts = data[format_key]["counts"]
                
            # if normalize_fill is true, divide counts by largest number and multiply by 100
            if normalize_fill is True and format == "db search":
                category_counts = [(i / max(category_counts)) * 100 for i in category_counts]
                     
            counts += category_counts      
            prefix, suffix = label_dict[format]
            
            threshold_names += [f"{prefix} {i} {suffix}" for i in data[format_key]["thresholds"]]
            ident_method += [format_key]*len(data[format_key]["counts"])
            sample_name += [sample]*len(data[format_key]["counts"])

        dist_data["ident method"] += ident_method
        dist_data["x axis"] += threshold_names
        dist_data["value"] += counts
        dist_data["sample"] += sample_name
        
    return pd.DataFrame(dist_data)


def reference_score_dist_peaks(db_search_psm: MetaPepDbSearch | None,
                               de_novo: MetaPepDeNovo | None,
                               ms2_number: int | None,
                               total_rt: float | None,
                               formats: List[str],
                               db_search_tresholds: List[int],
                               de_novo_tresholds: List[int],
                               normalize_psm: bool=False,
                               normalize_rt: bool=False,
                               normalize_fill: bool=False,
                               sample_name: str | None=None) -> pd.DataFrame:
    exclude_formats = []
    # Check if function arguments are valid given the supplied input data, else update
    if db_search_psm is None:
        normalize_psm = False
        exclude_formats += ["db search", "de novo only"]
    if de_novo is None:
        exclude_formats += ["de novo", "de novo only"]
    if ms2_number is None:
        normalize_psm = False
    if total_rt is None:
        normalize_rt = False
        
    if sample_name is None:
        sample_name = "Sample"
        
    formats = list(set(formats).difference(exclude_formats))
    
    if any([i not in ["db search", "de novo", "de novo only"] for i in formats]):
        raise ValueError("invalid format supplied")
            
    # initialize data lists and add format values
    counts = []
    ident_method = []
    threshold_names = []
    sample_names = []
    for data_format in formats:
        # compute values for specific format
        if data_format == "db search" and db_search_psm is not None:
            score_format = db_search_psm.confidence_format
            ax_suffix = ""
            format_key = "db search confidence dist"
            df = db_search_psm.data
            thresholds = db_search_tresholds
        elif data_format == "de novo" and de_novo is not None:
            score_format = de_novo.confidence_format
            ax_suffix = ""
            format_key = "de novo confidence dist"
            df = de_novo.data
            thresholds = de_novo_tresholds
        elif de_novo is not None and db_search_psm is not None:
            score_format = de_novo.confidence_format
            ax_suffix = "d-only"
            format_key = "de novo only confidence dist"
            df = filter_denovo_only(de_novo.data, db_search_psm.data)
            thresholds = de_novo_tresholds
        else:
            raise ValueError("incorrect format supplied or metaproteomics data missing.")
        
        format_counts, thresholds = count_threshold_values(df,
                                                           "Confidence",
                                                           thresholds,
                                                           True)
            
        # normalize db psm by the total number of psm's
        if normalize_psm is True:# and format_key == "-10lgP dist":
            category_counts = (np.array(format_counts) / ms2_number).tolist()
        # all sections can be normalized by retention time
        elif normalize_rt is True:
            category_counts = (np.array(format_counts) / total_rt).tolist()
        else:
            category_counts = format_counts
        
        # if normalize_fill is true, divide counts by largest number and multiply by 100
        if normalize_fill is True and data_format == "db search":
            category_counts = [(i / max(category_counts)) * 100 for i in category_counts]

        # add normalized counts to trace dataset
        counts += category_counts
        
        threshold_names += [f"{score_format} {i} {ax_suffix}" for i in thresholds]
        ident_method += [format_key]*len(format_counts)
        sample_names += [sample_name]*len(format_counts)
    
    # output variables
    dist_data = {"ident method": ident_method,
                 "x axis": threshold_names,
                 "value": counts,
                 "sample": sample_names}
    
    return pd.DataFrame(dist_data)
    