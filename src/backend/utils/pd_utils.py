"""This module contains helper functions related to pandas df processing,
Here, string processing is omitted as it is part of the string_utils
"""

from typing import List, Tuple, Dict, Any
from pathlib import Path
import json
import io

import pandas as pd
import numpy as np

from constants import *
from .string_utils import digest_proteins, wrangle_peptides
from .io_utils import decompress_string


# return the modus of a pandas series, only the most occuring value
mode_func = lambda x: pd.Series.mode(x).iat[0]


def convert_deprecated_metapeptable_naming(input_df: pd.DataFrame) -> pd.DataFrame:
    """Attempt to fix old metapeptable files where namings of columns have
    changed in the meantime.

    Args:
        input_df (pd.DataFrame): Input MetaPepView dataset.

    Returns:
        pd.DataFrame: Converted dataset.
    """
    out_df = input_df.rename(columns={
        "Superkingdom Id": "Domain Id", 
        "Superkingdom Name": "Domain Name"}).copy()
    return out_df


def get_gene_options(ko_map_df: pd.DataFrame) -> List[Dict[str, str]]:
    """Fetch symbol data from KEGG Orthology mapping dataset and return

    Args:
        ko_map_df (pd.DataFrame): _description_

    Returns:
        Tuple[List[Dict[str, str]], str]: _description_
    """
    options = [{'label': val['symbol'], 
                'value': str(idx) + "|" + val['symbol']} \
        for idx, val in ko_map_df.iterrows()]
    return options


def filter_taxonomy_clade(dataset: pd.DataFrame,
                          root_taxonomy: str,
                          root_rank: RankType,
                          taxonomy_format: Literal['Id', 'Name']) -> pd.DataFrame:
    """Filter all peptide records from the dataframe whose taxonomy
    does not belong to the specified root taxonomy.

    Args:
        dataset (pd.DataFrame): MetaPepView format dataset.
        root_taxonomy (str): Taxonomy name as root (or common ancestor)
        root_rank (RankType): Taxonomy rank of root taxonomy.
        taxonomy_format (Literal['Id', 'Name']): Filter based on root taxonomy
            id or name.

    Returns:
        pd.DataFrame: Filtered dataset
    """
    rank_suffix = ' Id' if taxonomy_format == 'Id' else ' Name'
    return dataset[dataset[root_rank + rank_suffix] == root_taxonomy]



def match_db_search_psm(spectral_df: pd.DataFrame,
                        db_search_psm: pd.DataFrame) -> pd.DataFrame:
    """Add column to spectral data db that shows if scan
    has peptide identification.

    Args:
        spectral_df (pd.DataFrame): spectra dataset with scan number.
        db_search_psm (pd.DataFrame): Dataframe from MetaPepDbSearch object
    """
    if 'scan number' not in spectral_df.columns:
        raise ValueError("spectral data does not have scan number.")
    
    ident_scans = set(db_search_psm['Scan'].to_list())
    spectral_df['db search identified'] = spectral_df['scan number'].isin(ident_scans)
    return spectral_df


def match_de_novo(spectral_df: pd.DataFrame,
                   de_novo: pd.DataFrame) -> pd.DataFrame:
    """Add column to spectral data db that shows if scan
    has peptide identification.

    Args:
        spectral_df (pd.DataFrame): spectra dataset with scan number.
        de_novo (pd.DataFrame): Dataframe from MetaPepDeNovo object.
    """
    if 'scan number' not in spectral_df.columns:
        raise ValueError("spectral data does not have scan number.")
    
    ident_scans = set(de_novo['Scan'].to_list())
    spectral_df['de novo identified'] = spectral_df['scan number'].isin(ident_scans)
    return spectral_df


def filter_denovo_only(de_novo_df: pd.DataFrame,
                       db_search_df: pd.DataFrame) -> pd.DataFrame:
    de_novo_df["Sequence"] = de_novo_df['Peptide'].apply(wrangle_peptides)
    db_search_df["Sequence"] = db_search_df['Peptide'].apply(wrangle_peptides)
    return de_novo_df[~de_novo_df["Sequence"].isin(db_search_df["Sequence"])]


def fetch_sort_column(df: pd.DataFrame,
                      df_col: str) -> pd.Series:
    """Extract column out of dataframe and sort descending.
    sorted column is returned as series object. Index is reset

    Args:
        df (pd.DataFrame): Input dataset
        df_col (str): Column name to extract from dataset.

    Returns:
        pd.Series: Column values sorted in descending order.
    """
    score_series = df[df_col].sort_values(ascending=False)\
        .reset_index(drop=True)
    
    return score_series

def value_to_nan(series: pd.Series,
                 value: Any) -> pd.Series:
    """Convert all elements of specified value inside series to nan.

    Args:
        series (pd.Series): Pandas series object.
        value (Any): Value to convert.

    Returns:
        pd.Series: Coverted series object.
    """
    series[series == value] = np.nan
    return series

def filter_crap(df: pd.DataFrame,
                pept_col: str,
                crap_peptides: pd.Series) -> pd.DataFrame:
    """Parse peptides from the input dataset and filter all peptides that match
    to the cRAP dataset. 
    
    The cRAP (common Repository of Adventitious Proteins)
    dataset contains a collection of proteins that are often encountered in the
    environment in which protein samples are prepared. As a result, these
    proteins will often be present in (meta)proteomics data despite not being
    part of the sample.

    Args:
        df (pd.DataFrame): Input peptide dataset.
        pept_col (str): Column name to parse peptide sequences.
        crap_peptides (pd.Series): peptides present in cRAP dataset.

    Returns:
        pd.DataFrame: Peptide dataset filtered from cRAP peptides
    """
    # filter peptides based on presence crap peptides, equalize Leu and Iso
    li_pept_col = df[pept_col].transform(
        lambda x: wrangle_peptides(x, ptm_filter=False)
    )
    crap_peptides = crap_peptides.transform(
        lambda x: wrangle_peptides(x, ptm_filter=False)
    )
    
    return df[~li_pept_col.isin(crap_peptides.values)]


def substitute_lineage_with_global_lineage(peptide_df: pd.DataFrame,
                                           supplement_prot_db_annotation: bool = True) -> pd.DataFrame:
    """Add global annotation lineage data to regular taxonomy lineage columns
    for samples that have no taxonomy db classification.

    Args:
        peptide_df (pd.DataFrame): peptide dataset from MetaPepTable.
        supplement_prot_db_annotation (bool, Optional): Add global taxonomy 
            annotation data to samples that already have taxonomy annotations 
            through protein matching. If false, only samples without any local
            taxonomy annotation are processed. If true, global annotation data
            is added to lineage data for peptides that do not have any taxonomy
            classification. Defaults to False

    Returns:
        pd.DataFrame: peptide dataset with lineage substituted
    """
    
    # lineage columns
    glob_tax_fields = GlobalConstants.metapep_table_global_taxonomy_lineage
    mg_tax_fields = GlobalConstants.metapep_table_taxonomy_lineage
    
    # if true, supplement any peptide without local annotation with global annotation
    if supplement_prot_db_annotation is True:
        no_tax_annot_idx = peptide_df[peptide_df["Taxonomy Id"].isnull()].index
    # otherwise, only substitute samples that have no taxonomy annotation
    else:
        no_tax_annot_idx = peptide_df[peptide_df["Taxonomy Annotation"] == False].index

    peptide_df.loc[no_tax_annot_idx, 
                    mg_tax_fields] = peptide_df.loc[no_tax_annot_idx, 
                                                    glob_tax_fields].values
    
    # add de novo matches to psm count only if no db search was observed
    peptide_df.loc[no_tax_annot_idx, 
                   'PSM Count'] = peptide_df.loc[no_tax_annot_idx, 
                                                 'PSM Count'].fillna(peptide_df['De Novo Match Count'])   
    
    return peptide_df


def reshape_taxonomy_df_to_denovo(input_df: pd.DataFrame,
                                  global_annot_de_novo_only: bool) -> pd.DataFrame:
    """Wrangle sample taxonomy classification dataset into a dataset that
    compares metagenome taxonomy classification against global classification.

    Args:
        input_df (pd.DataFrame): MetaPepTable dataframe from one sample.
        global_annot_de_novo_only (bool): Only classify de novo identified data
            for global classification.

    Returns:
        pd.DataFrame: Wrangled dataframe
    """
    rank_list = GlobalConstants.standard_lineage_ranks
    global_suffix = GlobalConstants.global_annot_suffix
    
    lin_ids = [i + " Id" for i in GlobalConstants.standard_lineage_ranks]
    lin_names = [i + " Name" for i in GlobalConstants.standard_lineage_ranks]
    global_lin_ids = [i + " Id" + global_suffix for i in rank_list]
    global_lin_names = [i + " Name" + global_suffix for i in rank_list]
    
    mg_annot_cols = lin_ids + lin_names + ['PSM Count', 'Area']
    metagenome_annot = input_df[mg_annot_cols]
    global_annot = input_df[global_lin_ids + \
                            global_lin_names + \
                            ['PSM Count', 'Area', 'De Novo Area', 'De Novo Match Count']]
    
    # if db search matches included in global annotation profile, db search
    # results get priority in the output over de novo results
    if global_annot_de_novo_only is False:
        area_nan = global_annot['Area'].isna()
        psm_nan = global_annot['PSM Count'].isna()
        
        # fill all de novo area and match count with db search area (and match count) if not nan
        global_annot.loc[~area_nan, 'De Novo Area'] = global_annot.loc[~area_nan, 'Area']
        global_annot.loc[~psm_nan, 'De Novo Match Count'] = global_annot.loc[~psm_nan, 'PSM Count']
    
    # move de novo columns to same name as psm columns
    global_annot['Area'] = global_annot['De Novo Area'] 
    global_annot['PSM Count'] = global_annot['De Novo Match Count']
    global_annot = global_annot.drop(columns=['De Novo Area', 'De Novo Match Count'])
    
    # align column names and concatenate global annotation and metagenome annotation datasets
    renamer = {i: j for i, j in zip(global_lin_ids + global_lin_names,
                                    lin_ids + lin_names)}
    global_annot = global_annot.rename(columns=renamer)
    metagenome_annot['Sample Name'] = 'DB search taxonomy'
    global_annot['Sample Name'] = 'De novo taxonomy (Unipept)'
    
    return pd.concat([metagenome_annot, global_annot])


def peptide_allocation_across_lineage(peptide_df: pd.DataFrame,
                                      lineage: List[str | float],
                                      quant_col: str) -> Tuple[Tuple[str, int | float],
                                                               List[Tuple[str,
                                                                          Tuple[int | float]]]]:
    """Given a taxonomy lineage, compute allocation of peptide abundance across
    lineage.

    Args:
        peptide_df (pd.DataFrame): MetaPepTable peptide dataset
        lineage (List[str  |  float]): Taxonomy lineage.
        quant_col (str): Column used for abundance quantification

    Returns:
        Tuple[Tuple[str, int | float],
              List[Tuple[str,
                   Tuple[int | float]]]]: Set of allocation counts across lineage:
            First value contains tuple of 'rank:tax_name' ids with total quantification
            for each taxa. Second value contains Tuple of *valid* 'rank:tax_name' id's
            with for each taxa a tuple showing allocation of peptides from the lower
            rank clade: {annotation limit (LCA), Branched to other taxa, 
            allocated to current taxa}.
    """
    # fetch name and column data for lineage
    rank_letters = GlobalConstants.lineage_ranks_short    
    rank_names = GlobalConstants.standard_lineage_ranks
    lineage_cols = [rank + ' Name' for rank in rank_names]
    
    # Make index list for lineage, drop clades used for lineage gap filling
    valid_ranks = []
    for i in range(len(lineage)):
        if lineage[i] != lineage[i] or lineage[i] == "-":
            continue
        elif i == len(lineage) - 1 or lineage[i] != lineage[i + 1]:
            valid_ranks.append(i)
        else:
            # if not valid, create gap in lineage
            lineage[i] = "-"
            
    
    # for lineage, count number of matches at every rank
    lineage_counts = [
        (
            "{}: {}".format(rank_letters[x], lineage[x]), 
            np.nan if x not in valid_ranks else \
                peptide_df.groupby(by=lineage_cols[x])[quant_col]\
                    .agg('sum')\
                    .loc[lineage[x]]
        )
        for x in range(len(lineage))
    ]
    
    
    # initialize empty lineage dropoff
    lineage_dropoff = [
        (
            "{}: {}".format(rank_letters[x + 1], 
                            lineage[x + 1]), 
            (np.nan, np.nan, np.nan, np.nan, np.nan)
        )
        for x in range(len(lineage) - 1)
    ]
    
    # Fill lineage_dropoff with dropoff values between valid ranks
    valid_dropoffs = [(valid_ranks[i], 
                       valid_ranks[i + 1]) for i in range(len(valid_ranks) - 1)]
    
    for low_rank_idx, high_rank_idx in valid_dropoffs:
        lineage_dropoff[high_rank_idx - 1] = (
            "{}: {}".format(rank_letters[high_rank_idx], 
                            lineage[high_rank_idx]), 
            compute_taxonomy_dropoff(
                peptide_df,
                quant_col,
                lineage_cols[low_rank_idx],
                lineage_cols[high_rank_idx],
                lineage[low_rank_idx],
                lineage[high_rank_idx],
                rank_letters[high_rank_idx]
            )
        )

    # dropoff of lowest rank is just counting empty values
    lineage_dropoff = [
        (
            "{}: {}".format(rank_letters[valid_ranks[0]], 
                            lineage[valid_ranks[0]]),
            compute_taxonomy_dropoff(
                peptide_df,
                quant_col,
                np.nan,
                lineage_cols[valid_ranks[0]],
                np.nan,
                lineage[valid_ranks[0]],
                rank_letters[0]
            )
        )
    ] + lineage_dropoff
    
    return (lineage_counts, lineage_dropoff)


def compute_taxonomy_dropoff(
    peptide_df: pd.DataFrame,
    quant_col: str,
    rank_lower: str | float, 
    rank_upper: str,
    rank_lower_name: str | float, 
    rank_upper_name: str | float,
    name_prefix: str | None = None) -> Tuple[float, float, float, float]:
    """Get the drop in annotation counts when comparing peptide quantification
    from a single lower rank clade to the quantification of a higher
    rank clade of interest. Dropoff may be due to either branching of peptide
    annotations across other higher rank clades, or due to taxonomic resolution
    limit (LCA at lower rank).

    Args:
        peptide_df (pd.DataFrame): MetaPepTable dataset (assumes single sample).
        quant_col (str): Column to use for quantification.
        rank_upper (str): Upper rank name (e.g. "Phylum").
        rank_lower_name (str): Name of lower rank clade to quantify.
        rank_upper_name (str): Name of upper rank clade of interest (only used
            to check if upper rank is still defined. If not, no dropoff is 
            calculated). 
        name_prefix (str | None): Prefix to taxonomy name used in branching taxa.
            Can be used to add rank information.

    Returns:
        Tuple[float, float]: Dropoff values. 
            {annotation resolution limit dropoff,
             clade branching dropoff,
             valid abundance upper clade,
             total abundance lower clade}
    """
    # Process upper rank dropoff over whole peptide dataset if no lower rank given (root)
    if rank_lower != rank_lower or rank_lower == "-":
        
        # compute dropoff from taxonomy resolution loss
        total_count = peptide_df[quant_col].sum()
        annot_dropoff = peptide_df[
            peptide_df[rank_upper].isna()][quant_col].sum()
        
        if rank_upper_name != rank_upper_name\
            or rank_upper_name == "-":
            valid_counts = np.nan
            branch_dropoff = np.nan
        else:
            # Count groups of upper taxonomy rank
            valid_counts = peptide_df[
                peptide_df[rank_upper] == rank_upper_name][quant_col].sum()
            branch_dropoff = total_count - valid_counts - annot_dropoff
        
        # generate dictionary of branching counts: {'tax_1': 3, 'tax_2': 5, ...}
        branch_counts = peptide_df[
            peptide_df[rank_upper] != rank_upper_name
            ][[rank_upper, quant_col]]\
            .dropna()\
            .groupby(rank_upper)[quant_col]\
            .sum()
            
        if name_prefix is not None:
            branch_counts.index = branch_counts.index.map(
                lambda x: name_prefix + ": " + x
            )
        
        return (annot_dropoff, 
                branch_dropoff, 
                valid_counts, 
                total_count,
                branch_counts.to_dict())
    
    # Dropoff is only of interest for valid lower rank taxa names
    # Also, check if upper and lower ranks are different
    if rank_lower_name != rank_lower_name or rank_lower_name == "-" or\
        rank_upper == rank_lower:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    else:
        # set initial dropoff and quantification values to 0
        annot_dropoff, valid_counts = 0.0, 0.0
        
        # fill missing values to allow indexing, then count groups
        peptide_df[rank_upper] = peptide_df[rank_upper].fillna(-1.0)
        rank_counts = peptide_df\
            .groupby(by=[rank_lower, rank_upper])[quant_col]\
            .agg('sum')
        
        # get total count of the lower rank clade. If it is not present, we 
        # cannot compute any dropoff
        if rank_lower_name in rank_counts.index.get_level_values(rank_lower):
            lower_rank_total = rank_counts\
                .groupby(level=rank_lower)\
                .sum()\
                .loc[rank_lower_name]
        else:
            # undo peptide df modification
            peptide_df[rank_upper] = value_to_nan(peptide_df[rank_upper], -1.0)
            return (np.nan, np.nan, np.nan, np.nan, np.nan)
            
        # calculate annotation dropoff if any observed
        if (rank_lower_name, -1.0) in rank_counts.index:
            annot_dropoff = rank_counts.loc[(rank_lower_name, -1.0)]
        
        # valid counts is only meaningful if an upper rank name is specified
        if rank_upper_name != rank_upper_name or rank_upper_name == "-":
            valid_counts = np.nan
        elif (rank_lower_name, rank_upper_name) in rank_counts.index:
            valid_counts = rank_counts.loc[(rank_lower_name, rank_upper_name)]
            
        # generate dictionary of branching counts: {'tax_1': 3, 'tax_2': 5, ...}
        upper_rank_group = rank_counts.loc[rank_lower_name]
        branch_counts = upper_rank_group[
            ~upper_rank_group.index.isin([rank_upper_name, -1.0])       # drop nan values
        ]
        
        if name_prefix is not None:
            branch_counts.index = branch_counts.index.map(
                lambda x: name_prefix + ": " + x)
        
        # calculate branching dropoff
        branch_dropoff = lower_rank_total - valid_counts - annot_dropoff
        
        # undo peptide df modification
        peptide_df[rank_upper] = value_to_nan(peptide_df[rank_upper], -1.0)
        
        return (annot_dropoff, 
                branch_dropoff, 
                valid_counts, 
                lower_rank_total,
                branch_counts.to_dict())
        

def compute_lineage_cumulative_annotation_drop(
    peptide_allocation: List[Tuple[str, Tuple[int | float]]],
    root_rank: RankType | Literal["Root"]) -> float:
    """Combine the drop of peptide annotations as a result of LCA limits across
    ranks into a global observed annotation drop for a given lineage.
    
    This approach differs from calculating the loss of annotations from a
    specified root clade down to a higher rank level as that would count
    peptide losses across all subclades. This method iteratively takes the
    peptide annotation loss from taxa x at its corresponding rank to the next
    valid rank, across the complete lineage. This means that any divergent clade
    group is ignored in the next annotation loss computation.
    
    Args:
        lineage_counts (Tuple[str, int  |  float]): _description_
        peptide_allocation (List[Tuple[str, Tuple[int  |  float]]]): _description_
        root_rank (RankType | Literal[&quot;Root&quot;]): _description_

    Returns:
        float: _description_
    """
    # fetch allocation categories into separate arrays
    annotation_dropoff = np.array([x[0] for i, x in peptide_allocation])
    sum_array = np.array([x[3] for i, x in peptide_allocation])

    # calculate percentage drop for each rank
    annotation_dropoff_frac = annotation_dropoff / sum_array
    
    if root_rank == "Root":
        start_idx = 0
    # species is the resolution limit, cannot drop beyond that
    elif root_rank == "Species":
        return np.nan
    else:
        start_idx = GlobalConstants.standard_lineage_ranks.index(root_rank) + 1
    
    # remove values from array that fall outside the root clade rank
    annotation_dropoff_frac = annotation_dropoff_frac[start_idx:]
    # remove potential missing values
    annotation_dropoff_frac = annotation_dropoff_frac[~np.isnan(annotation_dropoff_frac)]

    if len(annotation_dropoff_frac) == 0:
        return np.nan
    else:
        return (1 - np.multiply.reduce(1 - annotation_dropoff_frac)) * 100
    
    
def compute_global_cumulative_annotation_drop(
    peptide_df: pd.DataFrame,
    root_rank: RankType | Literal["Root"],
    root_clade: str | float,
    limit_rank: RankType,
    quant_col: str) -> float:
    """Compute the drop of peptide annotations as a result of LCA limits between
    a specified root clade up to the desired limit rank
    
    In contrast to `compute_lineage_cumulative_annotation_drop`, this method
    computes the fraction of peptide annotations lost directly from a root clade
    up to a specified higher rank. This serves as a global average annotation
    drop from the complete clade, as opposed to the drop from a specific lineage.
    
    Args:
        peptide_df (pd.DataFrame): MetaPepTable dataset.
        root_rank (RankType | Literal["Root"]): Rank start peptide annotation 
            drop computation.
        root_clade (str | float): Taxonomy clade to set as root clade for
            annotation drop computation.
        limit_rank (RankType): Taxonomy rank to stop peptide annotation drop.
        quant_col (str): Column used for peptide abundance quantification

    Returns:
        float: Annotation drop observed in percent.
    """
    
    if root_rank == "Root":
        root_rank = np.nan
    else:
        root_rank = root_rank + " Name"
    limit_rank = limit_rank + " Name"
    
    # get the first value for the taxonomy dropoff
    dropoff = compute_taxonomy_dropoff(
        peptide_df,
        quant_col,
        root_rank,
        limit_rank,
        root_clade,
        np.nan)
    
    total_loss, total_quant = dropoff[0], dropoff[3]
    
    return ((total_loss / total_quant)) * 100
