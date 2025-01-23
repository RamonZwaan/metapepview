"""This module contains helper functions related to pandas df processing,
Here, string processing is omitted as it is part of the string_utils
"""

from typing import List, Tuple, Dict
from pathlib import Path

import pandas as pd

from constants import *
from .string_utils import digest_proteins, wrangle_peptides


# return the modus of a pandas series, only the most occuring value
mode_func = lambda x: pd.Series.mode(x).iat[0]


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


def substitute_lineage_with_global_lineage(peptide_df: pd.DataFrame) -> pd.DataFrame:
    """Add global annotation lineage data to regular taxonomy lineage columns
    for samples that have no taxonomy db classification.

    Args:
        peptide_df (pd.DataFrame): peptide dataset from MetaPepTable

    Returns:
        pd.DataFrame: peptide dataset with lineage substituted
    """
    
    # lineage columns
    glob_tax_fields = GlobalConstants.metapep_table_global_taxonomy_lineage
    mg_tax_fields = GlobalConstants.metapep_table_taxonomy_lineage
    
    # only substitute samples that have no taxonomy annotation
    no_tax_annot_idx = peptide_df[peptide_df["Taxonomy Annotation"] == False].index
    peptide_df.loc[no_tax_annot_idx, 
                    mg_tax_fields] = peptide_df.loc[no_tax_annot_idx, 
                                                    glob_tax_fields].values
    
    # add de novo matches to psm count only if no db search was observed
    peptide_df['PSM Count'] = peptide_df['PSM Count'].fillna(peptide_df['De Novo Match Count'])   
    
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
    metagenome_annot['Sample Name'] = 'Metagenome annotation'
    global_annot['Sample Name'] = 'Global annotation'
    
    return pd.concat([metagenome_annot, global_annot])