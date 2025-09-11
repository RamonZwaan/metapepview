"""This module contains functions to match de novo peptide sequences
to the unipept database.
"""

import pandas as pd
import numpy as np

from metapepview.backend.types import NcbiTaxonomy
from metapepview.backend.type_operations import add_lineage
from metapepview.backend.utils import request_unipept_pept_to_lca
from metapepview.constants import GlobalConstants


def global_taxonomic_annotation(pept_df: pd.DataFrame,
                                taxonomy_db: NcbiTaxonomy) -> pd.DataFrame:
    """
    Perform taxonomic annotation of  peptide sequences by matching against
    a global reference protein database. The peptide sequences are processed by
    unipept.

    Args:
        pept_df (pd.DataFrame): peptide dataset.

    Raises:
        RuntimeError: Unipept request unsuccessfull

    Returns:
        pd.DataFrame: Peptide dataset with global LCA data included
    """
    # get dataset of peptide sequences, lca taxonomy and lca rank
    lca_df = request_unipept_pept_to_lca(pept_df, 'Sequence')
    
    # fetch lineage
    lca_lin_df = add_lineage(lca_df, taxonomy_db, tax_col="Global LCA")
    
    # rename lineage to separate from db search lineage
    lin = GlobalConstants.standard_lineage_ranks
    suffix = GlobalConstants.global_annot_suffix
    lca_lin_df = lca_lin_df.rename(
        columns={i + " Id": i + f" Id{suffix}" for i in lin} |\
            {i + " Name": i + f" Name{suffix}" for i in lin})
    
    # filter away any taxonomy branch that has less than the unique peptide count threshold
    def values_to_counts(series: pd.Series) -> pd.Series:
        unique_count = series.value_counts()
        return series.apply(lambda x: unique_count.get(x, np.nan))
       
    lin_cols = GlobalConstants.metapep_table_global_taxonomy_lineage 
    lin_counts = lca_lin_df[lin_cols].apply(values_to_counts, axis=0)
    lca_lin_df.loc[:, lin_cols] = lca_lin_df[lin_cols].where(
        lin_counts >= GlobalConstants.global_annot_min_pept_tax)
    
    # join global search data to peptide dataframe and return
    return pd.merge(pept_df, lca_lin_df, on='Sequence', how='left')
