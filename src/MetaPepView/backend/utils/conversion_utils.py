"""
This module contains functions that convert a variety of input formats into
custom Metapepview data formats.
"""

import pandas as pd
from datetime import datetime
import time
from typing import IO, Type, Dict, TypeVar, Callable, List


def custom_groupby(
    df: pd.DataFrame,
    groupby_col: str,
    custom_aggs: Dict[str, str | Callable] | None = None,
    match_idxmax: Dict[str, str | List[str]]| None = None,
    include_group_size: bool = True,
    group_size_name: str = 'group_size'
    ) -> pd.DataFrame:
    """Group a DataFrame by a specified column with custom aggregations and group size.

    Args:
        df (pd.DataFrame): Input DataFrame
        groupby_col (str): Column to group by
        custom_aggs (Dict[str, str  |  Callable] | None, optional): Dict of
            column names and custom aggregation functions. Defaults to None.
        match_idxmax (Dict[str, str | List[str]] | None, optional): Dict of
            columns whose value selection is bound to the maximum value of a
            different column. Key is the column whose idxmax is taken. Value
            the selection column('s). Defaults to None.
        include_group_size (bool, optional): Specify if group size should be
            added as separate column. Defaults to True.
        group_size_name (str, optional): Name of column storing group size, if
            included. Defaults to 'group_size'.
    Returns:
        pd.DataFrame: Grouped DataFrame
    """
    # Define default aggregation (first value for all columns except groupby_col)
    agg_dict: Dict[str, str | Callable]
    agg_dict = {col: 'first' for col in df.columns if col != groupby_col}
    
    # Update with custom aggregations
    if custom_aggs is not None:
        agg_dict.update(custom_aggs)
    
    # Perform groupby with aggregations
    groupby_df = df.groupby(groupby_col)
    result = groupby_df.agg(agg_dict)
    
    # Include group size if requested
    if include_group_size is True:
        result[group_size_name] = groupby_df.size()
    
    # Handle max value columns
    if match_idxmax is not None:
        for max_col, cols in match_idxmax.items():
            df_max_col = df[[groupby_col, max_col]]
            df_max_col.loc[:, max_col] = df_max_col[max_col].fillna(0)
            max_indices = df_max_col.groupby(groupby_col)[max_col].idxmax()
            result[cols] = df.loc[max_indices, cols].values # type: ignore
    
    return result.reset_index()


def quantity_to_fold_change(dataset: pd.DataFrame,
                            group_name: str,
                            group_col: str,
                            rescale_column: str,
                            value_column: str) -> pd.DataFrame:
    """Rescale values from a specified column to the values from a
    user specified group. This is used to gain insight into fold change data
    between groups of interest.
    
    The group_name defines the group whose values the dataset will rescale to.
    This way, all values from all groups represent a fold change to the group_name.
    
    The rescale_column defines the categories that are compared. Values
    belonging to the same category are rescaled to the category value of the
    group_name.
    
    Finally, the value_column represents the column whose values are rescaled.
    
    Example:
    >>> df = pd.DataFrame({'groups': ['a', 'a', 'b', 'b'], 'cats': ['x', 'y', 'x', 'y'], 'vals': [3, 2, 2, 4]})
      groups   cats   vals
    0      a      x      3
    1      a      y      2
    2      b      x      2
    3      b      y      4

    >>> quantity_to_fold_change(df, 'a', 'groups', 'cats', 'vals')
      groups   cats   vals
    0      a      x      1
    1      a      y      1
    2      b      x    2/3
    3      b      y      2  

    Args:
        dataset (pd.DataFrame): Input dataset.
        group_name (str): Name of group to set rescale factor. In this group,
            the values of the value_column are set to one (divided by itself).
        group_col (str): Column to search group name in.
        rescale_column (str): Categorical column to compare groups. Here, values
            in the value column will be rescaled by division of the group_name
            if the same category.
        value_column (str): Column of values to rescale.

    Returns:
        pd.DataFrame: Dataset with normalized values
    """
    # assert that value column is numeric
    if not pd.api.types.is_numeric_dtype(dataset[value_column]):
        raise TypeError("Non-numeric column supplied for rescaling.")
    
    # get psm data from sample that should be normalized
    sample_df = dataset[dataset[group_col] == group_name]
    
    # iterate through all protein types in the dataset
    for category in dataset[rescale_column].unique():
        # obtain category value from normalization sample
        norm_row = sample_df[sample_df[rescale_column] == category]
        
        # get indices of from dataset that belong to the category
        category_idx = dataset[dataset[rescale_column] == category].index
        
        # if sample has no value, set all values to 0
        if norm_row.shape[0] == 0:
            norm_val = 0.0
        # when multiple instances of category in group_name, sum values
        elif norm_row.shape[0] > 1:
            norm_val: float = norm_row[value_column].sum()
        else:
            norm_val: float = norm_row[value_column].item()
        
        # if category is missing from rescale group, set all to 0 as rescale is not possible
        if norm_val == 0.0:
            dataset.loc[category_idx, value_column] = 0.0
        else:
            # divide values from category by rescale value from group_name
            dataset.loc[category_idx, value_column] = dataset.loc[category_idx, value_column] / norm_val # type:ignore
    return dataset
