"""This module defines separate functions that are used by the class objects,
but are not part of the class methods itself.
"""

from typing import Type, TypeVar, Dict, Any
import pandas as pd
import json

from metapepview.backend.utils.io_utils import compress_string


def to_json(data: pd.DataFrame,
            metadata_dict: Dict[str, Any]) -> str:
    """Write Metapep object to json format and store at file location.

    Args:
        data (pd.DataFrame): Data from object.
        metadata_dict (Dict[str, Any]): Dictionary of metadata to add to json.

    Returns:
        str: json formatted data as string.
    """
    # Serialize the DataFrame to json
    # df_json = data.to_json()
    # Serialize the DataFrame and compress
    df_csv_comp = compress_string(data.to_csv())

    # Serialize metadata to json
    metadata_json = json.dumps(metadata_dict)
    
    # combine both data into single json
    # combined_json = json.dumps({'dataframe': json.loads(df_json),
    #                             'metadata': json.loads(metadata_json)})
    combined_json = json.dumps({'dataframe': df_csv_comp,
                                'metadata': json.loads(metadata_json)})
    
    return combined_json
