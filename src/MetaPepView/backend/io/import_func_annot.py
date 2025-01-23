import pandas as pd
from ..utils import *


def validate_eggnog_file(upload_contents: str,
                         archive_format: str | None = None) -> Tuple[bool, None | str]:

    # required columns for file
    emapper_cols = ["#query", "evalue", "score", "eggNOG_OGs", "COG_category",
                    "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
                    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction",
                    "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "PFAMs"]

    # import content and extract data if delivered in archive
    file_buffer = memory_to_stringio(upload_contents, archive_format)

    # read data into dataframe format, only first 100 rows
    try:
        df = pd.read_csv(file_buffer,
                        sep='\t', header=4, nrows=100)
    except:
        err_msg = "Failed to read input data, is it '*.emapper.annotation' file?"
        return (False, err_msg)
        
    # check if all expected columns are present
    for col in emapper_cols:
        if col not in df.columns:
            err_msg = f"Column '{col}' not in input data, correct file imported? See 'Data Format' above."
            return (False, err_msg)
    
    return (True, None)
    