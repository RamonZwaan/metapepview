from typing import IO, Tuple
import pandas as pd

from .base_class import AccessionTaxaMap


class AccessionTaxaMapNcbi(AccessionTaxaMap):
    
    TAXONOMY_FORMAT = "NCBI"
    ROOT_TAXONOMY_ID = 1
    
    def __init__(self,
                 accession_to_taxa: pd.DataFrame):
        super().__init__(accession_to_taxa)
        
    
    @staticmethod
    def validate_input(str_file_obj: IO[str],
                       acc_col=0,
                       tax_col=1,
                       delimiter=", "):
        # check correct format for ncbi taxonomy
        def check_numeric(prot_df: pd.DataFrame) -> Tuple[bool, str | None]:
            if not pd.api.types.is_numeric_dtype(prot_df.iloc[:, tax_col]):
                err_msg = "Invalid taxonomy format, expects ncbi id's (numeric)"
                return (False, err_msg)
            else:
                return (True, None)
            
        # perform regular validation
        success, err_msg = AccessionTaxaMap.validate_input(
            str_file_obj,
            acc_col,
            tax_col,
            delimiter,
            [check_numeric])
        
        return (success, err_msg)
