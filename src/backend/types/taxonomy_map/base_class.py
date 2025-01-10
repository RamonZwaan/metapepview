from pathlib import Path
from typing import IO, Self, Sequence

import numpy as np
import pandas as pd

from backend.utils import *
from backend.types import TaxonomyDatabase


class AccessionTaxaMap:

    def __init__(self,
                 accession_to_taxa: pd.DataFrame):
        # accession_to_taxa cols: accession, taxonomy id
        self.accession_set = set(accession_to_taxa.index.to_list())
        self.accession_tax_dict = accession_to_taxa.set_index('accession').to_dict()['taxonomy_id']


    @classmethod
    def from_string_buffer(cls,
                           str_file_obj: IO[str],
                           acc_col: int=0,
                           tax_col: int=1,
                           acc_regex: str | None=None,
                           delimiter: str | None=", ") -> Self:
        
        # it is assumed that the dataset starts with accession, with second column tax id
        prot_df = pd.read_csv(str_file_obj,
                              usecols=[acc_col, tax_col], 
                              names=["accession", "taxonomy_id"],
                              sep=delimiter)

        # apply pattern on accession if given        
        if acc_regex is not None:
            repl = lambda m: m.group(0) if m is not None else ""
            prot_df["accession"] = prot_df["accession"].str.replace(acc_regex, repl, regex=True)
        
        #TODO: extend functionality towards lca or exception throw in case of duplicates
        # drop duplicate items if identical accessions found
        prot_df = prot_df.drop_duplicates(subset="accession")

        return cls(prot_df)

    @staticmethod
    def validate_input(str_file_obj: IO[str],
                       acc_col=0,
                       tax_col=1,
                       delimiter=", ",
                       custom_validators: Sequence[Callable[[pd.DataFrame], Tuple[bool, str | None]]] | None = None) -> Tuple[bool, str | None]:
        """Function definition to check if input data is valid for the class
        object.

        Args:
            str_file_obj (IO[str]): Input file buffer
            acc_col (int, optional): Protein accession column id. Defaults to 0.
            tax_col (int, optional): Taxonomy id column id. Defaults to 1.
            delimiter (str, optional): Delimiter character of input file.
                Defaults to ", ".
        """
        # get first 100 rows
        try:
            prot_df = pd.read_csv(str_file_obj,
                                  sep=delimiter,
                                  nrows=100)
        except:
            err_msg = "Failed to read input, should be delimited text file like '.csv'"
            return (False, err_msg)
        
        # check that indices are different and inside table range
        max_idx = len(prot_df.columns) - 1

        # no delimiter found, incorrect delimiter
        if max_idx == 0:
            err_msg = "Only one column in input data, is delimiter correct?"
            return (False, err_msg)
        
        # index found outside range           
        if not all(i <= max_idx for i in (acc_col, tax_col)):
            err_msg = f"Supplied column indices exeed range input data, max index: {max_idx}"
            return (False, err_msg)
        
        if acc_col == tax_col:
            err_msg = "Accession and Taxonomy column index cannot be equal."
            return (False, err_msg)
        
        if custom_validators is not None:
            for validator in custom_validators:
                success, err_msg = validator(prot_df)
                if success is False:
                    return (False, err_msg)
                
        return (True, None)


    def accession_list_to_lca(self,
                              accessions: Sequence[Any] | pd.Series | float | None,
                              taxonomy_db: TaxonomyDatabase) -> int | str | float:
        """Return the last common ancestor from taxonomy id's related
        to the list of accession id's supplied.

        Args:
            accessions (List[Any] | Tuple[Any] | pd.Series | float | None): list of accession id's
            taxonomy_db (TaxonomyDatabase): Taxonomy database object
        Returns:
            int | str | float: Last common ancestor
        """
        # if nan or None in accessions, return that value
        if isinstance(accessions, float | int) or accessions is None:
            return np.nan
        
        # convert series object to list
        if isinstance(accessions, pd.Series):
            accessions = accessions.to_list()

        tax_ids = [self.accession_to_taxonomy(i) for i in accessions]
        
        return taxonomy_db.taxa_to_lca(tax_ids)
    
    
    def accession_to_taxonomy(self,
                              accession: str) -> int | str | float:
        """Return the taxonomy annotation for a given protein accession.

        Args:
            accession (str): Protein accession.

        Returns:
            int | str | float: Taxonomy annotation from protein accession.
        """
        # return nan if accession does not exist in the database
        return self.accession_tax_dict.get(accession, np.nan)