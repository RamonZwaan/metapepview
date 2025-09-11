from pathlib import Path
from typing import IO, Self, Sequence, TypeVar, Protocol

import numpy as np
import pandas as pd

from metapepview.backend.utils import *
from metapepview.backend.types import TaxonomyDatabase



class AccessionTaxaMapMethods(Protocol):
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
        ...
    
    
    def accession_to_taxonomy(self,
                              accession: str) -> int | str | float:
        """Return the taxonomy annotation for a given protein accession.

        Args:
            accession (str): Protein accession.

        Returns:
            int | str | float: Taxonomy annotation from protein accession.
        """
        ...


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
                           delimiter: str | None=",",
                           drop_duplicates: bool = True,
                           tax_name_to_id: bool = False,
                           taxonomy_obj: TaxonomyDatabase | None = None,
                           wrangle_peptide_accessions: bool = False) -> Self:
        """Import protein accession to taxonomy mapping data from string buffer
        generated from file.

        Args:
            str_file_obj (IO[str]): String buffer data
            acc_col (int, optional): Protein accession column index.
                Defaults to 0.
            tax_col (int, optional): Taxonomy column index. Defaults to 1.
            acc_regex (str | None, optional): Regex pattern to process accession
                column, storing only the pattern match portion. Defaults to None.
            delimiter (str | None, optional): Column delimiter. Defaults to ",".
            drop_duplicates (bool, optional): In case of duplicate protein
                accessions, keep only first occurrence. If false, duplicate
                accessions are merged by taking the last common ancestor for all
                taxa within the protein accession group. For this operation a
                TaxonomyDatabase object needs to be added. Defaults to True.
            tax_name_to_id (bool, optional): When elements in taxonomy column
                are names instead of id's, convert values to id. Defaults to
                False.
            taxonomy_obj (TaxonomyDatabase | None, optional): TaxonomyDatabase
                object for LCA processing of redundant protein accessions.
                Defaults to None.
            wrangle_peptide_accessions (bool, optional): If accession is peptide
                sequence column, perform removal of non-amino acid elements and
                equate Leucin and Isoleucin. Defaults to False.

        Returns:
            Self: AccessionTaxaMap instance.
        """
        
        
        # it is assumed that the dataset starts with accession, with second column tax id
        prot_df = pd.read_csv(str_file_obj,
                              usecols=[acc_col, tax_col], 
                              names=["accession", "taxonomy_id"],
                              sep=delimiter,
                              engine="python")
        
        # if accessions are peptides, wrangle sequences into consistent format
        if wrangle_peptide_accessions is True:
            prot_df.loc[:, "accession"] = prot_df["accession"].apply(wrangle_peptides)
        
        # if tax names in column, convert to id
        if tax_name_to_id is True and taxonomy_obj is not None:
            print("Taxonomy name to id conversion\n\nFailed conversions:")  
            prot_df.loc[:, "taxonomy_id"] = prot_df["taxonomy_id"].apply(
                taxonomy_obj.name_to_id, 
                print_fails=True)
        elif tax_name_to_id is True:
            raise ValueError("Taxonomy database required for name to id conversion...")
               
        # If regex provided, extract substring from each query using regex
        # if match, get first occurence of match
        if acc_regex != "" and acc_regex is not None:
            prot_df.loc[:, 'accession'] = regex_over_column(prot_df['accession'],
                                                            acc_regex,
                                                            no_match_to_nan=False)


        # manage duplicate protein names if present in dataset
        if not (prot_df['accession'].duplicated() == False).all():
            if drop_duplicates is True:
                prot_df = prot_df.drop_duplicates(subset="accession")
            elif drop_duplicates is False and\
                taxonomy_obj is not None:
                # full duplicates will be discarded anyway
                prot_df = prot_df.drop_duplicates(subset=["accession", "taxonomy_id"])

                # only compute LCA if duplicates still persist with diverging taxa
                dupl_prot = prot_df["accession"].duplicated(keep=False)
                if any(dupl_prot):
                    print("Aggregating protein duplicates...")
                    # create a copy of the protein db with only duplicates
                    dupl_prot_df = prot_df[dupl_prot].copy(deep=True)
                    prot_df = prot_df[~dupl_prot]
                    dupl_prot_agg = dupl_prot_df.groupby(by="accession")\
                                     .aggregate(taxonomy_obj.taxa_to_lca)\
                                     .reset_index()
                    
                    prot_df = pd.concat([prot_df, dupl_prot_agg])
            else:
                raise ValueError("Taxonomy database missing, add TaxonomyDatabase object or set drop_duplicates to True.")

        return cls(prot_df)


    @staticmethod
    def validate_input(str_file_obj: IO[str],
                       acc_col=0,
                       tax_col=1,
                       delimiter=", ",
                       custom_validators: Sequence[Callable[[pd.DataFrame], Tuple[bool, str | None]]] | None = None,
                       **kwargs) -> Tuple[bool, str | None]:
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
                                  engine='python',
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