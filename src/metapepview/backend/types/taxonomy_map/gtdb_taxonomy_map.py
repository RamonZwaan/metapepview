from typing import IO, Self, Any, Sequence
import pandas as pd
import numpy as np

from metapepview.backend.types.taxonomy_map.base_class import AccessionTaxaMap
from metapepview.backend.types import GtdbTaxonomy, GtdbGenomeToNcbi
from metapepview.backend.utils import wrangle_peptides, regex_over_column


class AccessionTaxaMapGtdb(AccessionTaxaMap):
    
    TAXONOMY_FORMAT = "GTDB"
    ROOT_TAXONOMY_ID = None # gtdb consists of ref genomes with lineage, no root tax id is present
    
    def __init__(self,
                 accession_to_taxa: pd.DataFrame):
        super().__init__(accession_to_taxa)


    @classmethod
    def from_string_buffer(cls,
                           str_file_obj: IO[str],
                           acc_col: int=0,
                           tax_col: int=1,
                           acc_regex: str | None=None,
                           delimiter: str | None=",",
                           drop_duplicates: bool = True,
                           tax_name_to_id: bool = False,
                           taxonomy_obj: GtdbTaxonomy | None = None,
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
            taxonomy_obj (GtdbTaxonomy | None, optional): GtdbTaxonomy
                object for LCA processing of redundant protein accessions.
                Defaults to None.
            wrangle_peptide_accessions (bool, optional): If accession is peptide
                sequence column, perform removal of non-amino acid elements and
                equate Leucin and Isoleucin. Defaults to False.

        Returns:
            Self: AccessionTaxaMapGtdb instance.
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
        
        
        # if genome id's given, convert to species id
        if tax_name_to_id is False and taxonomy_obj is not None:
            prot_df.loc[:, "taxonomy_id"] = prot_df["taxonomy_id"].apply(
                lambda x:taxonomy_obj.genome_id_to_species_id(x) if\
                    taxonomy_obj.genome_id_in_dataset(x) else x)
        
        # if tax names in column, convert to id
        if tax_name_to_id is True and taxonomy_obj is not None:
            print("Taxonomy name to id conversion\n\nFailed conversions:")  
            prot_df.loc[:, "taxonomy_id"] = prot_df["taxonomy_id"]\
                .apply(taxonomy_obj.name_to_id, print_fails=True)
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

     
    
    def accession_list_to_lca(self,
                              accessions: Sequence[Any] | pd.Series | float | None,
                              taxonomy_db: GtdbTaxonomy) -> int | str | float:
        """Return the last common ancestor from taxonomy id's related
        to the list of accession id's supplied.

        Args:
            accessions (List[Any] | Tuple[Any] | pd.Series | float | None): list of accession id's
            taxonomy_db (GtdbTaxonomy): Taxonomy database object
        Returns:
            int | str | float: Last common ancestor
        """
        # if nan or None in accessions, return that value
        if isinstance(accessions, float | int) or accessions is None:
            return np.nan
        
        # convert series object to list
        if isinstance(accessions, pd.Series):
            accessions = accessions.to_list()

        tax_ids = [self.accession_to_taxonomy(i, taxonomy_db) for i in accessions]
        
        return taxonomy_db.taxa_to_lca(tax_ids)
    
    
    def accession_to_taxonomy(self,
                              accession: str,
                              taxonomy_db: GtdbTaxonomy) -> str | float:
        """Return the taxonomy annotation for a given protein accession.

        Args:
            accession (str): Protein accession.

        Returns:
            str | float: Taxonomy annotation from protein accession.
        """
        # return nan if accession does not exist in the database
        return self.accession_tax_dict.get(accession, np.nan)
   
   
