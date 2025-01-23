from typing import IO, Self, Any, Sequence
import pandas as pd
import numpy as np

from .base_class import AccessionTaxaMap
from backend.types import GtdbTaxonomy


class AccessionTaxaMapGtdb(AccessionTaxaMap):
    
    TAXONOMY_FORMAT = "GTDB"
    ROOT_TAXONOMY_ID = None # gtdb consists of ref genomes with lineage, no root tax id is present
    
    def __init__(self,
                 accession_to_taxa: pd.DataFrame):
        super().__init__(accession_to_taxa)
    
    
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
        gtdb_genome = self.accession_tax_dict.get(accession, np.nan)
        return taxonomy_db.genome_id_to_species_id(gtdb_genome)
    
    
    def accession_to_genome(self, accession: str) -> str | float:
        """Return gtdb genome id from accession.

        Args:
            accession (str): Protein accession.

        Returns:
            str | float: Gtdb genome annotation from protein accession.
        """
        return self.accession_tax_dict.get(accession, np.nan)
    