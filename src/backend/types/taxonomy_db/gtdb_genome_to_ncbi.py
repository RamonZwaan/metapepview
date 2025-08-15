from pathlib import Path
from typing import Dict, Self

import pandas as pd
import numpy as np



class GtdbGenomeToNcbi:

    def __init__(self, genome_to_ncbi_dict: Dict[str, float]):
        self.genome_to_ncbi_dict = genome_to_ncbi_dict
    

    def genome_to_ncbi(self, genome_id: str) -> float:
        """Return NCBI taxonomy id for a given Genbank or Refseq genome id as
        present in the GTDB database.

        Args:
            genome_id (str): Genbank or Refseq genome id.

        Returns:
            float: NCBI taxonomy id. If no match, NaN is returned.
        """
        
        # remove genbank/refseq prefixes
        genome_id = genome_id.removeprefix("GB_").removeprefix("RS_")
        return self.genome_to_ncbi_dict.get(genome_id, np.nan)
     
     
    @classmethod
    def from_metadata_file(cls,
                           bac_metadata: str | Path,
                           ar_metadata: str | Path) -> Self:
        """Import GtdbGenomeToNcbi object from metadata files

        Args:
            bac_metadata (str | Path): Bacterial metadata.
            ar_metadata (str | Path): Archaeal metadata

        Returns:
            Self: Instance of GtdbGenomeToNcbi.
        """
        bac_metadata = Path(bac_metadata)
        ar_metadata = Path(ar_metadata)
        
        # check presence of files
        if not all([x.exists() for x in [bac_metadata, ar_metadata]]):
            raise FileNotFoundError("GTDB metadata files not found in path location...")
        
        # import bacteria and archaea into dataframe and concatenate them
        bac_df = pd.read_csv(bac_metadata,
                             sep="\t",
                             engine="python"
                             )
        arch_df = pd.read_csv(ar_metadata,
                              sep="\t",
                              engine="python")
        metadata_df = pd.concat([bac_df, arch_df], axis=0).reset_index(drop=True)

        # only the genome column and ncbi taxonomy id column are of interest
        metadata_df = metadata_df[["accession", "ncbi_taxid"]]
        metadata_df["ncbi_taxid"] = metadata_df["ncbi_taxid"].astype(float)
        
        # if present, remove genbank/refseq prefixes
        metadata_df["accession"] = metadata_df["accession"].apply(
            lambda x: x.removeprefix("GB_").removeprefix("RS_")
        )
        
        genome_to_ncbi_dict = dict(metadata_df.values)
        
        return cls(genome_to_ncbi_dict)
        
        
        
        