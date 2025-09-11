"""This module describes the object that manages gtdb taxonomies.

Author: Ramon van der Zwaan
Date: 03-04-2023
"""

from pathlib import Path
from typing import Dict, List, Tuple, overload, Literal, TypeVar, Sequence
from collections import defaultdict
from functools import lru_cache
import re

import pandas as pd
import numpy as np

from metapepview.backend.types.taxonomy_db.taxonomy_database import TaxonomyDatabase
from metapepview.constants import GlobalConstants


class GtdbTaxonomy(TaxonomyDatabase):
    
    FILE_NAMES = GlobalConstants.gtdb_taxonomy_files
    
    LINEAGE_COLUMNS = ("taxonomy_id", "lineage")
    
    # GTDB taxonomy lacks the Kingdom rank in the lineage
    RANK_LIST = [x.lower() for x in GlobalConstants.standard_lineage_ranks]
    
    RANK_ID_TO_NAME = {"d": "domain",
                       "p": "phylum",
                       "c": "class",
                       "o": "order",
                       "f": "family",
                       "g": "genus",
                       "s": "species"}
    
    ROOT_NAME = "Root"
    
    RANK_PREFIX = re.compile(r"[dpcofgs]__$")
    
    
    def __init__(self,
                 lineage_dict: Dict[str, List[str]],
                 child_dict: Dict[str, List[str]],
                 genome_species_map: Dict[str, str]):
        # map each taxa to lineage 
        self.lineage_dict = lineage_dict    
        # map each taxa to offspring taxa
        self.child_dict = child_dict
        
        # map genome id to species taxonomy
        self.genome_species_map = genome_species_map
    
    
    @classmethod
    def get_root_id(cls) -> str:
        """Return the global root id for the gtdb taxonomy database.
        Here, since there is no official name, the object returns 
        "Root".

        Returns:
            str: The global root for the gtdb taxonomy database.
        """
        return cls.ROOT_NAME
    
    
    @overload
    @classmethod
    def lineage_to_id(cls,
                      lineages: List[str | float],
                      root_on_empty: bool) -> str | float:
        ...

    @overload
    @classmethod
    def lineage_to_id(cls,
                      lineages: List[List[str | float]],
                      root_on_empty: bool) -> List[str | float]:
        ...
        
    @overload
    @classmethod
    def lineage_to_id(cls,
                      lineages: np.ndarray,
                      root_on_empty: bool) -> np.ndarray | str | float:
        ...
    
    @classmethod
    def lineage_to_id(cls,
                      lineages: List[str | float] | List[List[str | float]] | np.ndarray,
                      root_on_empty: bool = True) -> List[str | float] | np.ndarray | str | float:
        """Convert lineage vector or matrix of lineages to the highest valid taxonomy id for each lineage

        Args:
            lineages (List[str | float] | List[List[str | float]] | np.ndarray): Array or matrix of lineages.
            root_on_empty (bool, Optional): Specify if root taxonomy should be returned if an empty lineage
                is supplied. Else, nan will be returned. Defaults to True.
            root_value (Any, Optional): Root taxonomy id. Will be returned when empty lineage is supplied if
                specified from `root_on_empty`. Defaults to "Root".
                
        Returns:
            List[str | float] | np.ndarray | str | float: Taxonomy id for each lineage
        """
        # manage different formats of input data
        convert_list = False
        if isinstance(lineages, list):
            convert_list = True
            lineages = np.array(lineages)
        
        # specify output format (array if input = 2d, scalar if input = 1d)
        return_array = False
        if lineages.ndim == 2:
            return_array = True
        elif lineages.ndim == 1:
            lineages.reshape(1, -1)
        else:
            raise ValueError(f"Invalid lineage data, dimensionality of '{lineages.ndim}' not supported.")
        
        # for any lineage, take the last valid taxonomy value        
        tax_vector =  np.take_along_axis(lineages,
                                         (~(lineages != lineages)).cumsum(1).argmax(1).reshape(-1, 1),
                                         1).reshape(-1)
       
        # fill nan values with the root taxonomy value
        if root_on_empty is True:
            tax_vector[tax_vector != tax_vector] = cls.get_root_id()

        # return the taxa id's in the appropriate format
        if tax_vector.size == 1:
            return tax_vector[0]
        elif convert_list == True:
            return tax_vector.tolist()
        else:
            return tax_vector

    lineage_t = TypeVar('lineage_t',
                        List[str | float],
                        List[List[str | float]],
                        np.ndarray,
                        pd.DataFrame,
                        pd.Series)

    @staticmethod
    def fill_lineage_gaps(lineages: lineage_t) -> lineage_t:
        """Assign unannotated ranks in lineages with a higher rank annotation if present.
        This fills the gaps of lineages for which there is a high rank (e.g. genus) annotation,
        but no known annotation for some parent ranks.
        
        Note: This method is only present for compatibility between other taxonomy database
        objects. GTDB lineages have no gaps by definition.

        Args:
            lineages (List[str | float] | List[List[str | float]] | np.ndarray | pd.DataFrame | pd.Series):
                A list or matrix of lineages

        Returns:
            List[str | float] | List[List[str | float]] | np.ndarray | pd.DataFrame | pd.Series:
                The lineage list/matrix with filled gaps.
        """
        return lineages
    
    
    @lru_cache(maxsize=None)
    def id_to_standard_lineage(self, tax_id: str) -> Tuple[str | float, ...]:
        """Return the lineage from a given taxonomy id.

        Args:
            tax_id (str): Taxonomy id.

        Returns:
            Tuple[str | float]: Standard lineage array.
        """
        if not self.id_in_dataset(tax_id):
            return (np.nan,) * 7
        
        # if genome given, convert to species id
        if self.genome_id_in_dataset(tax_id):
            tax_id = self.genome_id_to_species_id(tax_id)
        
        # get lineage up to (and including) tax id
        lineage = self.lineage_dict[tax_id] + [tax_id]
        
        # if lineage does not go to species, add none values
        trailing_none = len(self.RANK_LIST) - len(lineage)
        return tuple(lineage + [np.nan] * trailing_none)
        
        
    def id_in_dataset(self, tax_id: str) -> bool:
        """Check if a given taxonomy id is present in the taxonomy database.

        Args:
            tax_id (str): Taxonomy id.

        Returns:
            bool: Test that id is in dataset.
        """
        # First check if genome given
        if self.genome_id_in_dataset(tax_id):
            return True
        else:
            return tax_id in self.lineage_dict.keys()
        
        
    def id_to_rank(self, tax_id: str | float) -> str | float:
        """Get the rank of a specified taxonomy id.

        Args:
            tax_id (str | float): Taxonomy id.

        Returns:
            str | float: Taxonomy rank.
        """
        # if genome given, convert to species id
        if self.genome_id_in_dataset(tax_id):
            tax_id = self.genome_id_to_species_id(tax_id)
        
        if not isinstance(tax_id, str):
            return np.nan
        
        # id is represented as x__abc, where x is the rank character ID
        rank_id = tax_id[0]
        
        return self.RANK_ID_TO_NAME[rank_id]
    
    
    def id_to_parent(self,
                     tax_id: str | float,
                     parent_rank: str) -> str | float:
        """Return a parent taxonomy id at a specified rank for a
        given taxonomy id.

        Args:
            tax_id (str | float): Taxonomy id.
            rank (str): Rank to fetch parent id from.

        Returns:
            str | float: Taxonomy id of parent.
        """
        # if genome given, convert to species id
        if self.genome_id_in_dataset(tax_id):
            tax_id = self.genome_id_to_species_id(tax_id)        

        # check that taxonomy id is valid
        if not isinstance(tax_id, str):
            return np.nan
        elif not self.id_in_dataset(tax_id):
            print(f"""
                Tax ID {tax_id} not recognized in gtdb taxonomy database.
                """)
            return np.nan
        
        # check that given taxonomy rank is valid
        if parent_rank not in self.RANK_LIST:
            raise ValueError(f"Invalid rank name supplied: '{parent_rank}'")
        
        # return the taxonomy name at the specified parent rank, if no annotation, return nan
        parent_rank_index = self.RANK_LIST.index(parent_rank)
        lineage = self.lineage_dict[tax_id]
        if len(lineage) < parent_rank_index:
            return np.nan
        elif len(lineage) == parent_rank_index:
            return tax_id
        else:
            return lineage[parent_rank_index]
    

    def id_to_parents(self, tax_id: str | float) -> np.ndarray:
        """Return array of parent taxa for a given taxonomy id.

        Args:
            tax_id (str | float): Taxonomy id.

        Returns:
            np.ndarray: Array of parent taxa.
        """
        # if genome given, convert to species id
        if self.genome_id_in_dataset(tax_id):
            tax_id = self.genome_id_to_species_id(tax_id)
        
        if not isinstance(tax_id, str) or not self.id_in_dataset(tax_id):
            return np.array([])
        
        return np.array(self.lineage_dict[tax_id] + [tax_id])
    
    
    def id_to_children(self, tax_id: str | float) -> np.ndarray:
        """Return array of all children that belong to a supplied taxonomy id.

        Args:
            tax_id (str | float): Taxonomy id.

        Returns:
            np.ndarray: Array of child taxonomy id's under a given tax_id.
        """
        # if genome given, convert to species id
        if self.genome_id_in_dataset(tax_id):
            tax_id = self.genome_id_to_species_id(tax_id)
            
        if not isinstance(tax_id, str) or not self.id_in_dataset(tax_id):
            return np.array([])
        
        return np.array(self.child_dict[tax_id])
        
        
    def taxa_to_lca(self,
                    tax_ids: List[str | float] | pd.Series | np.ndarray,
                    unknown_taxa: str = 'ignore') -> str | float:
        """Return the last common ancestor for a list of taxonomy id's.
        If invalid id's are passed, they can be ignored, an error may be
        raised or the root taxa can be returned. If no valid taxonomy id
        is present, an empty value is returned.

        Args:
            tax_ids (List[str | float] | np.ndarray | pd.Series): List,
                numpy array or series of tax ids.
            unknown_taxa (str, optional): Specify behavior if tax id is
                encountered that is absent from the database. Options:
                {'error', 'ignore', 'root', 'none'}. Defaults to 'ignore'.

        Raises:
            ValueError: Unknown taxa present in input
            ValueError: Invalid `unknown_taxa` value given.

        Returns:
            str | float: Last common ancestor for input taxa
        """
        # convert all array types into numpy array
        if isinstance(tax_ids, Sequence):
            tax_ids = np.array(tax_ids)
        elif isinstance(tax_ids, pd.Series):
            tax_ids = tax_ids.to_numpy()
            
        # if genomes given, convert to species id
        genome_to_id = lambda x: self.genome_id_to_species_id(x) if\
            self.genome_id_in_dataset(x) else x
        tax_ids = np.array([genome_to_id(i) for i in tax_ids])
            
        # remove nan values
        tax_ids = tax_ids[~np.array([i != i for i in tax_ids])]
        
        # check for unknown taxa. If present, either throw error, ignore or set lca to root
        known_taxa_array = np.array([self.id_in_dataset(i) for i in tax_ids])
        if False in known_taxa_array:
            if unknown_taxa.lower() == 'error':
                raise ValueError(f"Unknown taxa present in input")
            elif unknown_taxa.lower() == 'ignore':
                tax_ids = tax_ids[[i for i in known_taxa_array]]
            elif unknown_taxa.lower() == 'root':
                return self.ROOT_NAME
            elif unknown_taxa.lower() == 'none':
                return np.nan
            else:
                raise ValueError(f"Invalid handling of unknown taxa given: '{unknown_taxa}'")

        # if only one taxid valid, return this id. If no valid taxid remains, return nan
        if tax_ids.size == 0:
            return np.nan
        elif tax_ids.size == 1:
            return tax_ids[0]

        # retrieve lineages of all id's and store in list of lineage arrays
        lineages = [self.lineage_dict[tax_id] for tax_id in tax_ids]
        lineages = [lin + [tax] for lin, tax in zip(lineages, tax_ids)]

        # fetch highest rank id for which there is consensus among lineages
        return self.lineages_to_lca(lineages)
    
    
    def id_to_name(self, tax_id: str | float) -> str | float:
        """Convert taxonomy id to taxa name.

        Args:
            tax_id (str | float): Taxonomy id

        Returns:
            str | float: Taxonomy name
        """
        # if genome given, convert to species id
        if self.genome_id_in_dataset(tax_id):
            tax_id = self.genome_id_to_species_id(tax_id)
        
        if not isinstance(tax_id, str) or not self.id_in_dataset(tax_id):
            return np.nan
        else:
            return tax_id[3:]
        
        
    def lineage_id_to_name(self, id_lineage: List[str | float] | Tuple[str | float, ...]) -> Tuple[str | float, ...]:
        """Convert lineage of taxonomy id's to a lineage array of taxa names.

        Args:
            id_lineage (List[str | float] | Tuple[str | float, ...]): Lineage array of taxonomy id's.

        Returns:
            Tuple[str | float]: Lineage array of taxa names
        """
        return tuple([self.id_to_name(i) for i in id_lineage])


    @overload
    def name_to_id(self,
                   tax_name: str | float,
                   on_duplicates: Literal["all"],
                   print_fails: bool) -> str | float | List[str]:
        ...

    @overload
    def name_to_id(self,
                   tax_name: str | float,
                   on_duplicates: Literal["nan"],
                   print_fails: bool) -> str | float:
        ...
    
    def name_to_id(self,
                   tax_name: str | float,
                   on_duplicates: str="nan",
                   print_fails: bool=False) -> str | List[str] | float:
        """Convert organism name to taxonomy id. For the GTDB dataset, the 
        organism name is the same as the organism id, but without rank prefix.
        

        Args:
            tax_name (str | float): Taxonomy name.
            on_duplicates (str, optional): Specify function behavior if multiple
                id's encountered with same name. Options: {"nan", "all"}.
                'nan' will return nan when multiple id's encountered,
                'all' will return list of taxonomy id's.
            print_fails (bool, optional): Print cases where tax id retrieval
                fails to stdout. Defaults to False.

        Raises:
            ValueError: Invalid `on_duplicates` value given.

        Returns:
            str | List[str] | float: Taxonomy id('s) coupled to name
        """
        if print_fails is True:
            print("Taxonomy name to id conversion\n\nFailed conversions:")    
        
        if isinstance(tax_name, float):
            if print_fails is True:
                print("Empty input encountered...")
            return np.nan
        
        # if taxonomy name is NCBI genome id, return its own value
        if self.genome_id_in_dataset(tax_name):
            return self.genome_id_to_species_id(tax_name)
        
        potential_ids = np.array([i + f"__{tax_name}" for i in self.RANK_ID_TO_NAME.keys()])
        id_in_dataset = [self.id_in_dataset(i) for i in potential_ids]
        valid_ids = potential_ids[id_in_dataset]
        
        # return id's or nan based on selected parameter options
        if valid_ids.size == 0:
            if print_fails is True:
                print(f"No valid id's for {tax_name}")
            return np.nan
        elif valid_ids.size == [1]:
            return valid_ids[0]
        elif on_duplicates.lower() == "nan":
            if print_fails is True:
                print(f"Multiple id's for {tax_name}: {valid_ids}")
            return np.nan
        elif on_duplicates.lower() == "all":
            return valid_ids.tolist()
        else:
            raise ValueError(f"Invalid `on_duplicates` value given '{on_duplicates}'. Choose from {'nan', 'all'}.")

    
    @classmethod
    def lineages_to_lca(cls, lineages: Sequence[Sequence[str]]) -> str | float:
        """Retrieve the last common ancestor tax id from a list of taxonomy
        lineages. This method parses an array of lineage strings retrieved from
        the lineage dataset. It will return the highest ranked taxonomy for which
        there is consensus.

        Args:
            lineages (List[List[str | float] | Tuple[str | float]]): List of lineages.

        Returns:
            str | float: Last common ancestor of lineages.
        """
        # return nan if 0 lineages are given as input
        if len(lineages) == 0:
            return np.nan
        
        lca = cls.ROOT_NAME             # start with root
        min_len = min(len(x) for x in lineages)
        
        # parse lineages and update lca as long as there is full consensus
        for i in range(min_len):
            if any([x[i] != x[i] for x in lineages]):
                break
            if len({x[i] for x in lineages}) == 1:
                lca = lineages[0][i]
        
        return lca
    
    # NOTE: Only GtdbTaxonomy methods, not part of abstract base class
    def genome_id_in_dataset(self, genome_id: str | float) -> bool:
        """Check that genome id is present within the gtdb taxonomy dataset.

        Args:
            genome_id (str): Genome id name.

        Returns:
            bool: True if present in the dataset
        """
        if isinstance(genome_id, float):
            return False

        # remove genbank/refseq prefixes to be consistent with taxonomy db format
        genome_id = genome_id.removeprefix("GB_").removeprefix("RS_")
        
        return genome_id in self.genome_species_map.keys()
    
    def genome_id_to_species_id(self,
                                genome_id: str | float) -> str | float:
        """Return species id of genome id.

        Args:
            genome_id (str): Genome id name.

        Returns:
            str | None: Species id, if present in dataset
        """
        # remove genbank/refseq prefixes to be consistent with taxonomy db format
        genome_id = genome_id.removeprefix("GB_").removeprefix("RS_")
        
        return self.genome_species_map.get(genome_id, np.nan)
        
        
        
    ################################################################################ 
    # Constructor methods
    ################################################################################ 
        
    
    @classmethod
    def from_tsv_folder(cls, tsv_dir: str | Path) -> "GtdbTaxonomy":
        """Import gtdb taxonomy datasets from folder.
        imports "bac120_taxonomy.tsv", "ar53_taxonomy.tsv" files.

        Args:
            tsv_dir (str | Path): Location folder.

        Raises:
            FileNotFoundError: Database files not present in location.

        Returns:
            GtdbTaxonomy: GtdbTaxonomy object.
        """
        # convert to Path if str supplied
        tsv_dir = Path(tsv_dir)
        
        # check for presence of all files
        bact_file, arch_file = [Path(tsv_dir, x) for x in cls.FILE_NAMES]
        if all([x.exists() for x in [bact_file, arch_file]]):
            # call file constructor to import files
            return cls.from_tsv_files(bact_file, arch_file)
        else:
            raise FileNotFoundError(f"one of {cls.FILE_NAMES} not present in supplied directory")

        
    @classmethod
    def from_tsv_files(cls,
                       bacteria_file: str | Path,
                       archaea_file: str | Path) -> "GtdbTaxonomy":
        bacteria_file = Path(bacteria_file)
        archaea_file = Path(archaea_file)
        
        # import bacteria and archaea into dataframe and concatenate them
        bac_df = pd.read_csv(bacteria_file,
                             sep="\t",
                             engine="python", 
                             names=["genome", "lineage"])
        arch_df = pd.read_csv(archaea_file,
                              sep="\t", 
                              engine="python",
                              names=["genome", "lineage"])
        
        taxonomy_df = pd.concat([bac_df, arch_df], axis=0).reset_index(drop=True)
        lineage_df = cls.__import_lineage_data(taxonomy_df)
        
        # drop genbank/refseq prefixes
        taxonomy_df.loc[:, 'genome'] = taxonomy_df['genome']\
            .str.removeprefix("GB_")\
            .str.removeprefix("RS_")
                
        
        # construct lineage dictionary: Dict[tax_id: Lineage]
        lineage_dict = dict()
        child_dict = defaultdict(lambda: [])
        genome_dict = dict(
            zip(taxonomy_df['genome'],
                taxonomy_df['lineage'].apply(lambda x: x.split(';')[-1])
                )
        )
        
        # Parse dataframe row-by-row, species-to-superkingdom
        row_size = lineage_df.shape[1] - 1
        for row_index, row in lineage_df.iterrows():
            for elem_index, element in enumerate(row.iloc[::-1].values):
                row_index = row_size - elem_index
                
                # add offspring taxonomies to child_dict
                if row_index != row_size and \
                row.iat[row_index+1] not in child_dict[element]:
                    child_dict[element].append(row.iat[row_index+1])
                
                # if taxonomy id is present in the dict, the id's after are also present
                if element in lineage_dict.keys():
                    break
                
                # add lineage of selected taxonomy to dict
                lineage_dict[element] = row.iloc[:row_index].to_list()

        return cls(lineage_dict, dict(child_dict), genome_dict)        


    @classmethod
    def __import_lineage_data(cls,
                              dataset: pd.DataFrame) -> pd.DataFrame:
        """Process gtdb taxonomy file and return dataset of unique lineages.

        Args:
            dataset (pd.DataFrame): Dataframe of gtdb taxonomy tsv

        Returns:
            pd.DataFrame: Lineage dataset.
        """
        # split lineage string into rank cells, ignore 'domain' rank
        dataset = pd.DataFrame.from_records(dataset["lineage"].str.split(";"),
                                            columns=cls.RANK_LIST)
        return dataset.drop_duplicates()

    