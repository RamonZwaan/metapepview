"""This module describes the base class for the definition of classes
that manage different types of taxonomy databases (like ncbi or gtdb).

The base class will ensure that all derived classes share the same
interface so that new database types will integrate into existing
functions without code rewrites.
"""

from __future__ import annotations
from abc import ABCMeta, abstractmethod
from functools import lru_cache
from pathlib import Path
from typing import Tuple, List, Any, overload, TypeVar
from warnings import warn

import numpy as np
import pandas as pd



class TaxonomyDatabase(metaclass=ABCMeta):
    """Abstract base class definition that defines a common interface for different
    taxonomy databases to adhere to. This ensures that different subclasses can be
    used interchangeably in the functions that call this base class type.

    Args:
        metaclass (_type_, optional): _description_. Defaults to ABCMeta.
    """ 
    
    @classmethod
    @abstractmethod
    def get_root_id(cls) -> Any:
        ...

    
    @classmethod
    @abstractmethod
    def lineage_to_id(cls,
                      lineages: List[Any] | List[List[Any]] | np.ndarray,
                      root_on_empty: bool = False) -> List[Any] | np.ndarray | Any:
        """Convert lineage vector or matrix of lineages to the highest valid taxonomy id for each lineage

        Args:
            lineages (List[Any] | List[List[Any]] | np.ndarray): Array or matrix of lineages.

        Returns:
            List[Any] | np.ndarray | Any: Taxonomy id for each lineage
        """
        ...

    lineage_t = TypeVar('lineage_t',
                        List[str | float],
                        List[List[str | float]],
                        np.ndarray,
                        pd.DataFrame,
                        pd.Series)

    @staticmethod
    @abstractmethod
    def fill_lineage_gaps(lineages: lineage_t) -> lineage_t:
        """Assign unannotated ranks in lineages with a higher rank annotation if present.
        This fills the gaps of lineages for which there is a high rank (e.g. genus) annotation,
        but no known annotation for some parent ranks.

        Args:
            lineages (List[Any] | List[List[Any]] | np.ndarray): A list or matrix of lineages

        Returns:
            List[Any] | List[List[Any]] | np.ndarray: The lineage list/matrix with filled gaps.
        """
        ...
    
    @abstractmethod
    def id_to_standard_lineage(self, tax_id: Any) -> Tuple[Any]:
        """Return the lineage at the standard rank levels for a given taxonomy id.
        Any taxa inbetween valid ranks (sub-, clade, etc.) are filtered from the
        lineage.

        Args:
            tax_id (Any): Taxonomy id.

        Returns:
            Tuple[Any]: Standard lineage array.
        """
        ...
        
    @abstractmethod
    def id_in_dataset(self, tax_id: Any) -> bool:
        """Check if a given taxonomy id is present in the taxonomy database.

        Args:
            tax_id (Any): Taxonomy id.

        Returns:
            bool: Test that id is in dataset.
        """
        ...
        
    @abstractmethod
    def id_to_rank(self, tax_id: Any) -> str:
        """Get the rank of a specified taxonomy id.

        Args:
            tax_id (Any): Taxonomy id.

        Returns:
            str: Taxonomy rank.
        """
        ...
    
    @abstractmethod
    def id_to_parent(self, tax_id: Any, rank: str) -> Any:
        """Return a parent taxonomy id at a specified rank for a
        given taxonomy id.

        Args:
            tax_id (Any): Taxonomy id.
            rank (str): Rank to fetch parent id from.

        Returns:
            Any: Taxonomy id of parent.
        """
        ...
    
    @abstractmethod
    def id_to_parents(self, tax_id: Any) -> np.ndarray:
        """Return array of parent taxa for a given taxonomy id.

        Args:
            tax_id (Any): Taxonomy id.

        Returns:
            np.ndarray: Array of parent taxa.
        """
        ...
    
    @abstractmethod
    def id_to_children(self, tax_id: Any) -> np.ndarray:
        """Return array of all children that belong to a supplied taxonomy id.

        Args:
            tax_id (Any): Taxonomy id.

        Returns:
            np.ndarray: Array of child taxonomy id's under a given tax_id.
        """
        ...
        
    @abstractmethod
    def taxa_to_lca(self, tax_ids: List[Any] | pd.Series | np.ndarray) -> Any:
        """Return the last common ancestor from an array of taxonomy id's.

        Args:
            tax_ids (List[Any] | pd.Series | np.ndarray): Array of taxonomy id's

        Returns:
            Any: Last common ancestor.
        """
        ...
    
    @abstractmethod
    def id_to_name(self, tax_id: Any) -> str:
        """Convert taxonomy id to taxa name.

        Args:
            tax_id (Any): Taxonomy id

        Returns:
            str:  Taxa name
        """
        ...
        
    @abstractmethod
    def lineage_id_to_name(self, id_lineage: List[Any] | Tuple[Any, ...]) -> List[str | None] | Tuple[str | None]:
        """Convert lineage of taxonomy id's to a lineage array of taxa names.

        Args:
            id_lineage (List[Any] | Tuple[Any, ...]): Lineage array of taxonomy id's.

        Returns:
            Tuple[str | None]: Lineage array of taxa names
        """
        ...
    
    @abstractmethod
    def name_to_id(self, tax_name: str) -> Any:
        """Convert organism name to taxonomy id

        Args:
            tax_name (str): Organim or taxa name

        Returns:
            Any: Taxonomy id
        """
        ...
    
    @staticmethod
    @abstractmethod
    def lineages_to_lca(lin_series: List[Any] | Tuple[Any, ...] | pd.Series) -> Any:
        """Retrieve the last common ancestor tax id from a list of taxonomy
        lineages. This function does not require any external taxonomy
        datasets.

        Args:
            lin_series (Union[List[str], Tuple[str, ...]]): List of lineages.
            root_taxa (Any): root taxonomy.

        Returns:
            int: Last common ancestor of lineages.
        """
        ...




# class NcbiTaxonomy(TaxonomyDatabase):
    
#     FILE_NAMES = ('nodes.dmp', 'names.dmp', 'taxidlineage.dmp')
    
#     SERVER_URL = "ftp.ncbi.nih.gov"
#     SERVER_FS = "/pub/taxonomy/new_taxdump"
    
#     NODES_COLUMNS = ("taxonomy_id", "parent_id", "rank", "code", "div_id",
#                      "div_flag", "code_id", "gc_flag", "mit_code_id", "mgc_flag",
#                      "hid_flag", "root_flag", "comments", "plast_code_id",
#                      "pgc_flag", "spec_sp", "hydr_code_id", "hgc_flag")
#     NAMES_COLUMNS = ('taxonomy_id', 'name_txt', 'unique_name', 'name_class')
#     LINEAGE_COLUMNS = ('taxonomy_id', 'lineage')
    
#     STANDARD_RANK_LIST = ["superkingdom",
#                           "phylum",
#                           "class",
#                           "order",
#                           "family",
#                           "genus",
#                           "species"]
    
    
#     ROOT_NAME = 1.0
    
#     def __init__(self,
#                  node_df: pd.DataFrame,
#                  name_df: pd.DataFrame,
#                  lineage_df: pd.DataFrame):
#         self.node_df = node_df
#         self.name_df = name_df
#         self.lineage_df = lineage_df
        
#         # cache variables
#         self.valid_taxa = None      # set of tax id's to evaluate presence in O(1)
#         self.taxa_rank_dict = None  # dict that couples id to rank


#     @classmethod
#     def get_root_id(cls) -> float:
#         """Return the global root id for the ncbi taxonomy database.

#         Returns:
#             float: The global taxonomy root.
#         """
#         return cls.ROOT_NAME
    

#     @lru_cache(maxsize=None)
#     def id_to_standard_lineage(self,
#                                tax_id: int | float) -> Tuple[int | float, ...]:
#         """Return the lineage for a given taxonomy id. The lineage will only
#         contain taxa from the default ranks.

#         Args:
#             tax_id (int | None): Taxonomy id.

#         Returns:
#             Tuple[int | None, ...]: Lineage array
#         """
#         # return empty lineage if a missing value or the root taxonomy is given
#         if np.isnan(tax_id):
#             return (np.nan,)*7
#         if tax_id == 1 or not self.id_in_dataset(tax_id):
#             return (np.nan,)*7
        
#         # retrieve complete tax id lineage
#         complete_lin = self.__id_to_lineage(tax_id)
        
#         # resolve rank names for all lineage id's
#         rank_to_taxid = dict()
#         for lin_tax in complete_lin:
#             rank = self.id_to_rank(lin_tax)
#             if rank is not None:
#                 rank_to_taxid[rank] = lin_tax
        
#         # fetch all standard rank id's, for ranks without value, return nan
#         standard_lin = tuple(rank_to_taxid.get(i, np.nan)
#                         for i in self.STANDARD_RANK_LIST)
        
#         return standard_lin
    
#     @overload
#     @classmethod
#     def lineage_to_id(cls,
#                       lineages: List[int | float],
#                       root_on_empty: bool) -> int | float:
#         ...

#     @overload
#     @classmethod
#     def lineage_to_id(cls,
#                       lineages: List[List[int | float]],
#                       root_on_empty: bool) -> List[int | float]:
#         ...
        
#     @overload
#     @classmethod
#     def lineage_to_id(cls,
#                       lineages: np.ndarray,
#                       root_on_empty: bool) -> np.ndarray | int | float:
#         ...

#     @classmethod
#     def lineage_to_id(cls,
#                       lineages: List[int | float] | List[List[int | float]] | np.ndarray,
#                       root_on_empty: bool = True) -> List[int | float] | np.ndarray | int | float:
#         """Convert lineage array or matrix of lineages to the highest valid taxonomy id for each lineage.

#         Args:
#             lineages (List[int | float] | List[List[int | float]] | np.ndarray): Array or matrix of lineages.
#             root_on_empty (bool, Optional): Specify if root taxonomy should be returned if an empty lineage
#                 is supplied. Else, a NaN value will be returned. Defaults to True.
#             root_value (Any, Optional): Root taxonomy id. Will be returned when empty lineage is supplied if
#                 specified from `root_on_empty`. Defaults to 1.

#         Returns:
#             List[int, float] | np.ndarray | int | float: Taxonomy id for each lineage.
#         """
#         # manage different formats of input data
#         convert_list = False
#         if isinstance(lineages, list):
#             convert_list = True
#             lineages = np.array(lineages)
        
#         # specify output format (array if input is 2d, scalar if input is 1d)
#         return_array = False
#         if lineages.ndim == 2:
#             return_array = True
#         elif lineages.ndim == 1:
#             lineages.reshape(1, -1)
#         else:
#             raise ValueError(f"Invalid lineage data, dimensionality of '{lineages.ndim}' not supported.")
        
#         # for any lineage, take the last valid taxonomy value        
#         tax_vector =  np.take_along_axis(lineages,
#                                          (~np.isnan(lineages)).cumsum(1).argmax(1).reshape(-1, 1),
#                                          1).reshape(-1)
       
#         # fill nan values with the root taxonomy value
#         if root_on_empty is True:
#             tax_vector = np.nan_to_num(tax_vector, nan=cls.ROOT_NAME)

#         # return the taxa id's in the appropriate format
#         if tax_vector.size == 1 and return_array is False:
#             return tax_vector[0]
#         elif convert_list == True:
#             return tax_vector.tolist()
#         else:
#             return tax_vector
    
#     lineage_t = TypeVar('lineage_t',
#                         List[str | float],
#                         List[List[str | float]],
#                         np.ndarray,
#                         pd.DataFrame,
#                         pd.Series)

#     @staticmethod
#     def fill_lineage_gaps(lineages: lineage_t) -> lineage_t:
#         """Assign unannotated ranks in lineages with a higher rank annotation if present.
#         This fills the gaps of lineages for which there is a high rank (e.g. genus, species)
#         annotation, but no known annotation for some parent ranks.

#         Args:
#             lineages (List[int | float] | List[List[int | float]] | np.ndarray | pd.DataFrame | pd.Series):
#                 A list or matrix of lineages

#         Returns:
#             List[int | float] | List[List[int | float]] | np.ndarray | pd.DataFrame | pd.Series:
#                 The lineage list/matrix with filled gaps.
#         """
#         # variables to specify format of original parameter, used to convert back to expected value
#         convert_list = False
#         convert_dim = False
#         convert_df = False
#         convert_series = False
#         df_cols, df_index = None, None
#         series_name, series_index = None, None

#         if isinstance(lineages, pd.DataFrame):
#             convert_df = True
#             df_cols = lineages.columns
#             df_index = lineages.index
#             lineages_arr = lineages.to_numpy()    
#         elif isinstance(lineages, pd.Series):
#             convert_series = True
#             series_name = lineages.name
#             series_index = lineages.index
#             lineages_arr = np.array(lineages.to_list())
#         # if list supplied convert to array, remember to convert back to expected type
#         elif isinstance(lineages, list):
#             convert_list = True
#             lineages_arr = np.array(lineages)
#         elif isinstance(lineages, np.ndarray):
#             lineages_arr = lineages
        
#         # vector array should be converted to 2-d array, also here remember to convert back
#         if lineages_arr.ndim == 1:
#             convert_dim = True
#             lineages_arr = lineages_arr.reshape(1, -1)
        
#         # Create array that stores last valid taxa for each lineage
#         current_tax_id = np.empty(lineages_arr.shape[0])
#         current_tax_id[:] = np.nan
        
#         # iterate from species to superkingdom, replace nan with valid taxa if present
#         for rank_index in range(lineages_arr.shape[1] - 1, -1, -1):
#             column = lineages_arr[:, rank_index]  

#             notnan = ~np.isnan(column)
#             current_tax_id[notnan] = column[notnan]
            
#             # fill nan with higher ranked taxa       
#             lineages_arr[:, rank_index] = current_tax_id
        
#         # convert to same format as input parameter
#         if convert_dim == True:
#             lineages_arr = lineages_arr.reshape(-1)
#         if convert_list == True:
#             lineages_arr = lineages_arr.tolist()
#         if convert_df == True:
#             lineages_arr = pd.DataFrame(lineages_arr, index=df_index, columns=df_cols)
#         if convert_series == True:
#             lineages_arr = pd.Series(list(lineages_arr), index=series_index, name=series_name)
         
#         return lineages_arr # type:ignore (output type does not match TypeVar exactly)
        
        
#     def id_in_dataset(self, tax_id: int | float) -> bool:
#         """Assert that taxonomy id is present in the Ncbi Taxonomy database.

#         Args:
#             tax_id (int | float): Taxonomy id.

#         Returns:
#             bool: True if id is present in the database
#         """
#         # if an empty value is given (NaN) return False
#         if np.isnan(tax_id):
#             return False
        
#         # store tax id's in set for quick comparison
#         if self.valid_taxa is None:
#             self.valid_taxa = set(self.node_df.index.to_list())
        
#         return tax_id in self.valid_taxa
        
        
#     def id_to_rank(self, tax_id: int | float) -> str | None:
#         """Return the taxonomy rank to which a taxonomy id belongs.
#         If taxonomy id is not present in dataset, or if empty value
#         is given. Return None.

#         Args:
#             tax_id (int): Taxonomy id.

#         Returns:
#             str | None: Taxonomy rank if present.
#         """
#         if np.isnan(tax_id):
#             return None
        
#         try:
#             return self.node_df.at[tax_id, "rank"]
#         except KeyError:
#             return None

        
#     def id_to_parent(self,
#                      tax_id: int | float,
#                      parent_rank: str,
#                      absent_rank_behavior: str = "nan") -> int | float:
#         """Retrieve taxonomy id at specified parent rank for a given tax id.
#         If no id is present at the desired rank for a given taxonomy, the
#         method will do the following, based on user settings specified in
#         `absent_rank_behavior`:
        
#             - take the first lower rank that has a known annotation ("lower").
#             - the first higher rank that has a known annotation ("upper").
#             - return nan ("nan").

#         Args:
#             tax_id (int): taxonomy id
#             parent_rank (str): Desired rank to retrieve taxonomy id
#             absent_rank_behavior (str, optional): Behavior setting if rank
#             is absent for lineage. Options: {"upper", "lower", "nan}. Defaults to "nan".

#         Raises:
#             ValueError: Invalid rank name supplied.

#         Returns:
#             int | float: Taxonomy id of parent.
#         """
#         # if empty value given, return nan
#         if np.isnan(tax_id):
#             return np.nan
        
#         # check rank valid
#         if parent_rank not in self.STANDARD_RANK_LIST:
#             raise ValueError(f"Invalid rank name supplied: '{parent_rank}'")
        
#         # list of ranks below the specified rank
#         lower_ranks = self.STANDARD_RANK_LIST[:self.STANDARD_RANK_LIST.index(parent_rank)]
        
#         # traverse lineage
#         last_standard_id = tax_id
#         current_id = tax_id
#         while current_id != 1:
#             # taxonomy id not in dataset
#             if self.id_in_dataset(tax_id) is False:
#                 print(f"""
#                     Tax ID {tax_id} not recognized in ncbi taxonomy database,
#                     input data might match to a different ncbi taxonomy version.
#                     """)
#                 return np.nan
            
#             row = self.node_df.loc[current_id]

#             # Correct taxonomy id encountered
#             if row["rank"] == parent_rank:
#                 return current_id
            
#             # behavior if lower rank encountered
#             if row["rank"] in lower_ranks:
#                 if absent_rank_behavior == "lower":
#                     return current_id
#                 elif absent_rank_behavior == "upper":
#                     return last_standard_id
#                 else:
#                     return np.nan
                    
#             # update standard id if standard rank encountered
#             if row["rank"] in self.STANDARD_RANK_LIST:
#                 last_standard_id = current_id
            
#             # move back in lineage
#             current_id = row["parent_id"]

#         # no tax id encountered through lineage
#         return np.nan
    
    
#     def id_to_parents(self, tax_id: int | float) -> np.ndarray:
#         """Retrieve array of parent taxa for a given tax id.

#         Args:
#             tax_id (int | float): Taxonomy id

#         Returns:
#             np.ndarray: Array of parent taxa
#         """
#         # for empty value, or unknown value, return empty array
#         if np.isnan(tax_id) or not self.id_in_dataset(tax_id):
#             return np.array([])
        
#         # get complete lineage and return up to the tax id
#         lineage = self.__id_to_lineage(tax_id)
#         index_current_tax = lineage.index(tax_id)
        
#         return np.array(lineage[:index_current_tax + 1])
    

#     def id_to_children(self, root_id: int | float | List[int | float]) -> np.ndarray:
#         """Retrieve all taxonomy children present in the database from
#         a list of taxonomy id's.

#         Args:
#             root_id (int | List[int]): List of taxonomy id's to retrieve
#                 child id's from

#         Returns:
#             np.ndarray: Array of child id's
#         """
#         # return empty array for empty value
#         if isinstance(root_id, int | float) and np.isnan(root_id):
#             return np.ndarray([])
        
#         if isinstance(root_id, int | float):
#             root_id = [root_id]
        
#         # give warning if invalid tax id's are encountered
#         valid_ids = [self.id_in_dataset(id) for id in root_id]
#         if not all(valid_ids):
#             print("Warning, invalid taxonomy id's supplied. Will be ignored")
        
#         local_nodes = np.array(root_id)[valid_ids]
#         local_taxa = local_nodes
#         # parse iteratively through child tax_id's
#         while True:
#             node_set = self.node_df[self.node_df["parent_id"].isin(local_taxa)]
            
#             # stop loop if there are no offspring taxonomy ids
#             if node_set.shape[0] == 0:
#                 break
            
#             # get child taxa from the set and add to nodes list
#             local_taxa = node_set.index.to_list()
#             local_nodes = np.append(local_nodes, np.array(local_taxa))
        
#         return np.unique(local_nodes)
    
       
#     @lru_cache(maxsize=None)    
#     def taxa_to_lca(self,
#                     tax_ids: List[int | float] | np.ndarray | pd.Series,
#                     unknown_taxa: str='ignore') -> int | float:
#         """Return the last common ancestor for a list of taxonomy id's.
#         If invalid id's are passed, they can be ignored, an error may be
#         raised or the root taxa can be returned. If no valid taxonomy id
#         is present, an empty value is returned.

#         Args:
#             tax_ids (List[int | float] | np.ndarray | pd.Series): List,
#                 numpy array or series of tax ids.
#             unknown_taxa (str, optional): Specify behavior if tax id is
#                 encountered that is absent from the database. Options:
#                 {'error', 'ignore', 'root', 'none'}. Defaults to 'ignore'.

#         Raises:
#             ValueError: Unknown taxa present in input
#             ValueError: Invalid `unknown_taxa` value given.

#         Returns:
#             int: Last common ancestor for input taxa
#         """
#         # convert all array types into pandas series
#         if not isinstance(tax_ids, pd.Series):
#             tax_ids = pd.Series(tax_ids)
            
#         # remove nan and convert to int
#         tax_ids = tax_ids.dropna().astype('Int64')
        
#         # check for unknown taxa. If present, either throw error, ignore or set lca to root
#         unknown_taxa_array = (~tax_ids.transform(self.id_in_dataset)).to_list()
#         if any(unknown_taxa_array):
#             if unknown_taxa.lower() == 'error':
#                 raise ValueError(f"Unknown taxa present in input")
#             elif unknown_taxa.lower() == 'ignore':
#                 tax_ids = tax_ids[[not i for i in unknown_taxa_array]]
#             elif unknown_taxa.lower() == 'root':
#                 return 1
#             elif unknown_taxa.lower() == 'none':
#                 return np.nan
#             else:
#                 raise ValueError(f"Invalid handling of unknown taxa given: '{unknown_taxa}'")

#         # if only one taxid valid, return this id. If no valid taxid remains, return nan
#         if tax_ids.size == 0:
#             return np.nan
#         elif tax_ids.size == 1:
#             return tax_ids.iloc[0]

#         # retrieve lineages of all id's and store in list of lineage arrays
#         lineages = self.lineage_df.loc[tax_ids, 'lineage'].to_list()
#         lineages = [lin + f"{tax}"
#                     if len(lin) > 0
#                     else str(tax)
#                     for lin, tax in zip(lineages, tax_ids)]

#         # fetch highest rank id for which there is consensus among lineages
#         return self.lineages_to_lca(lineages)
    
    
#     @staticmethod
#     def lineages_to_lca(lin_series: List[str] | Tuple[str, ...] | pd.Series) -> int | float:
#         """Retrieve the last common ancestor tax id from a list of taxonomy
#         lineages. This method parses an array of lineage strings retrieved from
#         the lineage dataset. It will return the highest ranked taxonomy for which
#         there is consensus.

#         Args:
#             lin_series (Union[List[Any], Tuple[Any, ...]]): List of lineages.

#         Returns:
#             Any: Last common ancestor of lineages.
#         """
#         if not isinstance(lin_series, pd.Series):
#             lin_series = pd.Series(lin_series)
        
#         # If no lineage given at all, return nan
#         if lin_series.size == 0:
#             return np.nan
        
#         # split string with spaced taxa into list
#         lineages = [lin.split(" ") for lin in lin_series]

#         lca = 1             # start with root
#         min_len = min(len(x) for x in lineages)
        
#         # parse lineages and update lca as long as there is full consensus
#         for i in range(min_len):
#             if len({x[i] for x in lineages}) == 1:
#                 lca = lineages[0][i]
        
#         return int(lca)
    
#     @lru_cache(maxsize=10000)
#     def id_to_name(self, tax_id: int | float) -> str:
#         """Return the scientific taxonomy name for a given taxonomy id.

#         Args:
#             tax_id (int | float): taxonomy id

#         Returns:
#             str | None: Scientific name for taxonomy id.
#         """
#         # may receive np.nan values, directly return NoneType
#         if np.isnan(tax_id) == True:
#             return "undefined"
        
#         # return undefined if id not in dataset
#         if self.id_in_dataset(int(tax_id)) is False:
#             return 'undefined'
        
#         organism = self.name_df[self.name_df.index == int(tax_id)]
        
#         # when multiple rows present, get scientific name
#         if organism.shape[0] >= 1:
#             organism = organism[organism['name_class'] == "scientific name"]
        
#         if organism.shape[0] == 0:
#             return "undefined"
#         else:
#             return organism["name_txt"].item()
            
            
#     def lineage_id_to_name(self, id_lineage: Tuple[Any, ...]) -> Tuple[str]:
#         """Convert array of taxonomy id's to array of taxonomy names.
#         This method is used to convert complete lineages.

#         Args:
#             id_lineage (Tuple[Any, ...]): Lineage array of taxonomy id's.

#         Returns:
#             Tuple[str]: Lineage array of taxonomy names.
#         """
#         return tuple(self.id_to_name(i) for i in id_lineage)
    
    
#     @lru_cache(maxsize=None)
#     def name_to_id(self, tax_name: str | None) -> int | float:
#         """Convert taxonomy name to taxonomy id.

#         Args:
#             tax_name (str): Taxonomy name.

#         Returns:
#             int | float: Taxonomy id.
#         """
#         if tax_name is None:
#             return np.nan
        
#         # get tax_id of given name
#         tax_id = self.name_df[self.name_df['name_txt'] == tax_name]
        
#         if tax_id.empty is True:
#             return np.nan
#         elif tax_id.shape[0] > 1:
#             warn(f"Multiple tax id's observed for taxonomy name '{tax_name}'.\nLCA taken...",
#                  UserWarning)
#             return self.taxa_to_lca(tax_id.index.to_series())
#         else:
#             return tax_id.index.item()
    
    
#     # constructor methods
    
#     @classmethod
#     def from_dmp_folder(cls, dmp_dir: str | Path) -> "NcbiTaxonomy":
#         """Import ncbi taxonomy datasets from dmp folder. NcbiTaxonomy
#         imports 'nodes.dmp', 'names.dmp', and 'taxidlineage.dmp' files.

#         Args:
#             dmp_dir (str | Path): Location dmp folder.

#         Raises:
#             FileNotFoundError: Database files not present in location.

#         Returns:
#             NcbiTaxonomy: NcbiTaxonomy object.
#         """
#         # convert to Path if str supplied
#         dmp_dir = Path(dmp_dir)
        
#         # check for presence of all files
#         nodes, names, lineages = [Path(dmp_dir, x) for x in cls.FILE_NAMES]
#         if all([x.exists() for x in [nodes, names, lineages]]):
#             # call file constructor to import files
#             return cls.from_dmp_files(nodes, names, lineages)
#         else:
#             raise FileNotFoundError(f"one of {cls.FILE_NAMES} not present in supplied directory")
            
#     @classmethod
#     def from_dmp_files(cls,
#                        nodes_file: str | Path,
#                        names_file: str | Path,
#                        lineage_file: str | Path) -> "NcbiTaxonomy":
#         """Import ncbi taxonomy database from dataset files. NcbiTaxonomy
#         imports 'nodes.dmp', 'names.dmp', and 'taxidlineage.dmp' files.

#         Args:
#             nodes_file (str | Path): Location 'nodes.dmp' file.
#             names_file (str | Path): Location 'names.dmp' file.
#             lineage_file (str | Path): Location 'taxidlineage.dmp' file.

#         Returns:
#             NcbiTaxonomy: NcbiTaxonomy object.
#         """
        
#         # import and process nodes dataset
#         nodes_df = cls.__import_dmp(Path(nodes_file), cls.NODES_COLUMNS, encoder="utf8")
#         nodes_df = nodes_df[["taxonomy_id", "parent_id", "rank"]]       # omit unneeded columns
#         nodes_df = nodes_df.apply(pd.to_numeric, errors='ignore')
#         nodes_df["rank"] = nodes_df["rank"].astype("category")
#         nodes_df.set_index(["taxonomy_id"], inplace=True)
        
#         # import names dataset
#         names_df = cls.__import_dmp(Path(names_file), cls.NAMES_COLUMNS, encoder='latin-1') 
#         names_df = names_df[names_df['name_class'] == 'scientific name']    # keep only scientific names
#         names_df["taxonomy_id"] = names_df["taxonomy_id"].astype(int)
#         names_df.set_index(["taxonomy_id"], inplace=True)
        
#         # import lineage dataset
#         lineage_df = cls.__import_dmp(Path(lineage_file), cls.LINEAGE_COLUMNS, encoder='latin-1')
#         lineage_df["taxonomy_id"] = pd.to_numeric(lineage_df["taxonomy_id"])
#         lineage_df.set_index(["taxonomy_id"], inplace=True)
        
#         return cls(nodes_df, names_df, lineage_df)
    
    
#     # private methods
#     @staticmethod
#     def __import_dmp(dmp_file: Path, columns: Tuple[str, ...], encoder: str="utf8") -> pd.DataFrame:
#         """Parse dmp file formats and store into dataframe.

#         Args:
#             dmp_file (Path): Location of dmp file
#             columns (Tuple[str, ...]): Columns of file contents.
#             encoder (str, optional): Encoder to read dmp file with. Defaults to "utf8".

#         Raises:
#             FileNotFoundError: File does not exist

#         Returns:
#             pd.DataFrame: Imported dataset.
#         """
#         # check if path location is valid
#         if not dmp_file.exists():
#             raise FileNotFoundError("Path to nodes file does not exist.")
            
#         data = list()
        
#         # parse dmp file
#         with open(dmp_file, 'r', encoding=encoder) as read_file:
#             for line in read_file.readlines():
#                 line = line.replace("\t|\n", "") # remove trailing |
#                 line = line.split("\t|\t")
#                 data.append(line)
        
#         return pd.DataFrame(data, columns=columns)

    
#     def __id_to_lineage(self, tax_id: int | float) -> Tuple[int, ...]:
#         """Return the full lineage for a given taxonomy id. The lineage
#         includes the given taxonomy id.

#         Args:
#             tax_id (int | float): Taxonomy id.

#         Returns:
#             Tuple[int, ...]: Taxonomy id lineage.
#         """
#         # tax_id's 1 (root) and 131567 (cellular organisms) do not have lineages in dataset
#         if tax_id == 1:
#             return (1,)
#         elif tax_id == 131567:
#             return (1, 131567)
#         elif not self.id_in_dataset(tax_id):
#             return tuple()
#         else:
#             lineages = self.lineage_df.loc[tax_id, ['lineage']].item()
#             parents = lineages.split(" ")[:-1]
#             lineage_list = [1] + parents + [tax_id]
            
#         return tuple(int(_) for _ in lineage_list)
