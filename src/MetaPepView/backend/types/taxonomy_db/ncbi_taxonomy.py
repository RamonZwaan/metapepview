"""This module describes the object that manages the ncbi taxonomy database.
"""

import io
from pathlib import Path
from typing import Dict, List, Tuple, overload, Literal, TypeVar, Sequence, IO
from collections import defaultdict
from functools import lru_cache

from zipfile import ZipFile
from tarfile import TarFile

import pandas as pd
import numpy as np

from metapepview.backend.types.taxonomy_db.taxonomy_database import TaxonomyDatabase
from metapepview.constants import GlobalConstants



class NcbiTaxonomy(TaxonomyDatabase):
    
    FILE_NAMES = GlobalConstants.ncbi_taxonomy_files
    SERVER_URL = GlobalConstants.ncbi_taxonomy_url
    
    NODES_COLUMNS = ("taxonomy_id", "parent_id", "rank", "code", "div_id",
                     "div_flag", "code_id", "gc_flag", "mit_code_id", "mgc_flag",
                     "hid_flag", "root_flag", "comments", "plast_code_id",
                     "pgc_flag", "spec_sp", "hydr_code_id", "hgc_flag")
    NAMES_COLUMNS = ('taxonomy_id', 'name_txt', 'unique_name', 'name_class')
    LINEAGE_COLUMNS = ('taxonomy_id', 'lineage')
    
    RANK_LIST = [x.lower() for x in GlobalConstants.standard_lineage_ranks]
    RANK_DICT = dict(zip(RANK_LIST, range(7)))
    
    ROOT_NAME = 1.0
    
    
    def __init__(self,
                 taxonomy_dict: Dict[int | float, 
                                     Tuple[List[int],
                                           List[int],
                                           str | float,
                                           str | float]],
                 name_dict: Dict[str, List[int]]):
        # key: taxonomy id, value: (standard) lineage + 
        #                          child tax id + 
        #                          rank (if valid) +
        #                          taxonomy name
        self.taxonomy_dict = taxonomy_dict
        self.name_dict = name_dict
    
    
    @classmethod
    def get_root_id(cls) -> float:
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
                      lineages: List[int | float],
                      root_on_empty: bool) -> int | float:
        ...

    @overload
    @classmethod
    def lineage_to_id(cls,
                      lineages: List[List[int | float]],
                      root_on_empty: bool) -> List[int | float]:
        ...
        
    @overload
    @classmethod
    def lineage_to_id(cls,
                      lineages: np.ndarray,
                      root_on_empty: bool) -> np.ndarray | int | float:
        ...
    
    @classmethod
    def lineage_to_id(cls,
                      lineages: List[int | float] | List[List[int | float]] | np.ndarray,
                      root_on_empty: bool = True) -> List[int | float] | np.ndarray | int | float:
        """Convert lineage vector or matrix of lineages to the highest valid taxonomy id for each lineage

        Args:
            lineages (List[int | float] | List[List[int | float]] | np.ndarray): Array or matrix of lineages.
            root_on_empty (bool, Optional): Specify if root taxonomy should be returned if an empty lineage
                is supplied. Else, nan will be returned. Defaults to True.
            root_value (Any, Optional): Root taxonomy id. Will be returned when empty lineage is supplied if
                specified from `root_on_empty`. Defaults to "Root".
                
        Returns:
            List[int | float] | np.ndarray | int | float: Taxonomy id for each lineage
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
                                         (~np.isnan(lineages)).cumsum(1).argmax(1).reshape(-1, 1),
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
        This fills the gaps of lineages for which there is a high rank (e.g. genus, species)
        annotation, but no known annotation for some parent ranks.

        Args:
            lineages (List[int | float] | List[List[int | float]] | np.ndarray | pd.DataFrame | pd.Series):
                A list or matrix of lineages

        Returns:
            List[int | float] | List[List[int | float]] | np.ndarray | pd.DataFrame | pd.Series:
                The lineage list/matrix with filled gaps.
        """
        # variables to specify format of original parameter, used to convert back to expected value
        convert_list = False
        convert_dim = False
        convert_df = False
        convert_series = False
        df_cols, df_index = None, None
        series_name, series_index = None, None

        if isinstance(lineages, pd.DataFrame):
            convert_df = True
            df_cols = lineages.columns
            df_index = lineages.index
            lineages_arr = lineages.to_numpy()    
        elif isinstance(lineages, pd.Series):
            convert_series = True
            series_name = lineages.name
            series_index = lineages.index
            lineages_arr = np.array(lineages.to_list())
        # if list supplied convert to array, remember to convert back to expected type
        elif isinstance(lineages, list):
            convert_list = True
            lineages_arr = np.array(lineages)
        elif isinstance(lineages, np.ndarray):
            lineages_arr = lineages
        
        # vector array should be converted to 2-d array, also here remember to convert back
        if lineages_arr.ndim == 1:
            convert_dim = True
            lineages_arr = lineages_arr.reshape(1, -1)
        
        # Create array that stores last valid taxa for each lineage
        current_tax_id = np.empty(lineages_arr.shape[0])
        current_tax_id[:] = np.nan
        
        # iterate from species to superkingdom, replace nan with valid taxa if present
        for rank_index in range(lineages_arr.shape[1] - 1, -1, -1):
            column = lineages_arr[:, rank_index]  

            notnan = ~np.isnan(column)
            current_tax_id[notnan] = column[notnan]
            
            # fill nan with higher ranked taxa       
            lineages_arr[:, rank_index] = current_tax_id
        
        # convert to same format as input parameter
        if convert_dim == True:
            lineages_arr = lineages_arr.reshape(-1)
        if convert_list == True:
            lineages_arr = lineages_arr.tolist()
        if convert_df == True:
            lineages_arr = pd.DataFrame(lineages_arr, index=df_index, columns=df_cols)
        if convert_series == True:
            lineages_arr = pd.Series(list(lineages_arr), index=series_index, name=series_name)
         
        return lineages_arr # type:ignore (output type does not match TypeVar exactly)
    
    
    @lru_cache(maxsize=None)
    def id_to_standard_lineage(self, tax_id: int | float) -> Tuple[int | float, ...]:
        """Return the lineage from a given taxonomy id, including taxid.

        Args:
            tax_id (int | float): Taxonomy id.

        Returns:
            Tuple[int | float]: Standard lineage array.
        """
        if np.isnan(tax_id) or not self.id_in_dataset(tax_id):
            return (np.nan,)*7
        # convert to int in case float is supplied
        tax_id = int(tax_id)
        
        # get complete lineage up to (and including) tax id
        lineage = self.taxonomy_dict[tax_id][0] + [tax_id]
        
        # if lineage does not go to species, add none values
        lineage_ranks = [self.taxonomy_dict[i][2] for i in lineage]
        
        # build standard lineage (with gaps)
        standard_lineage = [np.nan] * 7
        for lin_index, rank_name in enumerate(lineage_ranks):
            if rank_name in self.RANK_DICT.keys():
                rank_index = self.RANK_DICT[rank_name] # type: ignore
                standard_lineage[rank_index] = lineage[lin_index]
        
        return tuple(standard_lineage)
        
        
    def id_in_dataset(self, tax_id: int | float) -> bool:
        """Check if a given taxonomy id is present in the taxonomy database.

        Args:
            tax_id (str): Taxonomy id.

        Returns:
            bool: Test that id is in dataset.
        """
        if np.isnan(tax_id):
            return False
        else:
            return int(tax_id) in self.taxonomy_dict.keys()
        
        
    def id_to_rank(self, tax_id: int | float) -> str | float:
        """Get the rank of a specified taxonomy id.

        Args:
            tax_id (int | float): Taxonomy id.

        Returns:
            int | float: Taxonomy rank.
        """
        if np.isnan(tax_id) or not self.id_in_dataset(tax_id):
            return np.nan
        
        return self.taxonomy_dict[int(tax_id)][2]
    
    
    def id_to_parent(self,
                     tax_id: int | float,
                     parent_rank: str,
                     absent_rank_behavior: str = "nan") -> int | float:
        """Retrieve taxonomy id at specified parent rank for a given tax id.
        If no id is present at the desired rank for a given taxonomy, the
        method will do the following, based on user settings specified in
        `absent_rank_behavior`:
        
            - take the first lower rank that has a known annotation ("lower").
            - the first higher rank that has a known annotation ("upper").
            - return nan ("nan").

        Args:
            tax_id (int): taxonomy id
            parent_rank (str): Desired rank to retrieve taxonomy id
            absent_rank_behavior (str, optional): Behavior setting if rank
            is absent for lineage. Options: {"upper", "lower", "nan"}. Defaults to "nan".

        Raises:
            ValueError: Invalid rank name supplied.

        Returns:
            int | float: Taxonomy id of parent.
        """
        # check that taxonomy id is valid
        if np.isnan(tax_id):
            return np.nan
            
        elif not self.id_in_dataset(tax_id):
            print(f"""
                Tax ID {int(tax_id)} not recognized in ncbi taxonomy database.
                """)
            return np.nan
        
        # consider potential capitalized
        parent_rank = parent_rank.lower()
        
        # check that given taxonomy rank is valid
        if parent_rank not in self.RANK_LIST:
            raise ValueError(f"Invalid rank name supplied: '{parent_rank}'")
        
        # return the taxonomy name at the specified parent rank, if no annotation, return nan
        parent_rank_index = self.RANK_DICT[parent_rank]
        standard_lineage = self.id_to_standard_lineage(tax_id)
        
        # if nan behavior, return taxid at parent rank, regardless if there is a taxid or not
        if absent_rank_behavior.lower() == "nan":
            return standard_lineage[parent_rank_index]
        elif absent_rank_behavior.lower() == "upper":
            while True:
                parent_taxid = standard_lineage[parent_rank_index]
                if np.isnan(parent_taxid) and parent_rank_index < len(self.RANK_LIST) - 1:
                    parent_rank_index += 1
                    continue
                return parent_taxid
        elif absent_rank_behavior.lower() == "lower":
            while True:
                parent_taxid = standard_lineage[parent_rank_index]
                if np.isnan(parent_taxid) and parent_rank_index > 0:
                    parent_rank_index -= 1
                    continue
                # if complete lineage is empty, return root: no nan since tax id in dataset, therefore valid
                elif np.isnan(parent_taxid):
                    return 1
                return parent_taxid
        else:
            raise ValueError("invalid rank behavior value given.")        


    def id_to_parents(self, tax_id: int | float) -> np.ndarray:
        """Return array of parent taxa for a given taxonomy id.

        Args:
            tax_id (int | float): Taxonomy id.

        Returns:
            np.ndarray: Array of parent taxa.
        """
        if np.isnan(tax_id) or not self.id_in_dataset(tax_id):
            return np.array([])
        
        # return array of lineage including taxonomy id
        return np.array(self.taxonomy_dict[tax_id][0] + [tax_id])
    
    
    def id_to_children(self, tax_id: int | float) -> np.ndarray:
        """Return array of all children that belong to a supplied taxonomy id.

        Args:
            tax_id (str | float): Taxonomy id.

        Returns:
            np.ndarray: Array of child taxonomy id's under a given tax_id.
        """
        if np.isnan(tax_id) or not self.id_in_dataset(tax_id):
            return np.array([])
        
        # store all children in list
        children = [tax_id]
        current_iter = [tax_id]
        
        while len(current_iter) > 0:
            new_iter = []
            # fetch all offspring in list
            for i in current_iter:
                new_iter += self.taxonomy_dict[i][1] 
            # update current iteration id's with offspring
            current_iter = new_iter
            # add offspring to total children list
            children += new_iter
        
        return np.array(children)
        
        
    def taxa_to_lca(self,
                    tax_ids: Sequence[str | float] | pd.Series | np.ndarray,
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
            
        # remove nan values
        tax_ids = tax_ids[~np.isnan(tax_ids)]
        
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
        lineages = [self.taxonomy_dict[tax_id][0] + [tax_id] for tax_id in tax_ids]

        # fetch highest rank id for which there is consensus among lineages
        return self.lineages_to_lca(lineages)
    
    
    def id_to_name(self, tax_id: int | float) -> str | float:
        """Convert taxonomy id to taxa name.

        Args:
            tax_id (int | float): Taxonomy id

        Returns:
            str | float:  Taxa name
        """
        if np.isnan(tax_id) or not self.id_in_dataset(tax_id):
            return np.nan
        else:
            return self.taxonomy_dict[int(tax_id)][3]
        
        
    def lineage_id_to_name(self, id_lineage: List[int | float] | Tuple[int | float, ...]) -> Tuple[str | float, ...]:
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
                   print_fails: bool) -> int | float | List[int]:
        ...

    @overload
    def name_to_id(self,
                   tax_name: str | float,
                   on_duplicates: Literal["nan"],
                   print_fails: bool) -> int | float:
        ...
    
    def name_to_id(self,
                   tax_name: str | float,
                   on_duplicates: str="nan",
                   print_fails: bool=False) -> int | List[int] | float:
        """Convert organism name to taxonomy id.
        
        Args:
            tax_name (str | float): Taxonomy name.
            on_duplicates (str, Optional): Specify function behavior if multiple
                id's encountered with same name. Options: {"nan", "all"}.
                'nan' will return nan when multiple id's encountered,
                'all' will return list of taxonomy id's.
            print_fails (bool, optional): Print cases where tax id retrieval
                fails to stdout. Defaults to False.

        Raises:
            ValueError: Invalid `on_duplicates` value given.

        Returns:
            int | List[int] | float: Taxonomy id('s) coupled to name
        """
        if isinstance(tax_name, float):
            if print_fails is True:
                print("Empty input encountered...")
            return np.nan
        
        if tax_name not in self.name_dict.keys():
            if print_fails is True:
                print(f"No valid id's for {tax_name}")
            return np.nan
    
        tax_ids = self.name_dict[tax_name]  # type: ignore
        
        # return id's or nan based on selected parameter options
        if len(tax_ids) == 1:
            return tax_ids[0]
        elif on_duplicates.lower() == "nan":
            if print_fails is True:
                print(f"Multiple id's for {tax_name}: {tax_ids}")
            return np.nan
        elif on_duplicates.lower() == "all":
            return tax_ids
        else:
            raise ValueError(f"Invalid `on_duplicates` value given '{on_duplicates}'. Choose from {'nan', 'all'}.")
        
        
    
    @classmethod
    def lineages_to_lca(cls, lineages: Sequence[Sequence[int | float]]) -> int | float:
        """Retrieve the last common ancestor tax id from a list of taxonomy
        lineages. This method parses an array of lineage strings retrieved from
        the lineage dataset. It will return the highest ranked taxonomy for which
        there is consensus.

        Args:
            lineages (List[List[int | float] | Tuple[int | float]]): List of lineages.

        Returns:
            int | float: Last common ancestor of lineages.
        """
        # return nan if 0 lineages are given as input
        if len(lineages) == 0:
            return np.nan
        
        lca = cls.ROOT_NAME             # start with root
        min_len = min(len(x) for x in lineages)
        
        # parse lineages and update lca as long as there is full consensus
        for i in range(min_len):
            # if gaps encountered in lineage, continue to next rank by definition
            if any([np.isnan(x[i]) for x in lineages]):
                continue
            if len({x[i] for x in lineages}) == 1:
                lca = lineages[0][i]
        
        return lca
        
        
    ################################################################################ 
    # Constructor methods
    ################################################################################ 

    @classmethod
    def from_dmp_archive(cls, dmp_arch: str | Path) -> "NcbiTaxonomy":
        """Import ncbi taxonomy datasets from taxdump archive file.
        Extracts 'nodes.dmp', 'names.dmp', and 'taxidlineage.dmp' files and
        imports them.

        Args:
            dmp_arch (str | Path): Location of taxdump archive in ".zip" or 
                ".tar.gz" format.

        Returns:
            NcbiTaxonomy: NcbiTaxonomy object.
        """
        # convert to Path if str supplied
        dmp_arch = Path(dmp_arch)
        # check for presence archive
        if not dmp_arch.exists():
            raise FileNotFoundError(f"{dmp_arch.as_posix()} not present in supplied directory")
        
        # load archive file
        dmp_arch.suffixes
        if len(dmp_arch.suffixes) == 0:
            raise ValueError(f"Invalid archive format provided: '{dmp_arch.name}'")
        elif dmp_arch.suffixes[-1] == ".zip":
            arch_data = ZipFile(dmp_arch)
        elif dmp_arch.suffixes[-1] == ".tar":
            arch_data = TarFile.open(dmp_arch)
        elif dmp_arch.suffixes[-2:] == [".tar", ".gz"]:
            arch_data = TarFile.open(dmp_arch)
        else:
            raise ValueError(f"Invalid archive format provided: '{dmp_arch.name}'")

        # Encapsulate divergent interfaces to extract single file into buffer
        def extract_file(archive: TarFile | ZipFile, 
                         file_name: str) -> IO[str]:
            if isinstance(archive, TarFile):
                file_data = archive.extractfile(file_name)
                if file_data is None:
                    raise ValueError("Member not found in archive...")
            else:
                file_data = archive.open(file_name)
            return io.TextIOWrapper(file_data, encoding="utf-8")
        
        # initialize data structures to store taxonomy dataset in
        taxonomy_dict, name_dict = cls.__init_object_attribute_data()
        
        nodes_file, names_file, lineage_file = cls.FILE_NAMES
        with extract_file(arch_data, nodes_file) as nodes_file_data:
            taxonomy_dict = cls.__import_nodes(nodes_file_data, taxonomy_dict)
        
        with extract_file(arch_data, names_file) as names_file_data:
            taxonomy_dict, name_dict = cls.__import_names(names_file_data, 
                                                          taxonomy_dict, 
                                                          name_dict)
        
        with extract_file(arch_data, lineage_file) as lineage_file_data:
            taxonomy_dict = cls.__import_lineage(lineage_file_data, 
                                                 taxonomy_dict)
            
        return cls(dict(taxonomy_dict), dict(name_dict))
    
    
    @classmethod
    def from_dmp_folder(cls, dmp_dir: str | Path) -> "NcbiTaxonomy":
        """Import ncbi taxonomy datasets from folder.
        imports 'nodes.dmp', 'names.dmp', and 'taxidlineage.dmp' files.

        Args:
            dmp_dir (str | Path): Location folder.

        Raises:
            FileNotFoundError: Database files not present in location.

        Returns:
            NcbiTaxonomy: NcbiTaxonomy object.
        """
        # convert to Path if str supplied
        dmp_dir = Path(dmp_dir)
        
        # check for presence of all files
        nodes_file, names_file, lineage_file = [Path(dmp_dir, x) for x in cls.FILE_NAMES]
        if all([x.exists() for x in [nodes_file, names_file, lineage_file]]):
            # call file constructor to import files
            return cls.from_dmp_files(nodes_file, names_file, lineage_file)
        else:
            raise FileNotFoundError(f"one of {cls.FILE_NAMES} not present in supplied directory")

        
    @classmethod
    def from_dmp_files(cls,
                       nodes_file: str | Path,
                       names_file: str | Path,
                       lineage_file: str | Path) -> "NcbiTaxonomy":
        """Import ncbi taxonomy database from dataset files. NcbiTaxonomy
        imports 'nodes.dmp', 'names.dmp', and 'taxidlineage.dmp' files.

        Args:
            nodes_file (str | Path): Location 'nodes.dmp' file.
            names_file (str | Path): Location 'names.dmp' file.
            lineage_file (str | Path): Location 'taxidlineage.dmp' file.

        Returns:
            NcbiTaxonomy: NcbiTaxonomy object.
        """
        # check if path location is valid
        if not all([x.exists() for x in [Path(nodes_file),
                                         Path(names_file),
                                         Path(lineage_file)]]):
            raise FileNotFoundError("Provided files missing...")
        
        # initialize data structures to store taxonomy dataset in
        taxonomy_dict, name_dict = cls.__init_object_attribute_data()
        
        with Path(nodes_file).open() as nodes_file_data:
            taxonomy_dict = cls.__import_nodes(nodes_file_data, taxonomy_dict)
        with Path(names_file).open() as names_file_data:
            taxonomy_dict, name_dict = cls.__import_names(names_file_data, 
                                                          taxonomy_dict, 
                                                          name_dict)
        with Path(lineage_file).open() as lineage_file_data:
            taxonomy_dict = cls.__import_lineage(lineage_file_data, 
                                                 taxonomy_dict)
            
        return cls(dict(taxonomy_dict), dict(name_dict))
    
    @staticmethod
    def __init_object_attribute_data():
        """Constructor to initialize data structures used for storage of the
        taxonomy dataset. It generates a taxonomy dict and a name dict.
        """
        # data will be stored as dictionaries
        # taxonomy dict: {tax_id: lineage, child ids, rank, tax name}
        taxonomy_dict = defaultdict(lambda: [(), [], np.nan, np.nan])
        # name dict: {tax_name: list[tax_id]}
        name_dict = defaultdict(lambda: [])

        return (taxonomy_dict, name_dict)

    @staticmethod
    def __import_nodes(read_file: IO[str], taxonomy_dict: Dict) -> Dict:
        for line in read_file.readlines():
            # split row in separate cells
            line = line.replace("\t|\n", "") # remove trailing |
            line = line.split("\t|\t")
            
            # set rank name of taxid in dict
            taxonomy_dict[int(line[0])][2] = line[2]
            
            # add taxid as offspring to parent tax in dict, ignore for root
            if line[1] == line[0]:
                continue
            taxonomy_dict[int(line[1])][1].append(int(line[0]))
        
        # return dictionary
        return taxonomy_dict

    @staticmethod
    def __import_names(read_file: IO[str], 
                       taxonomy_dict: Dict, 
                       name_dict: Dict) -> Tuple[Dict, Dict]:
        for line in read_file.readlines():
            # split row in separate cells
            line = line.replace("\t|\n", "") # remove trailing |
            line = line.split("\t|\t")
            
            # update dict with scientific name of tax id, ignore other names
            if line[3] == "scientific name":
                # set name of taxid in dict
                taxonomy_dict[int(line[0])][3] = line[1]
                name_dict[line[1]].append(int(line[0]))
            # for non-scientific names, only add name-to-id link, might be overwritten
            else:
                if line[1] not in name_dict.keys():
                    name_dict[line[1]].append(int(line[0]))
                    
        # return dictionary
        return (taxonomy_dict, name_dict)

    @staticmethod
    def __import_lineage(read_file: IO[str], 
                         taxonomy_dict: Dict) -> Dict:
        for line in read_file.readlines():
            # split row in separate cells
            line = line.replace("\t|\n", "") # remove trailing |
            line = line.split("\t|\t")
            
            # convert lineage string to list of integers
            lineage = line[1]
            # if non-empty lineage, remove trailing whitespace
            if lineage == "":
                lineage = []
            else:
                lineage = lineage[:-1]
                lineage = lineage.split(" ")
                lineage = [int(i) for i in lineage]
            
            # set lineage to taxonomy id in dict
            taxonomy_dict[int(line[0])][0] = lineage
        
        # return dictionary
        return taxonomy_dict


    # private methods
    @staticmethod
    def __import_dmp(dmp_file: Path, columns: Tuple[str, ...], encoder: str="utf8") -> pd.DataFrame:
        """Parse dmp file formats and store into dataframe.

        Args:
            dmp_file (Path): Location of dmp file
            columns (Tuple[str, ...]): Columns of file contents.
            encoder (str, optional): Encoder to read dmp file with. Defaults to "utf8".

        Raises:
            FileNotFoundError: File does not exist

        Returns:
            pd.DataFrame: Imported dataset.
        """
        # check if path location is valid
        if not dmp_file.exists():
            raise FileNotFoundError("Path to dmp file does not exist.")
            
        data = list()
        
        # parse dmp file
        with open(dmp_file, 'r', encoding=encoder) as read_file:
            for line in read_file.readlines():
                line = line.replace("\t|\n", "") # remove trailing |
                line = line.split("\t|\t")
                data.append(line)
        
        return pd.DataFrame(data, columns=columns)
