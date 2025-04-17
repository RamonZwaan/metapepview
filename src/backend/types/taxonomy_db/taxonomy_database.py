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
    def name_to_id(self,
                   tax_name: str,
                   **kwargs) -> Any:
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
