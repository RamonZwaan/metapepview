from typing import List, Protocol, Self, Sequence, Any, Type, TypeVar, Dict, ClassVar, IO, Tuple
import numpy as np
import pandas as pd

from metapepview.backend.types.base_classes import DataValidator
from metapepview.backend.types.taxonomy_db import TaxonomyDatabase
from metapepview.backend.types.taxonomy_map.base_class import AccessionTaxaMapMethods
from metapepview.constants import GlobalConstants
from metapepview.backend.utils import regex_over_column


T = TypeVar('T', bound='FunctionDbMapper')

class FunctionDbMapper(Protocol):
    NANFILL: ClassVar[str | float]
    LIST_DELIM: ClassVar[str]

    idx_dict: Dict[str, int]
    df: pd.DataFrame
    
    @classmethod
    def read_file_buffer(cls: Type[T],
                         file_buffer: IO[str],
                         max_evalue: float=1e-6,
                         acc_regex: str | None=None) -> T:
         ...
    
    @classmethod
    def from_dataframe(cls: Type[T], 
                       input_df: pd.DataFrame,
                       acc_regex: str | None=None) -> T:
        ...
    
    def data_from_accession(self, accession_list: str | List[str]) -> pd.DataFrame | float:
        """Retrieve rows from dataset related to the supplied accession list.

        Args:
            accession_list (str | List[str]): list of protein accession

        Returns:
            pd.DataFrame: data from dataset that have these accessions
        """
        if isinstance(accession_list, str):
            accession_list = [accession_list]
            
        # filter accessions not in dataset
        valid_accessions = []
        for acc in accession_list:
            if acc in self.idx_dict.keys():
                valid_accessions.append(acc)
        
        # when no valid data in list, return nan
        if len(valid_accessions) == 0:
            return np.nan
        
        # fetch all data from accessions
        output = self.df.iloc[[self.idx_dict[i] for i in valid_accessions]]
        
        # return rows from dataset
        return output



class EggnogMapper(FunctionDbMapper, DataValidator):
    
    REQUIRED_FIELDS: List[str] = ['query', 'evalue', 'eggNOG_OGs', 'COG_category',
                                  'Preferred_name', 'EC', 'KEGG_ko', 'CAZy']
    NUMERIC_FIELDS: List[str] = ['evalue']
    
    NANFILL = GlobalConstants.func_db_combine_nan_fill
    FIELD_DELIM = "\t"
    LIST_DELIM = GlobalConstants.func_db_combine_delimiter
    
    def __init__(self,
                 idx_dict: Dict,
                 df: pd.DataFrame) -> None:
        # map accession to row index for quick data retrieval
        self.idx_dict = idx_dict
        self.df = df
    
    
    @classmethod
    def read_file_buffer(cls,
                         file_buffer: IO[str],
                         max_evalue: float=1e-6,
                         acc_regex: str | None=None) -> Self:
        
        # assert presence of metadata fields
        while True:
            line_txt = file_buffer.readline()
            if not line_txt:
                raise ValueError("No data stored; end of file reached...")
            if not line_txt.startswith("##"):
                header = line_txt.split(cls.FIELD_DELIM)
                break

        # read data into dataframe format
        df = pd.read_csv(file_buffer,
                         names=header, 
                         sep=cls.FIELD_DELIM,
                         engine="python")
        
        # rename columns and set index to accession
        df.rename(columns={"#query": "query"}, inplace=True)
        #df.set_index("query", drop=True, inplace=True)

        df = df[df["evalue"] < max_evalue]

        # convert '-' to NaN
        df.replace("-", np.nan, inplace=True)

        # remove prefix from KEGG KO annotations
        df["KEGG_ko"] = df["KEGG_ko"].str.replace(r"ko:", "")
        
        # filter columns away to conserve space
        df = df[['query', 'evalue', 'eggNOG_OGs', 'COG_category', 'Preferred_name', 'EC', 'KEGG_ko', 'CAZy']]
        
        return cls.from_dataframe(df, acc_regex)
    
    
    @classmethod
    def from_dataframe(cls, 
                       input_df: pd.DataFrame,
                       acc_regex: str | None=None) -> Self:
        success, msg = cls.validate_input(input_df)
        if success is False:
            raise ValueError(msg)
        
        # If regex provided, extract substring from each query using regex
        # if match, get first occurence of match
        if acc_regex != "" and acc_regex is not None:
            input_df.loc[:, 'query'] = regex_over_column(input_df['query'],
                                                         acc_regex,
                                                         no_match_to_nan=False)

        # drop rows without accession
        input_df.dropna(subset='query', inplace=True)

        # set query as index, while keeping query column in dataset
        input_df = input_df.set_index('query', drop=False)
        
        #acc_dict = dict(input_df.iterrows())
        # couple df index values to numeric index
        acc_dict = {j:i for i, j in enumerate(input_df.index)}
        
        return cls(acc_dict, input_df)
    

class KeggMapper(FunctionDbMapper, DataValidator):
    REQUIRED_FIELDS: List[str] = ['query', 'KEGG_ko']
    NUMERIC_FIELDS: List[str] = []
    
    NANFILL = GlobalConstants.func_db_combine_nan_fill
    LIST_DELIM = GlobalConstants.func_db_combine_delimiter    

    def __init__(self,
                 idx_dict: Dict,
                 df: pd.DataFrame) -> None:
        # map accession to row index for quick data retrieval
        self.idx_dict = idx_dict
        self.df = df   
    

    @classmethod
    def read_file_buffer(cls,
                         file_buffer: IO[str],
                         *args,
                         acc_regex: str | None=None,
                         query_col: int = 0,
                         ko_col: int = 1,
                         **kwargs) -> Self:
        kegg_df = pd.read_csv(file_buffer,
                              sep='\t',
                              engine="python",
                              names=["query", "KEGG_ko"])
        return cls.from_dataframe(kegg_df, acc_regex)
    
    
    @classmethod
    def from_dataframe(cls, 
                       input_df: pd.DataFrame,
                       acc_regex: str | None=None) -> Self:
        success, msg = cls.validate_input(input_df)
        if success is False:
            raise ValueError(msg)

        # If regex provided, extract substring from each query using regex
        # if match, get first occurence of match
        if acc_regex != "" and acc_regex is not None:
            input_df.loc[:, 'query'] = regex_over_column(input_df['query'],
                                                         acc_regex,
                                                         no_match_to_nan=False)

        # drop rows without accession and without kegg ko
        input_df = input_df.dropna(subset=['query', 'KEGG_ko'])
        
        # set query as index, while keeping query column in dataset
        input_df = input_df.set_index('query', drop=False)
        
        #acc_dict = dict(input_df.iterrows())
        acc_dict = {j:i for i, j in enumerate(input_df.index)}


        return cls(acc_dict, input_df)    


class GhostkoalaMapper(FunctionDbMapper, DataValidator, AccessionTaxaMapMethods):
    REQUIRED_FIELDS: List[str] = ['query', 'KEGG_ko', 'Genus', 'GHOSTX score']
    NUMERIC_FIELDS: List[str] = ['GHOSTX score']
    
    NANFILL = GlobalConstants.func_db_combine_nan_fill
    LIST_DELIM = GlobalConstants.func_db_combine_delimiter    

    FILE_DELIM = '\t'
    FILE_COLS = ["query", 
                 "KEGG_ko", 
                 "Kingdom", 
                 "Phylum",
                 "Genus",
                 "Genes ID",
                 "GHOSTX score"]
    ACC_PREFIX = "user:"

    def __init__(self,
                 idx_dict: Dict,
                 func_df: pd.DataFrame,
                 tax_map: Dict[str, str]) -> None:
        # map accession to row index for quick data retrieval
        self.idx_dict = idx_dict
        self.df = func_df 
        self.accession_tax_dict = tax_map  


    @classmethod
    def validate_buffer(cls, file_buffer: IO[str]) -> IO[str]:
        # validate format by parsing first 10 lines
        for _ in range(0, 10):
            line = file_buffer.readline()
            ncols = line.count(cls.FILE_DELIM) + 1
            ko_prefix = line.startswith(cls.ACC_PREFIX)
            if ncols != len(cls.FILE_COLS):
                raise ValueError(f"Invalid GhostKOALA file format, expected {len(cls.FILE_COLS)} fields, got {ncols}...")
            if not ko_prefix:
                raise ValueError("Invalid accession format encountered, expected 'user:' prefix")
        file_buffer.seek(0)
        return file_buffer


    @classmethod
    def read_file_buffer(cls,
                         file_buffer: IO[str],
                         acc_regex: str | None=None,
                         **kwargs) -> Self:
        # perform validation of file buffer, required as the file has no headers
        file_buffer = cls.validate_buffer(file_buffer)

        kegg_df = pd.read_csv(file_buffer,
                              sep='\t',
                              engine="python",
                              names=cls.FILE_COLS)
        return cls.from_dataframe(kegg_df, acc_regex)
    
    
    @classmethod
    def from_dataframe(cls, 
                       input_df: pd.DataFrame,
                       acc_regex: str | None=None) -> Self:
        success, msg = cls.validate_input(input_df)
        if success is False:
            raise ValueError(msg)

        # remove the prefix that KEGG adds to accessions in the file
        input_df.loc[:, 'query'] = input_df['query'].str.removeprefix("user:")

        # If regex provided, extract substring from each query using regex
        # if match, get first occurence of match
        if acc_regex != "" and acc_regex is not None:
            input_df.loc[:, 'query'] = regex_over_column(input_df['query'],
                                                         acc_regex,
                                                         no_match_to_nan=False)

        # drop rows without accession
        input_df = input_df.dropna(subset=['query'])
        
        # set query as index, while keeping query column in dataset
        input_df = input_df.set_index('query', drop=False)
        
        # create function dataset containing only function information for accessions
        func_df = input_df[["query", "KEGG_ko", "Genes ID", "GHOSTX score"]]\
            .dropna(subset="KEGG_ko")

        acc_dict = {j:i for i, j in enumerate(func_df.index)}

        tax_map = input_df["Genus"].to_dict()

        return cls(acc_dict, input_df, tax_map)    


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

        tax_names = [self.accession_to_taxonomy(i) for i in accessions]

        tax_ids = [taxonomy_db.name_to_id(i) for i in tax_names]
        
        return taxonomy_db.taxa_to_lca(tax_ids)
    
    
    def accession_to_taxonomy(self,
                              accession: str) -> int | str | float:
        """Return the taxonomy annotation for a given protein accession.

        Args:
            accession (str): Protein accession.

        Returns:
            int | str | float: Taxonomy annotation from protein accession.
        """
        return self.accession_tax_dict.get(accession, np.nan)

