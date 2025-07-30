from typing import List, Protocol, Self, Type, TypeVar, Dict, ClassVar, IO
import numpy as np
import pandas as pd

from .base_classes import DataValidator
from constants import GlobalConstants


T = TypeVar('T', bound='FunctionDbMapper')

class FunctionDbMapper(Protocol):
    NANFILL: ClassVar[str | float]
    LIST_DELIM: ClassVar[str]

    idx_dict: Dict[str, int]
    df: pd.DataFrame
    
    @classmethod
    def read_file_buffer(cls: Type[T],
                         file_buffer: IO[str],
                         max_evalue: float=1e-6) -> T:...
    
    @classmethod
    def from_dataframe(cls: Type[T], input_df: pd.DataFrame) -> T:...
    
    def data_from_accession(self, accession_list: str | List[str]) -> pd.DataFrame | float:
        """Retrieve rows from eggnog dataset related to the supplied accession list.

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
        
        # data_list = []
        # retrieve data from dict
        # for acc in valid_accessions:
            # data_list.append(self.idx_dict[acc])
        
        # fetch all data from accessions
        output = self.df.iloc[[self.idx_dict[i] for i in valid_accessions]]
        
        # return rows from dataset
        # return pd.concat(data_list, axis=1).T
        return output



class EggnogMapper(FunctionDbMapper, DataValidator):
    
    REQUIRED_FIELDS: List[str] = ['query', 'evalue', 'eggNOG_OGs', 'COG_category',
                                  'Preferred_name', 'EC', 'KEGG_ko', 'CAZy']
    NUMERIC_FIELDS: List[str] = ['evalue']
    
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
                         max_evalue: float=1e-6) -> Self:
        # read data into dataframe format
        df = pd.read_csv(file_buffer,
                        sep='\t', header=4)
        
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
        
        return cls.from_dataframe(df)
    
    
    @classmethod
    def from_dataframe(cls, input_df: pd.DataFrame) -> Self:
        success, msg = cls.validate_input(input_df)
        if success is False:
            raise ValueError(msg)
        
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
                         query_col: int = 0,
                         ko_col: int = 1,
                         **kwargs) -> Self:
        kegg_df = pd.read_csv(file_buffer,
                              sep='\t',
                              names=["query", "KEGG_ko"])
        return cls.from_dataframe(kegg_df)
    
    
    @classmethod
    def from_dataframe(cls, input_df: pd.DataFrame) -> Self:
        success, msg = cls.validate_input(input_df)
        if success is False:
            raise ValueError(msg)
        
        # drop rows without accession and without
        input_df.dropna(subset=['query', 'KEGG_ko'], inplace=True)
        # input_df.dropna(subset='KEGG_ko', inplace=True)
        
        # set query as index, while keeping query column in dataset
        input_df = input_df.set_index('query', drop=False)
        
        #acc_dict = dict(input_df.iterrows())
        acc_dict = {j:i for i, j in enumerate(input_df.index)}


        return cls(acc_dict, input_df)    

