from typing import IO, Tuple, Self
import pandas as pd

from .base_class import AccessionTaxaMap
from backend.types.definitions import TaxonomyElementFormat
from backend.types.taxonomy_db import NcbiTaxonomy, GtdbGenomeToNcbi


class AccessionTaxaMapNcbi(AccessionTaxaMap):
    
    TAXONOMY_FORMAT = "NCBI"
    ROOT_TAXONOMY_ID = 1
    
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
                           taxonomy_obj: NcbiTaxonomy | None = None) -> Self:
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
            taxonomy_obj (NcbiTaxonomy | None, optional): NcbiTaxonomy
                object for LCA processing of redundant protein accessions.
                Defaults to None.

        Returns:
            Self: AccessionTaxaMapNcbi instance.
        """
        
        # it is assumed that the dataset starts with accession, with second column tax id
        prot_df = pd.read_csv(str_file_obj,
                              usecols=[acc_col, tax_col], 
                              names=["accession", "taxonomy_id"],
                              sep=delimiter)
        
        return cls.__from_dataframe(
            prot_df,
            acc_regex,
            drop_duplicates,
            tax_name_to_id,
            taxonomy_obj
        )
    
    @classmethod
    def from_gtdb_genome_ids(cls,
                             genome_to_ncbi_map: GtdbGenomeToNcbi,
                             str_file_obj: IO[str],
                             acc_col: int = 0,
                             tax_col: int = 1,
                             acc_regex: str | None = None,
                             delimiter: str | None = ",",
                             drop_duplicates: bool = True,
                             tax_name_to_id: bool = False,
                             taxonomy_obj: NcbiTaxonomy | None = None) -> Self:
        """Create AccessionTaxaMapNcbi object from GTDB genome mapping file. 
        It maps Genbank/Refseq genome id's to their corresponding NCBI taxonomy 
        id's.
        
        Note:
            Only GTDB genome id's can be mapped to NCBI taxonomy id's, GTDB 
            taxonomy id's cannot be matched to NCBI taxa as their tree structures
            differ. Therefore, any non-genome id in the dataset is discarded.

        Args:
            genome_to_ncbi_map (GtdbGenomeToNcbi): Object that maps genome id
                to NCBI taxonomy id using the metadata file from GTDB.
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
            taxonomy_obj (NcbiTaxonomy | None, optional): NcbiTaxonomy
                object for LCA processing of redundant protein accessions.
                Defaults to None.

        Returns:
            Self: Instance of AccessionTaxaMapNcbi.
        """
        # it is assumed that the dataset starts with accession, with second column tax id
        prot_df = pd.read_csv(str_file_obj,
                              usecols=[acc_col, tax_col], 
                              names=["accession", "taxonomy_id"],
                              sep=delimiter)
        
        prot_df["taxonomy_id"] = prot_df["taxonomy_id"].apply(
            genome_to_ncbi_map.genome_to_ncbi
        )
        prot_df = prot_df.dropna(subset="taxonomy_id")
        
        return cls.__from_dataframe(
            prot_df,
            acc_regex,
            drop_duplicates,
            tax_name_to_id,
            taxonomy_obj
        )

    @classmethod
    def __from_dataframe(cls,
                         prot_df: pd.DataFrame,
                         acc_regex: str | None=None,
                         drop_duplicates: bool = True,
                         tax_name_to_id: bool = False,
                         taxonomy_obj: NcbiTaxonomy | None = None) -> Self:
        """Generate class instance from mapping data loaded in dataframe.
        Private function that processes data from the public constructor
        functions.

        Args:
            prot_df (pd.DataFrame): Mapping data as dataframe.
            acc_regex (str | None, optional): Regex pattern to process accession
                column, storing only the pattern match portion. Defaults to None.
            drop_duplicates (bool, optional): In case of duplicate protein
                accessions, keep only first occurrence. If false, duplicate
                accessions are merged by taking the last common ancestor for all
                taxa within the protein accession group. For this operation a
                TaxonomyDatabase object needs to be added. Defaults to True.
            tax_name_to_id (bool, optional): When elements in taxonomy column
                are names instead of id's, convert values to id. Defaults to
                False.
            taxonomy_obj (NcbiTaxonomy | None, optional): NcbiTaxonomy
                object for LCA processing of redundant protein accessions.
                Defaults to None.

        Returns:
            Self: AccessionTaxaMapNcbi instance.
        """
        # if tax names in column, convert to id
        if tax_name_to_id is True and taxonomy_obj is not None:
            print("Taxonomy name to id conversion\n\nFailed conversions:")  
            prot_df.loc[:, "taxonomy_id"] = prot_df["taxonomy_id"]\
                .apply(taxonomy_obj.name_to_id, print_fails=True)\
                .apply(pd.to_numeric)
        elif tax_name_to_id is True:
            raise ValueError("Taxonomy database required for name to id conversion...")
        
        if not pd.api.types.is_numeric_dtype(prot_df["taxonomy_id"]):
            raise ValueError("Invalid taxonomy id's encountered. Ensure that the complete tax id column is numeric.")
        
        # apply pattern on accession if given        
        if acc_regex is not None:
            repl = lambda m: m.group(0) if m is not None else ""
            prot_df["accession"] = prot_df["accession"].str.replace(acc_regex, repl, regex=True)
        
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
                       acc_col: int=0,
                       tax_col: int=1,
                       delimiter: str=", ",
                       element_format: TaxonomyElementFormat = "taxonomy id") -> Tuple[bool, str]:

        # check correct format for ncbi taxonomy
        def check_numeric(prot_df: pd.DataFrame) -> Tuple[bool, str | None]:
            if not pd.api.types.is_numeric_dtype(prot_df.iloc[:, tax_col]):
                err_msg = "Invalid taxonomy format, expects ncbi id's (numeric)"
                return (False, err_msg)
            else:
                return (True, None)
            
        # only include numeric check for ncbi id's
        if element_format == "taxonomy id":
            validator_funcs = [check_numeric]
        else:
            validator_funcs = []
        # perform regular validation
        success, err_msg = AccessionTaxaMap.validate_input(
            str_file_obj,
            acc_col,
            tax_col,
            delimiter,
            validator_funcs)
        
        return (success, err_msg)
