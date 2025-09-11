from dataclasses import dataclass
from typing import List, Dict, Any, Sequence, IO
from pathlib import Path
import pandas as pd
import numpy as np

from metapepview.backend import *
from metapepview.backend.annotation.data_wrangling import *
# add importer objects into scope
from metapepview.backend.type_operations import import_func_map
from metapepview.backend.io import fasta_to_peptides
from metapepview.backend.utils import *
from metapepview.backend.type_operations import *
from metapepview.backend.annotation.unipept_search import global_taxonomic_annotation
from metapepview.constants import *
from metapepview.backend.exceptions import AnnotationError

# add type classes
from metapepview.backend.types import *


@dataclass
class AnnotationOptions():
    db_search_format: DbSearchSource
    db_search_acc_pattern: str | None
    min_db_search_score: int | float
    db_search_filter_crap: bool
    de_novo_format: DeNovoSource
    min_de_novo_score: int | float
    de_novo_filter_crap: bool
    min_pept_length: int
    tax_db_delimiter: str
    tax_db_name: str | None
    tax_db_accession_format: Literal["Sequence", "Accession"]
    tax_annot_format: TaxonomyMapFormat | None
    tax_db_format: TaxonomyDbFormat | None
    tax_db_element_format: TaxonomyElementFormat | None
    gtdb_to_ncbi: bool
    ncbi_db_loc: str
    global_taxonomy_annotation: bool
    func_db_name: str | None
    func_db_format: FuncAnnotFormat | None
    tax_db_acc_pattern: str
    tax_db_acc_idx: int
    tax_db_tax_idx: int
    func_annot_combine: bool
    func_annot_pattern: str | None
    merge_psms: bool


def annotate_peptides(sample_name: str,
                      db_search_list: Sequence[str | IO[str]],
                      db_search_names: List[str],
                      de_novo_list: List[str | IO[str]],
                      taxonomy_map: str | None,
                      func_annot_map: str | None,
                      tax_db_loc: str,
                      options: AnnotationOptions) -> MetaPepTable:
    """Process db search and de novo data from single or multiple samples into 
    taxonomic and functional annotated peptides dataset. Input metaproteomics
    data are imported, wrangled and annotated with taxonomy and function, before
    being concatenated into a single peptides dataset.

    Args:
        sample_name (str): Name of sample (only if data is merged).
        psm_list (Sequence[str | IO[str]]): List of db search psm files.
        psm_names (List[str]): Names of db search psm files
        de_novo_list (Dict[str, IO[str]] | None): Dictionary of raw MS
            file names and de novo data processed from that run.
        taxonomy_db (str | None): Taxonomy protein mapping database.
        func_annot_map (str | None): Protein to function mapping database.
        ncbi_loc (str): File location of ncbi taxonomy database.
        options (AnnotationOptions): Dataclass object with annotation options

    Returns:
        pd.DataFrame: Peptides dataset.
    """
    # set db format to ncbi if no taxonomy map added
    options.tax_db_format = "NCBI" if options.tax_db_format is None else options.tax_db_format

    if options.tax_db_format == "NCBI":
        tax_db_loc = Path(tax_db_loc, GlobalConstants.ncbi_taxonomy_archive)

    # Convert annotation in GTDB format to NCBI format
    if options.gtdb_to_ncbi is True and \
        taxonomy_map is not None and \
        options.tax_db_format == "GTDB" and \
        options.tax_db_name is not None:
        print("import taxonomy database...")
        try:
            taxonomy_db = import_taxonomy_db(Path(options.ncbi_db_loc,
                                             GlobalConstants.ncbi_taxonomy_archive),
                                            "NCBI")
        except Exception as err:
            raise AnnotationError(f"Failed to import taxonomy database: {err}")
        
        print("import protein to taxonomy map...")
        tax_db_archive_format = determine_archive_format(options.tax_db_name)
        
        # load gtdb to ncbi mapping object
        bac_metadata, ar_metadata = [Path(tax_db_loc, x) for x in GlobalConstants.gtdb_metadata_files]
        try:
            gtdb_to_ncbi_obj = GtdbGenomeToNcbi.from_metadata_file(bac_metadata, 
                                                                   ar_metadata)
        except Exception as err:
            raise AnnotationError(f"Failed to load GTDB to NCBI mapping: {err}")
        
        try:
            tax_df = import_acc_tax_map(taxonomy_map,
                                        options.tax_db_acc_idx,
                                        options.tax_db_tax_idx,
                                        options.tax_db_acc_pattern,
                                        options.tax_db_delimiter,
                                        taxonomy_db,
                                        options.tax_annot_format,
                                        options.tax_db_element_format,
                                        gtdb_to_ncbi_obj,
                                        tax_db_archive_format,
                                        options.tax_db_accession_format == "Sequence")
        except Exception as err:
            raise AnnotationError(f"Failed to import taxonomy annotations: {err}")
        # Change taxonomy db format after conversion to NCBI
        options.tax_db_format = "NCBI"
        
    elif taxonomy_map is not None and \
        options.tax_db_format is not None and \
        options.tax_db_name is not None:
        print("import taxonomy database...")
        try:
            taxonomy_db = import_taxonomy_db(Path(tax_db_loc),
                                            options.tax_db_format)
        except Exception as err:
            raise AnnotationError(f"Failed to import taxonomy database: {err}")
        
        print("import protein to taxonomy map...")
        tax_db_archive_format = determine_archive_format(options.tax_db_name)
        try:
            tax_df = import_acc_tax_map(taxonomy_map,
                                        options.tax_db_acc_idx,
                                        options.tax_db_tax_idx,
                                        options.tax_db_acc_pattern,
                                        options.tax_db_delimiter,
                                        taxonomy_db,
                                        options.tax_annot_format,
                                        options.tax_db_element_format, 
                                        None,
                                        tax_db_archive_format,
                                        options.tax_db_accession_format == "Sequence")
        except Exception as err:
            raise AnnotationError(f"Failed to import taxonomy annotations: {err}")
    else:
        try:
            taxonomy_db = import_taxonomy_db(Path(tax_db_loc),
                                             options.tax_db_format)
        except Exception as err:
            raise AnnotationError(f"Failed to import taxonomy database: {err}")
        tax_df = None
        
    
    if func_annot_map is not None and \
        options.func_db_format is not None and \
        options.func_db_name is not None:
        print("import functional annotation db...")
        # define archive format, get file name by removing archive suffix
        func_annot_archive_format = determine_archive_format(options.func_db_name)
        try:
            func_annot_db = import_func_map(func_annot_map,
                                            options.func_db_format,
                                            accession_pattern=options.func_annot_pattern,
                                            archive_format=func_annot_archive_format)
        except Exception as err:
            raise AnnotationError(f"Failed to import functional annotations: {err}")
    else:
        func_annot_db = None

    # load crap dataset if present
    if Path(GlobalConstants.crap_fasta_loc).exists():
        try:
            crap_pept = fasta_to_peptides(Path(GlobalConstants.crap_fasta_loc))
        except Exception as err:
            raise AnnotationError(f"Failed to import cRAP dataset: {err}")
    else:
        crap_pept = None

    # import de novo datasets into MetaPep format, files that fail to import are ignored
    if len(de_novo_list) != 0:
        print("import de novo files...")
        
        de_novo_crap = crap_pept if options.de_novo_filter_crap is True else None
        try:
            de_novo_dict = match_source_denovo(de_novo_list,
                                               options.de_novo_format,
                                               de_novo_crap)
        except Exception as err:
            raise AnnotationError(f"failed to import de novo data: {err}")
    else:
        de_novo_dict = None
    
    # import data from psm and protein db
    if db_search_list is not None and len(db_search_list) != 0:
        print("import db search files...")
        db_search_crap = crap_pept if options.db_search_filter_crap is True else None
        # assign identical names to peptide blocks if merger desired
        if options.merge_psms is True:
            try:
                psm_df_list = [load_metapep_db_search(db_search_psm,
                                                      sample_name,
                                                      options.db_search_format,
                                                      db_search_crap)
                            for db_search_psm in db_search_list]
            except Exception as err:
                raise AnnotationError(str(err))
            
            # psm_merge_df = pd.concat(psm_df_list, axis=0)]
            try:
                psm_merge_df = MetaPepDbSearch.concat_tables(psm_df_list)
            except Exception as err:
                raise AnnotationError(f"Failed to merge supplied db search files: {err}")

            try:
                metapep_table = build_metapep_table(psm_merge_df,
                                                    tax_df,
                                                    func_annot_db,
                                                    de_novo_dict,
                                                    sample_name,
                                                    taxonomy_db,
                                                    options)
            except Exception as err:
                raise AnnotationError(f"Failed to combine imports into new sample: {err}")
        else:
            metapep_table_list = []
            for i, db_search_psm in enumerate(db_search_list):
                print(f"process {db_search_names[i]}...")
                
                try:
                    psm_df = load_metapep_db_search(db_search_psm,
                                                    db_search_names[i],
                                                    options.db_search_format,
                                                    db_search_crap)
                except Exception as err:
                    raise AnnotationError(str(err))
                
                try:
                    sample_metapep_table = build_metapep_table(psm_df,
                                                               tax_df,
                                                               func_annot_db,
                                                               de_novo_dict,
                                                               db_search_names[i],
                                                               taxonomy_db,
                                                               options)
                except Exception as err:
                    raise AnnotationError(f"Failed to combine imports into new sample: {err}")
                metapep_table_list.append(sample_metapep_table)
                
            try:
                metapep_table = MetaPepTable.concat_tables(metapep_table_list)
            except Exception as err:
                raise AnnotationError(f"Failed to add sample to project table: {err}")
            
    # without db search data, build metapep table with only global taxonomy
    elif de_novo_dict is not None:
        try:
            metapep_table = build_metapep_table(None,
                                                None,
                                                None,
                                                de_novo_dict,
                                                sample_name,
                                                taxonomy_db,
                                                options)
        except Exception as err:
            raise AnnotationError(f"Failed to combine imports into new sample: {err}")
    else:
        raise AnnotationError("No valid db search or de novo data supplied")
    
    return metapep_table


# process the db search psm dataset with taxonomic (and in the future functional) annotation
def taxonomic_annotation(peptides: pd.DataFrame,
                         taxonomy_map: AccessionTaxaMapMethods,
                         taxonomy_db: TaxonomyDatabase,
                         accession_column: str) -> pd.DataFrame:
    """Add taxonomic annotation to peptide dataset
    """
    
    print("annotate taxonomy id to peptides...")
    # perform taxonomic annotation on peptides dataset
    acc_delim = GlobalConstants.peptides_accession_delimiter
    peptides["Taxonomy Id"] = peptides[accession_column].str.split(acc_delim) \
        .apply(taxonomy_map.accession_list_to_lca, taxonomy_db=taxonomy_db)
    
    print("get taxonomy names from id's...") 
    # perform taxonomic annotation on peptides dataset
    peptides["Taxonomy Name"] = peptides["Taxonomy Id"] \
        .apply(taxonomy_db.id_to_name)

    print("add lineage to psm data...")
    # add lineage to dataset
    peptides = add_lineage(peptides, taxonomy_db)
    
    print("finished taxonomy annotation")
    # return taxonomically annotated peptide dataset
    return peptides
    

def accession_list_to_function(acc_list: str | List[str],
                               func_annot_db: FunctionDbMapper,
                               combine: bool=False) -> pd.Series | float:
    """Return functional annotation data for a given list of accessions. In case of
    multiple accessions, either combine (unique) outputs of all accessions or only keep
    annotations present in all accessions.

    Args:
        acc_list (str | List[str]): List of protein accessions.
        func_annot_db (FunctionDbMapper): Functional annotation database
        combine (bool, optional): in case of multiple accessions, specify if 
            annotations should be combined. Defaults to False
            
    Returns:
        pd.Series | float: Series object of eggnog data for that accession list.
            If no data found, return nan.
    """
    # retrieve functional annotation data for all accessions in the list
    func_annot = func_annot_db.data_from_accession(acc_list)
    
    if not isinstance(func_annot, pd.DataFrame):
         return np.nan
        
    # if multiple accessions in list, combine all annotations from the accessions if combine is true
    if func_annot.shape[0] > 1 and combine is True:
        # multiple accessions are joined with ';' inbetween. This gives full information in the output
        # data, but makes it less useful for data analysis, as any combination will register as a different value.
        delim = func_annot_db.LIST_DELIM
        nanfill = func_annot_db.NANFILL
        merged_func_annot = func_annot.apply(
            lambda x: delim.join(x.fillna(nanfill).astype(str)),
            axis=0
        )
    # if combine is false, only include annotations with full consensus
    elif func_annot.shape[0] > 1:
        merged_func_annot = func_annot.apply(
            lambda x: x.iloc[0] if x.nunique(dropna=False) == 1 else np.nan,
            axis=0)
    # if no match found at all, return nan
    elif func_annot.shape[0] == 0:
        return np.nan
    # finally if only one accession present/matched, squeeze single row to series
    else:
        merged_func_annot: pd.Series = func_annot.squeeze()
    
    return merged_func_annot
    
    
def functional_annotation(peptides: pd.DataFrame,
                          function_db: FunctionDbMapper,
                          combine: bool = False) -> pd.DataFrame:
    """Annotate psm data with functional annotation dataset.
    """     
    acc_delim = GlobalConstants.peptides_accession_delimiter
    for i, acc_list in enumerate(peptides["Accession"].str.split(acc_delim)):
        # some peptides in psm dataset may have no accession annotation, these cannot be functionally annotated
        if acc_list is None:
            continue
        
        func_data = accession_list_to_function(acc_list,
                                               function_db,
                                               combine=combine)

        # For first value, check if columns exist. If not, initialize them with correct dtype        
        if i == 0:
            for column in func_data.index:
                if column not in peptides.columns and column == "evalue":
                    peptides.loc[:, column] = pd.Series(dtype="float64")
                elif column not in peptides.columns:
                    peptides.loc[:, column] = pd.Series(dtype="object")
        
        # check if series is returned, then add to dataset
        if isinstance(func_data, pd.Series):
            # to add (new column) data to df in one row, only loc method is compatible. It requires to specify the row index name (loc, not iloc)
            # and the series index as new column names. To keep compatibility if peptides dataset has different index format, row index is taken directly
            # from peptides index array
            idx = peptides.index[i]
            peptides.loc[idx, func_data.index] = func_data

    print("finished function annotation")
    # return taxonomically annotated peptide dataset
    return peptides


def build_metapep_table(metapep_db_search: MetaPepDbSearch | None,
                        taxonomy_mapper: AccessionTaxaMapMethods | None,
                        functional_mapper: FunctionDbMapper | None,
                        denovo_source_map: Dict[str, MetaPepDeNovo] | None,
                        sample_name: str,
                        taxonomy_db: TaxonomyDatabase,
                        options: AnnotationOptions) -> MetaPepTable:
    print("wrangle db search dataset to metapepview format...")
    if metapep_db_search is not None:
        try:
            peptides = process_db_search_data(metapep_db_search, options)
        except Exception as err:
            raise AnnotationError("Failed db search data wrangling.")

        if taxonomy_mapper is not None and options.tax_db_format is not None:
            print("annotate taxonomy...")
            try:
                peptides = taxonomic_annotation(peptides,
                                                taxonomy_mapper,
                                                taxonomy_db,
                                                options.tax_db_accession_format)

            except Exception as err:
                raise AnnotationError("Failed taxonomy annotation of db search peptides.")
        else:
            peptides[GlobalConstants.metapep_table_taxonomy_fields] = np.nan
        
        if functional_mapper is not None:
            print("annotate function...")
            # perform functional annotation
            try:
                peptides = functional_annotation(peptides,
                                                functional_mapper,
                                                combine=options.func_annot_combine)
            except:
                raise AnnotationError("Failed functional annotation of db search peptides.")
        else:
            peptides[GlobalConstants.metapep_table_function_fields] = np.nan
            
        db_search_source_files = metapep_db_search.data["Source File"].unique()
    else:
        peptides = pd.DataFrame(columns=['Peptide', 'Sequence'] +
            GlobalConstants.metapep_table_db_search_fields + 
            GlobalConstants.metapep_table_taxonomy_fields +
            GlobalConstants.metapep_table_function_fields
            )
        
        
        db_search_source_files = None
    
    # add de novo information to peptide table
    try:
        peptides, de_novo_format, de_novo_confidence_format = include_de_novo(
            peptides,
            denovo_source_map, 
            db_search_source_files, # type: ignore
            options)
    except Exception as err:
        raise AnnotationError("Failed to supplement de novo peptides to sample data.")
    
    # perform de novo sequence taxonomy annotation, only supports ncbi format
    if options.global_taxonomy_annotation is True and isinstance(taxonomy_db, NcbiTaxonomy):
        try:
            peptides = global_taxonomic_annotation(peptides, taxonomy_db)
        except Exception as err:
            raise AnnotationError(f"Unipept error: {err}")
        peptides['Global Taxonomy Annotation'] = True
    else:
        peptides['Global Taxonomy Annotation'] = False

    # add general metadata to the dataset 
    peptides['Sample Name'] = sample_name
    peptides['Taxonomy Annotation'] = False if options.tax_db_name is None else True
    peptides['Taxonomy DB Name'] = options.tax_db_name
    peptides['Functional Annotation'] = False if options.func_db_name is None else True
    peptides['Functional Annotation DB Name'] = options.func_db_name
    peptides['DB Search Imported'] = False if metapep_db_search is None else True
    
    if metapep_db_search is not None:
        db_search_data_source = metapep_db_search.data_source
        db_search_confidence_format = metapep_db_search.confidence_format
    else:
        db_search_data_source = None
        db_search_confidence_format = None
    
    # initialize into MetaPepTable object
    try:
        return MetaPepTable(
            peptides,
            options.tax_db_format,
            options.func_db_format,
            db_search_data_source,
            db_search_confidence_format,
            de_novo_format,
            de_novo_confidence_format,
            None)
    except Exception as err:
        raise AnnotationError(f"Generated invalid project table: {err}")


def process_db_search_data(
    metapep_db_search: MetaPepDbSearch,
    options: AnnotationOptions) -> pd.DataFrame:
    """Perform filtering and data wrangling on db search data. The dataset
    is grouped by peptide sequences.

    Args:
        metapep_db_search (MetaPepDbSearch): MetaPepDbSearch object.
        options (AnnotationOptions): Specific annotation options.

    Returns:
        pd.DataFrame: Wrangled and sequence grouped db search dataset.
    """
    # filter psm by score threshold
    if not np.isnan(options.min_db_search_score):
        df = metapep_db_search.data
        metapep_db_search.data = df[df["Confidence"] > options.min_db_search_score]

    # apply regex to protein id's to extract substrings if provided
    re_pattern = options.db_search_acc_pattern
    metapep_db_search = metapep_db_search_accession_regex(metapep_db_search,
                                                          "Accession",
                                                          re_pattern)
        
    # combine group table by peptide sequences, aggregate fields and count matches per peptide
    return metapep_table_to_peptides(metapep_db_search)


def include_de_novo(db_search_peptides: pd.DataFrame | MetaPepDbSearch,
                    de_novo_source_map: Dict[str, MetaPepDeNovo] | None,
                    source_files: Sequence[str],
                    options: AnnotationOptions) -> Tuple[pd.DataFrame,
                                                         DeNovoSource | None,
                                                         DeNovoConfFormat | None]:
    """Extend (annotated) db search dataset with de novo data. The merged
    peptide dataset adheres to the MetaPepTable format. If no de novo data 
    supplied, then the dataset will be extended with empty de novo columns.

    Args:
        db_search_peptides (pd.DataFrame): Db search peptide data.
        de_novo_source_map (Dict[str, MetaPepDeNovo] | None): De novo data
            mapped by corresponding Raw file name.
        source_files (Sequence[str]): Raw file names processed within db search
            dataset.
        options (AnnotationOptions): Annotation options object.

    Returns:
        Tuple[pd.DataFrame, DeNovoSource | None, DeNovoConfFormat | None]:
            Merged dataset, with de novo format and de novo confidence format.
    """
    # if present within metapep db search object, fetch dataset
    if isinstance(db_search_peptides, MetaPepDbSearch):
        # group match based dataset by peptide sequences
        peptides_df = metapep_table_to_peptides(db_search_peptides)
    else:
        peptides_df = db_search_peptides
    
    # add denovo data to sample if present
    if de_novo_source_map is not None and len(de_novo_source_map.keys()) > 0:
        # fetch correct de novo data from source map, only taking source files in de novo data
        if options.merge_psms is False:
            de_novo_source_files = np.intersect1d(source_files, 
                                                  list(de_novo_source_map.keys()))
            denovo_data = [de_novo_source_map[source_file] for source_file in de_novo_source_files]
        # if all data under single sample name, process all files
        else:
            denovo_data = list(de_novo_source_map.values())
        
        if len(denovo_data) == 0:
            peptides_df['De Novo Imported'] = False
            return peptides_df, None, None
        
        print('add de novo data...')
        metapep_de_novo = MetaPepDeNovo.concat_tables(denovo_data)

        # wrangle de novo data and merge into peptides dataset
        de_novo_peptides = process_de_novo_data(metapep_de_novo, options)
        
        peptides_df = merge_de_novo_db_search(peptides_df, 
                                              de_novo_peptides)
        
        # get denovo metadata
        de_novo_format = denovo_data[0].data_source
        de_novo_confidence_format = denovo_data[0].confidence_format
        peptides_df['De Novo Imported'] = True
    else:
        peptides_df[GlobalConstants.metapep_table_de_novo_fields] = np.nan
        peptides_df['De Novo Imported'] = False
        de_novo_format = None
        de_novo_confidence_format = None
        
    return peptides_df, de_novo_format, de_novo_confidence_format
    

def process_de_novo_data(
    metapep_de_novo: MetaPepDeNovo,
    options: AnnotationOptions) -> pd.DataFrame:
    """Perform filtering and data wrangling on de novo data. The dataset is
    grouped by peptide sequences.

    Args:
        metapep_de_novo (Sequence[MetaPepDeNovo]): MetaPepDeNovo object
        options (AnnotationOptions): Options settings object.

    Returns:
        pd.DataFrame: Wrangled and sequence grouped de novo dataset.
    """
    print("match de novo peptides to dataset...")
    # filter denovo to only peptides above score threshold
    if not np.isnan(options.min_de_novo_score):
        df = metapep_de_novo.data
        metapep_de_novo.data = df[df['Confidence'] > options.min_de_novo_score]
    if not np.isnan(options.min_pept_length):
        df = metapep_de_novo.data
        metapep_de_novo.data = df[df['Length'] > options.min_pept_length]
    
    # group de novo identification table to peptides
    return metapep_table_to_peptides(metapep_de_novo)


def merge_de_novo_db_search(db_search_peptides: pd.DataFrame,
                            de_novo_peptides: pd.DataFrame,
                            match_source_file: bool = False) -> pd.DataFrame:
    """Add de novo dataset information to db search peptides dataset.
    The concatenated dataset adheres to the format laid out for the MetaPepTable
    object.

    Args:
        de_novo_peptides (pd.DataFrame): de novo peptides dataset.
        db_search_peptides (pd.DataFrame): db search peptides dataset

    Returns:
        pd.DataFrame: Merged dataset
    """
    # filter out unrelevant de novo columns, rename to differentiate denovo from db search
    de_novo_fields = GlobalConstants.metapep_de_novo_fields
    de_novo_peptides = de_novo_peptides[['Peptide', 'Sequence'] + de_novo_fields]
    de_novo_peptides = de_novo_peptides.rename(
        columns={field: 'De Novo ' + field for field in de_novo_fields + ['Peptide']})
    
    # merge denovo data into peptides dataset
    if match_source_file is True:
        de_novo_peptides["Source File"] = de_novo_peptides["De Novo Source File"]
        merge_df = pd.merge(db_search_peptides, de_novo_peptides, how="outer", on=["Sequence", "Source File"])
    else:
        merge_df = pd.merge(db_search_peptides, de_novo_peptides, how="outer", on=["Sequence"])
        # if no db search, peptide field is empty, in that case, fill the values
        merge_df.loc[:, "Peptide"] = merge_df["Peptide"].fillna(merge_df["De Novo Peptide"])
    
    merge_df = merge_df.drop(columns='De Novo Peptide')
    return merge_df
        
    