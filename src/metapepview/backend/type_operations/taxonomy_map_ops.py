from metapepview.backend.types import AccessionTaxaMapGtdb,\
    AccessionTaxaMapNcbi,\
    AccessionTaxaMap,\
    TaxonomyDbFormat,\
    TaxonomyMapFormat,\
    TaxonomyElementFormat,\
    TaxonomyDatabase,\
    GtdbGenomeToNcbi,\
    GhostkoalaMapper
from metapepview.backend.types.taxonomy_map.base_class import AccessionTaxaMapMethods
from metapepview.backend.utils import *
from metapepview.backend.types.object_mappings import taxonomy_db_formats



def import_acc_tax_map(upload_contents: str | IO[str],
                       acc_col: int,
                       tax_col: int,
                       acc_regex: str | None,
                       delimiter: str | None,
                       taxonomy_db: TaxonomyDatabase,
                       taxonomy_db_format: TaxonomyMapFormat,
                       taxonomy_element_format: TaxonomyElementFormat,
                       gtdb_to_ncbi_obj: GtdbGenomeToNcbi | None = None,
                       archive_format: str | None = None,
                       wrangle_peptide_accessions = False) -> AccessionTaxaMapMethods:
    # if file buffer has not been extracted from raw string data, extract data into buffer
    if isinstance(upload_contents, str):
        str_file_obj = memory_to_stringio(upload_contents, archive_format)
    else:
        str_file_obj = upload_contents

    # construct correct class object based on taxonomy db format
    tax_format = taxonomy_db_formats.get(taxonomy_db_format, None)

    # if tax names in mapping data, convert to id
    name_to_id = True if taxonomy_element_format == "taxonomy name" else False
    
    if taxonomy_db_format == "GTDB" and gtdb_to_ncbi_obj is not None:
        return AccessionTaxaMapNcbi.from_gtdb_genome_ids(
            gtdb_to_ncbi_obj,
            str_file_obj,
            acc_col,
            tax_col,
            acc_regex,
            delimiter,
            drop_duplicates=False,
            tax_name_to_id=name_to_id,
            taxonomy_obj=taxonomy_db,
            wrangle_peptide_accessions=wrangle_peptide_accessions
            )
    elif taxonomy_db_format == "GhostKOALA":
        return tax_format.read_file_buffer(
            str_file_obj,
            acc_regex
        )
    elif tax_format is not None:
        return tax_format.from_string_buffer(
            str_file_obj,
            acc_col,
            tax_col,
            acc_regex,
            delimiter,
            drop_duplicates=False,
            tax_name_to_id=name_to_id,
            taxonomy_obj=taxonomy_db,
            wrangle_peptide_accessions=wrangle_peptide_accessions
            )
    else:
        raise ValueError("Invalid taxonomy format")


def validate_acc_tax_map(upload_contents: str | IO[str],
                         acc_col: int,
                         tax_col: int,
                         delimiter: str,
                         taxonomy_db_format: TaxonomyMapFormat,
                         taxonomy_element_format:TaxonomyElementFormat,
                         archive_format: str | None = None) -> Tuple[bool, str | None]:
    # if file buffer has not been extracted from raw string data, extract data into buffer
    if isinstance(upload_contents, str):
        str_file_obj = memory_to_stringio(upload_contents, archive_format)
    else:
        str_file_obj = upload_contents
        # construct correct class object based on taxonomy db format
    tax_format = taxonomy_db_formats.get(taxonomy_db_format, None)
    
    if issubclass(tax_format, AccessionTaxaMap):
        return tax_format.validate_input(
            str_file_obj,
            acc_col,
            tax_col,
            delimiter,
            element_format=taxonomy_element_format)
    elif issubclass(tax_format, GhostkoalaMapper):
        try:
            tax_format.read_file_buffer(str_file_obj)
            return True, None
        except Exception as e:
            return False, repr(e)
    else:
        raise ValueError("Invalid taxonomy format")