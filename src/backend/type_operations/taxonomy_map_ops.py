from typing import TypeVar

from ..types import AccessionTaxaMapGtdb,\
    AccessionTaxaMapNcbi,\
    AccessionTaxaMap,\
    TaxonomyFormat,\
    TaxonomyDatabase
from backend.utils import *
from .object_mappings import taxonomy_db_formats



def import_acc_tax_map(upload_contents: str | IO[str],
                       acc_col: int,
                       tax_col: int,
                       acc_regex: str | None,
                       delimiter: str | None,
                       taxonomy_db: TaxonomyDatabase,
                       taxonomy_db_format: TaxonomyFormat,
                       archive_format: str | None = None) -> AccessionTaxaMap:
    # if file buffer has not been extracted from raw string data, extract data into buffer
    if isinstance(upload_contents, str):
        str_file_obj = memory_to_stringio(upload_contents, archive_format)
    else:
        str_file_obj = upload_contents

    # construct correct class object based on taxonomy db format
    tax_format = taxonomy_db_formats.get(taxonomy_db_format, None)
    
    if tax_format is not None:
        return tax_format.from_string_buffer(
            str_file_obj,
            acc_col,
            tax_col,
            acc_regex,
            delimiter,
            drop_duplicates=False,
            taxonomy_obj=taxonomy_db
            )
    else:
        raise ValueError("Invalid taxonomy format")


def validate_acc_tax_map(upload_contents: str | IO[str],
                         acc_col: int,
                         tax_col: int,
                         delimiter: str,
                         taxonomy_db_format: TaxonomyFormat,
                         archive_format: str | None = None) -> Tuple[bool, str | None]:
    # if file buffer has not been extracted from raw string data, extract data into buffer
    if isinstance(upload_contents, str):
        str_file_obj = memory_to_stringio(upload_contents, archive_format)
    else:
        str_file_obj = upload_contents
        # construct correct class object based on taxonomy db format
    tax_format = taxonomy_db_formats.get(taxonomy_db_format, None)
    
    if tax_format is not None:
        return tax_format.validate_input(
            str_file_obj,
            acc_col,
            tax_col,
            delimiter)
    else:
        raise ValueError("Invalid taxonomy format")