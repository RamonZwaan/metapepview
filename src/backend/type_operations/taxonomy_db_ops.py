from ..types import TaxonomyFormat, TaxonomyDatabase
from backend.utils import *
from .object_mappings import taxonomy_db_importers


def import_taxonomy_db(folder_loc: Path | str,
                       taxonomy_db_format: TaxonomyFormat) -> TaxonomyDatabase:
    if isinstance(folder_loc, str):
        folder_loc = Path(folder_loc)

    # construct correct class object based on taxonomy db format
    obj_constructor = taxonomy_db_importers.get(taxonomy_db_format, None)
    
    if obj_constructor is not None:
        return obj_constructor(folder_loc)
    else:
        raise ValueError("Invalid taxonomy format")
    

def add_lineage(dataset: pd.DataFrame,
                taxonomy_db: TaxonomyDatabase,
                fill_gaps: bool = True,
                tax_col: str = "Taxonomy Id"):
    # retrieve lineages for all peptide sequences and convert tax ids to tax names
    lineage = dataset[tax_col].transform(taxonomy_db.id_to_standard_lineage)
    
    # fill undefined rank values with names defined higher in the lineage
    if fill_gaps is True:
        lineage = taxonomy_db.fill_lineage_gaps(lineage)
    
    lineage_names = lineage.transform(taxonomy_db.lineage_id_to_name)

    # Add both lineage of ids and names to grouped dataset
    lin = GlobalConstants.standard_lineage_ranks
    lin_id_df = pd.DataFrame(lineage.tolist(),
                                index=lineage.index,
                                columns=[i + " Id" for i in lin])
    lin_name_df = pd.DataFrame(lineage_names.tolist(),
                                index=lineage_names.index,
                                columns=[i + " Name" for i in lin])
    
    # synchronize lineage indices to indices of dataset
    lin_id_df.index = dataset.index
    lin_name_df.index = dataset.index
    return pd.concat([dataset, lin_id_df, lin_name_df], axis=1)
    