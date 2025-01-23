from typing import Dict, Any

from backend.types import MetaPepDbSearch, MetaPepDeNovo

def assert_db_search_compatibility(
    spectral_ref_data: Dict[str, Dict[str, Any]],
    db_search: MetaPepDbSearch) -> bool:
    """Check that db search file is compatible with the reference
    dataset. Different file sources or processing methods impact
    the results of the metrics processing.

    Args:
        spectral_ref_data (Dict[str, Dict[str, Any]]): Reference dataset
        db_search (MetaPepDbSearch): Db search object

    Returns:
        bool: True if db search file is compatible with ref data.
    """
    ref_metadata = spectral_ref_data['metadata']
    
    # assert format
    return ref_metadata['db search format'] == db_search.data_source
    

def assert_de_novo_compatibility(
    spectral_ref_data: Dict[str, Dict[str, Any]],
    de_novo: MetaPepDeNovo) -> bool:
    """Check that de novo file is compatible with the reference
    dataset. Different file sources or processing methods impact
    the results of the metrics processing.

    Args:
        spectral_ref_data (Dict[str, Dict[str, Any]]): Reference dataset.
        de_novo (MetaPepDeNovo): De novo object.

    Returns:
        bool: True if de novo file is compatible with ref data.
    """
    ref_metadata = spectral_ref_data['metadata']
    
    # assert format
    return ref_metadata['de novo format'] == de_novo.data_source
    