from typing import Dict, Type

from ..types import *
    

taxonomy_db_formats: Dict[str, Type[AccessionTaxaMap]] = {
    'NCBI': AccessionTaxaMapNcbi,
    'GTDB': AccessionTaxaMapGtdb
}

taxonomy_db_importers = {
    'NCBI': NcbiTaxonomy.from_dmp_folder,
    'GTDB': GtdbTaxonomy.from_tsv_folder
}

functional_db_importers: Dict[str, Type[FunctionDbMapper]] = {
    'EggNOG': EggnogMapper,
    'gKOALA': KeggMapper
}


db_search_importers: Dict[str, Type[DbSearchMethods]] = {
    'Peaks 11': PeaksDbSearchPsm11,
    'Peaks 10': PeaksDbSearchPsm10,
    'MaxQuant': MaxQuantDbSearch,
    'Sage': SageDbSearch
    #'ProteomeDiscoverer': ...
}

de_novo_importers: Dict[str, Type[DeNovoMethods]] = {
    'Peaks 11': PeaksDeNovo11,
    'Peaks 10': PeaksDeNovo10,
    'Novor': NovorDeNovo
}