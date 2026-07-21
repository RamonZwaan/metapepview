"""This module stores type definitions that are used throughout the program.
"""
from __future__ import annotations

from typing import Literal, TypeAlias


# Literal types

# db search format options
DbSearchSource: TypeAlias = Literal['Peaks 13', 
                                    'Peaks 11', 
                                    'Peaks 10', 
                                    'MaxQuant', 
                                    'FragPipe',
                                    #'ProteomeDiscoverer', 
                                    'Sage']
DeNovoSource: TypeAlias = Literal['Peaks 13', 
                                  'Peaks 11', 
                                  'Peaks 10', 
                                  'Novor', 
                                  'Casanovo']

# de novo format options
DbSearchConfFormat: TypeAlias = Literal['-10lgp', 'Hyperscore']
DeNovoConfFormat: TypeAlias = Literal['ALC', 'Score']

# annotation db format options
TaxonomyDbFormat: TypeAlias = Literal['NCBI', 'GTDB']
TaxonomyMapFormat: TypeAlias = Literal['NCBI', 'GTDB', 'GhostKOALA']
TaxonomyElementFormat: TypeAlias = Literal["taxonomy id", "taxonomy name"]
FuncAnnotFormat: TypeAlias = Literal['GhostKOALA', 'EggNOG']



