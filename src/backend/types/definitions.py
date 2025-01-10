"""This module stores type definitions that are used throughout the program.
"""

from typing import Literal


# Literal types

# db search format options
DbSearchSource = Literal['Peaks', 'Peaks10', 'MaxQuant', 'ProteomeDiscoverer', 'Sage']
DeNovoSource = Literal['Peaks', 'Peaks10', 'Novor']

# de novo format options
DbSearchConfFormat = Literal['-10lgp', 'Hyperscore']
DeNovoConfFormat = Literal['ALC', 'Score']

# annotation db format options
TaxonomyFormat = Literal['NCBI', 'GTDB']
FuncAnnotFormat= Literal['gKOALA', 'EggNOG']



