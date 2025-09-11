#!/usr/bin/env python

import sys
from pathlib import Path
import argparse
from textwrap import dedent
import pandas as pd

from metapepview.backend.spectral_ref_builder import build_ref_data



def main():
    # set downcasting behavior to manage FutureWarning in `replace` function
    pd.set_option('future.no_silent_downcasting', True)

    # set CLI arguments
    parser = argparse.ArgumentParser(
        prog='buildref',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=dedent(
            """Create benchmark dataset of metaproteomics experiments from set
            of experimental data.
            
            It parses a supplied root directory for all relevant experimental 
            data files (mzML, peptide identification data from db search and de 
            novo). All files are matched by their common source file name and
            performance metrics are extracted from the data to construct a 
            performance statistics dataset.
            """
            )
        )
    
    # define arguments for the tool
    parser.add_argument('-d', 
                        '--db-search-format', 
                        required=True,
                        choices=['peaks10', 'peaks11', 'maxquant', 'sage']
                        )
    parser.add_argument('-n', 
                        '--de-novo-format', 
                        required=True,
                        choices=['peaks10', 'peaks11', 'novor', 'casanovo']
                        )
    parser.add_argument('-o', 
                        '--output', 
                        default="./out.json"
                        )
    parser.add_argument('-D', 
                        '--db-search-thresholds', 
                        default='30 50 80',
                        help="Threshold values used when splitting db search matches into confidence brackets."
                        )
    parser.add_argument('-N', 
                        '--de-novo-thresholds', 
                        default='50 80 90',
                        help="Threshold values used when splitting de novo identifications into confidence brackets."
                        )
    parser.add_argument('-I', 
                        '--intensity-percentiles', 
                        default='50 90 99',
                        help="MS signal intensity percentiles to be calculated."
                        )
    parser.add_argument('-T', 
                        '--transmission-loss-percentiles', 
                        default='50 90 99',
                        help='Ion transmission loss percentiles to be calculated.'
                        )
    parser.add_argument('-r', 
                        '--db-search-file-pattern',
                        help="Custom regex pattern used to fetch db search match files by file name."
                        )
    parser.add_argument('-x', 
                        '--de-novo-file-pattern',
                        help="Custom regex pattern used to fetch de novo identification files by file name."
                        )
    parser.add_argument('directory')
    
    # parse arguments
    args = parser.parse_args()
    
    output_path = Path(args.output)
    filename = output_path.name
    if not filename.endswith('.json'): filename += '.json'
    
    # convert threshold specification to list of floats
    try:
        db_search_thresholds = [float(x) for x in args.db_search_thresholds.split()]
    except ValueError:
        sys.exit(f"error: invalid value for --db-search-thresholds: {args.db_search_thresholds}. It must contain a string of spaced numbers, e.g. '50 80 90'")
    try:
        de_novo_thresholds = [float(x) for x in args.de_novo_thresholds.split()]
    except ValueError:
        sys.exit(f"error: invalid value for --de-novo-thresholds: {args.db_search_thresholds}. It must contain a string of spaced numbers, e.g. '50 80 90'")
    try:
        intensity_percentiles = [int(x) for x in args.intensity_percentiles.split()]
    except ValueError:
        sys.exit(f"error: invalid value for --intensity-percentiles: {args.db_search_thresholds}. It must contain a string of spaced numbers, e.g. '50 80 90'")
    try:
        transmission_loss_percentiles = [int(x) for x in args.transmission_loss_percentiles.split()]
    except ValueError:
        sys.exit(f"error: invalid value for --transmission-loss-percentiles: {args.db_search_thresholds}. It must contain a string of spaced numbers, e.g. '50 80 90'")

    # assign correct db search format from input
    match args.db_search_format:
        case 'peaks10':
            db_search_format = 'Peaks 10'
        case 'peaks11':
            db_search_format = 'Peaks 11'
        case 'maxquant':
            db_search_format = 'MaxQuant'
        case 'sage':
            db_search_format = 'Sage'
        case _:
            sys.exit("Invalid de novo format provided.")

    # assign correct db search format from input
    match args.de_novo_format:
        case 'peaks10':
            de_novo_format = 'Peaks 10'
        case 'peaks11':
            de_novo_format  = 'Peaks 11'
        case 'novor':
            de_novo_format = 'Novor'
        case 'casanovo':
            de_novo_format = 'Casanovo'
        case _:
            sys.exit("Invalid de novo format provided.")

    # build the reference dataset
    build_ref_data(
        Path(args.directory),
        output_path.parent,
        filename,
        db_search_format,
        de_novo_format,
        db_search_thresholds,
        de_novo_thresholds,
        intensity_percentiles,
        transmission_loss_percentiles,
        args.db_search_file_pattern,
        args.de_novo_file_pattern,
    )
    
    
if __name__ == "__main__":
    sys.exit(main())