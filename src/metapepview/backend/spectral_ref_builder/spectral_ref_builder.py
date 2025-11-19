# ms functions
from typing import List, Dict, Any
from pathlib import Path
import pandas as pd

from metapepview.backend.type_operations import * # load_metapep_db_search, load_metapep_de_novo, db_search_importers, de_novo_importers
from metapepview.backend.utils import *

from metapepview.backend.spectral_ref_builder.utils import *
from metapepview.backend.spectral_ref_builder.definitions import *  


def fetch_file_locations(options: RefBuilderOptions) -> Dict[str, Dict[str, Path | None]]:
    """Parse supplied parent dictionary and fetch all files that adhere to the
    data formats of interest. From the collected files, build a mapping dataset
    that maps raw file names to spectral data and metaproteomics output data.
    
    Note:
        Each raw file must be accompanied by single files for each format. that
        is, only one db search file, one de novo file, and one spectral file.
        db search and de novo identification files may contain multiple source
        files (for example files of multiple fractions). In that case,
        the data file will be mapped to all dictionary keys (raw file names) 
        that are sources from the file. In the case that multiple files match
        to the same source file name, the data will be overwritten and only the
        last analyzed file will be mapped to the raw source file name.

    Args:
        options (RefBuilderOptions): Reference dataset builder options.

    Returns:
        Dict[str, Dict[str, Path | None]]: Mapping dataset matching raw file name
            to spectral and metaproteomics output files.
    """
    # fetch all raw files, db search files and de novo files
    root_dir = options.root_dir
    raw_files = list(root_dir.rglob("*.raw"))
    mzxml_files = root_dir.rglob("*.mzXML")
    mzml_files = root_dir.rglob("*.mzML")
    db_search_files = root_dir.rglob(options.db_search_file_pattern)
    de_novo_files = root_dir.rglob(options.de_novo_file_pattern)

    # create dict that combines raw files with db search and denovo files
    output_dict = {key.stem: {"raw": key,
                              "db search": None, 
                              "de novo": None,
                              "mzxml": None,
                              "mzml": None} for key in raw_files}
    
    # parse file lists and update dict
    for mzxml_path in mzxml_files:
        source_names = [mzxml_path.stem]
        output_dict = add_to_source_dict(output_dict, source_names, mzxml_path, "mzxml")
    for mzml_path in mzml_files:
        source_names = [mzml_path.stem]
        output_dict = add_to_source_dict(output_dict, source_names, mzml_path, "mzml")
    for db_search_path in db_search_files:
        db_search_obj = read_ident_file(db_search_path, "db search", options)
        if check_valid_ident_file(db_search_obj, options) is True:
            source_names = ident_file_source(db_search_obj)
            output_dict = add_to_source_dict(output_dict, source_names, db_search_path, "db search")
    for de_novo_path in de_novo_files:
        de_novo_obj = read_ident_file(de_novo_path, "de novo", options)
        if check_valid_ident_file(de_novo_obj, options) is True:
            source_names = ident_file_source(de_novo_obj)
            output_dict = add_to_source_dict(output_dict, source_names, de_novo_path, "de novo")
    
    return output_dict


def build_spectral_reference_data(
    file_loc_dict: Dict[str, Dict[str, Path | None]],
    options: RefBuilderOptions) -> Dict[str, Dict[str, Any]]:
    """Fetch relevant data from large dataset of raw spectral files and
    (meta)proteomics (db search + de novo) files and combine them in a json
    formatted reference dataset containing key metrics across a large set of
    proteomics experiments. This method enables processing of many GB's of data
    from 10's - 100's of experiments into a MB's sized metrics dataset suitable
    for visualization.

    Args:
        file_loc_dict (Dict[str, Dict[Path | None]]): Dictionary of file locations
            (raw data, db search and de novo) grouped by analyzed sample.
        db_search_threhsholds (List[int | float]): List of threshold values that
            partition db search rows into groups by confidence.
        de_novo_thresholds (List[int | float]): List of threshold values that
            partition de novo rows into groups by confidence.
        intensity_percentiles (List[int | float]): List of percentile thresholds
            that partition scan rows by intensity.

    Returns:
        Dict[str, Dict[str, Any]]: Dataset of key metrics grouped by experiment
            extracted from the input datasets. The output dataset is written in
            the following convention:

            {
                "metadata": {
                    "sample size": val,
                    "db search format": val,
                    "db search confidence format": val,
                    "de novo format": val,
                    "de novo confidence format": val,
                    "db search thresholds": [...],
                    "de novo thresholds": [...],
                    "intensity percentiles": [...]
                },
                "global": {
                    "db search confidence dist": {
                        "mean": [...],
                        "std": [...],
                    },
                    "db search confidence dist norm": {
                        "mean": [...],
                        "std": [...],
                        "x vals": [...],
                    },
                    "db search confidence dist norm ms2": {
                        "mean": [...],
                        "std": [...],
                        "x vals": [...],
                    },
                    "de novo confide dist": {
                        "mean": [...],
                        "std": [...],
                    },
                    "de novo confidence dist norm": {
                        "mean": [...],
                        "std": [...],
                        "x vals": [...],
                    },
                    "de novo confidence dist norm ms2": {
                        "mean": [...],
                        "std": [...],
                        "x vals": [...],
                    }
                },
                "samples": {
                    'x.raw': {
                        'timestamp': val,
                        'ms1 count': val
                        'ms2 count': val,
                        'total rt: val
                        'mean pept len': val,
                        'de novo mean pept len': val,
                        'mean mass': val,
                        'de novo mean mass': val,
                        'db search matches: val,
                        'de novo matches': val,
                        'de novo mean mass': val,
                        'mean transmission loss': val,
                        'ion injection time': val,
                        'de novo confidence dist' : {
                            'thresholds': [a, b, c],
                            'counts': [x1, x2, x3]
                        },
                        'db search confidence dist' : {
                            'thresholds': [a, b, c],
                            'counts': [y1, y2, y3]
                        },
                        'de novo only confidence dist': {
                            'thresholds': [a, b, c],
                            'counts': [z1, z2, z3]
                        },
                        'charge dist': {
                            'charge': [1, 2, 3, 4],
                            'counts': [z1, z2, z3, z4]
                        },
                        'miscleave dist': {
                            'miscleavage': [0, 1, 2, 3],
                            'counts': [m0, m1, m2, m3]
                        },
                        'ms1 intensity': {
                            'percentiles': [d, e, f],
                            'values': [t1, t2, t3]
                        },
                        'ms2 intensity': {
                            'percentiles': [d, e, f],
                            'values': [t1, t2, t3]
                        },
                        'transmission loss': {
                            'percentiles: [a, b, c],
                            'values: [s1, s2, s3]
                        },
                        'transmission loss ion injection time scaled': {
                            'percentiles': [a, b, c],
                            'values': [sc1, sc2, sc3]
                        },
                    },
                    ...
                }
            }
    """
    # initialize output dict
    statistics_dict = {
        'global': {
            'db search confidence dist': {
                'mean': [],
                'std': []    
            },
            'db search confidence dist norm': {
                'mean': [],
                'std': [],
                'x vals': []    
            },
            'db search confidence dist norm ms2': {
                'mean': [],
                'std': [],
                'x vals': []    
            },
            'de novo confidence dist': {
                'mean': [],
                'std': []
            },
            'de novo confidence dist norm': {
                'mean': [],
                'std': [],
                'x vals': [] 
            },
            'de novo confidence dist norm ms2': {
                'mean': [],
                'std': [],
                'x vals': [] 
            },
        },
        'samples': { 
        },
        'metadata': {
            "sample size": len(file_loc_dict.keys()),
            "db search format": options.db_search_format,
            "db search confidence format": db_search_importers[options.db_search_format].CONFIDENCE_FORMAT,
            "de novo format": options.de_novo_format,
            "de novo confidence format": de_novo_importers[options.de_novo_format].CONFIDENCE_FORMAT,
            "db search thresholds": options.db_search_thresholds,
            "de novo thresholds": options.de_novo_thresholds,
            "intensity percentiles": options.intensity_percentiles,
        }
    }
    
    
    # first, compute lgp and alc distributions, add them to the output dict
    db_search_conf_mean, db_search_conf_std = score_rank_dist(file_loc_dict,
                                                              'db search',
                                                              options.db_search_format,
                                                              'Confidence')
    de_novo_conf_mean, de_novo_conf_std = score_rank_dist(file_loc_dict,
                                                          'de novo', 
                                                          options.de_novo_format,
                                                          'Confidence')

    db_search_norm_mean, db_search_norm_std, db_search_norm_x_vals = score_rank_dist_norm(
        file_loc_dict,
        'db search',
        options.db_search_format,
        'Confidence')
    de_novo_norm_mean, de_novo_norm_std, de_novo_norm_x_vals = score_rank_dist_norm(
        file_loc_dict,
        'de novo',
        options.de_novo_format,
        'Confidence')
    
    statistics_dict['global']['db search confidence dist'] = {
        'mean': db_search_conf_mean.tolist(),
        'std': db_search_conf_std.tolist()
    }
    statistics_dict['global']['de novo confidence dist'] = {
        'mean': de_novo_conf_mean.tolist(),
        'std': de_novo_conf_std.tolist()
    }
    statistics_dict['global']['db search confidence dist norm'] = {
        'mean': db_search_norm_mean.tolist(),
        'std': db_search_norm_std.tolist(),
        'x vals': db_search_norm_x_vals.tolist()
    }
    statistics_dict['global']['de novo confidence dist norm'] = {
        'mean': de_novo_norm_mean.tolist(),
        'std': de_novo_norm_std.tolist(),
        'x vals': de_novo_norm_x_vals.tolist()
    }

    # counter
    c = 0
    # parse each file separately to compute sample specific statistics
    for name, data in file_loc_dict.items():
        # initialize empty sample field
        sample_output = {
            'timestamp': np.nan,
            'ms1 count': np.nan,
            'ms2 count': np.nan,
            'total rt': np.nan,
            'mean pept len': np.nan,
            'de novo mean pept len': np.nan,
            'mean mass': np.nan,
            'de novo mean mass': np.nan,
            'db search matches': np.nan,
            'de novo matches': np.nan,
            'de novo confidence dist' : {
                'thresholds': [],
                'counts': []
            },
            'db search confidence dist' : {
                'thresholds': [],
                'counts': []
            },
            'de novo only confidence dist': {
                'thresholds': [],
                'counts': []
            },
            'charge dist': {
                'charge': [],
                'counts': []
            },
            'miscleave dist': {
                'miscleavage': [],
                'counts': []
            },
            'ms1 intensity': {
                'percentiles': [],
                'values': []
            },
            'ms2 intensity': {
                'percentiles': [],
                'values': []
            },
            'transmission loss': {
                'percentiles': [],
                'values': []
            },
            'transmission loss ion injection time scaled': {
                'percentiles': [],
                'values': []
            },
        }
        # Open psm and de novo datasets, if they exist
        if data["db search"] is None:
            db_search, db_search_data = None, None
        else:
            db_search = load_metapep_db_search(data["db search"].open('r'), 
                                               name, 
                                               options.db_search_format)
            # if file contains data from other source files, omit them
            if len(db_search.source_files) > 1:
                db_search = db_search.filter_spectral_name(name)
            db_search_data = db_search.data
            
        if data["de novo"] is None:
            de_novo, de_novo_data = None, None
        else:
            de_novo = load_metapep_de_novo(data['de novo'].open('r'),
                                           name,
                                           options.de_novo_format)
            if len(de_novo.source_files) > 1:
                de_novo = de_novo.filter_spectral_name(name)
            de_novo_data = de_novo.data
        
        # wrangle datasets
        # fetch data from db search
        if db_search_data is not None:
            sample_output['mean pept len'] = db_search_data['Sequence'].apply(len).mean()
            sample_output['mean mass'] = db_search_data['Mass'].mean()
            sample_output['db search matches'] = db_search_data.shape[0]
            db_search_thres_counts, db_search_thres_names = count_threshold_values(
                db_search_data,
                "Confidence",
                options.db_search_thresholds
            )
            sample_output['db search confidence dist'] = {
                'thresholds': db_search_thres_names,
                'counts': db_search_thres_counts
            }

            # add miscleavage distribution to sample
            miscleave_groups, miscleave_counts = calculate_miscleavages(db_search_data)
            sample_output['miscleave dist'] = {
                'miscleavage': miscleave_groups,
                'counts': list(miscleave_counts)
            }

        # fetch data from de novo
        if de_novo_data is not None:
            sample_output['de novo mean pept len'] = de_novo_data['Sequence']\
                .apply(len)\
                .mean()
            sample_output['de novo mean mass'] = de_novo_data['Mass'].mean()
            sample_output['de novo matches'] = de_novo_data.shape[0]
            de_novo_thres_counts, de_novo_thres_names = count_threshold_values(
                de_novo_data,
                "Confidence",
                options.de_novo_thresholds)
            sample_output['de novo confidence dist'] = {'thresholds': de_novo_thres_names,
                                         'counts': de_novo_thres_counts}
        
        # fetch de novo only, which requires both de novo and db search
        if de_novo is not None and db_search is not None:
            de_novo_only = de_novo.filter_de_novo_only(db_search)
            de_novo_only_counts, de_novo_only_names = count_threshold_values(
                de_novo_only.data,
                "Confidence",
                options.de_novo_thresholds
                )
            sample_output['de novo only confidence dist'] = {
                'thresholds': de_novo_only_names,
                'counts': de_novo_only_counts
            }
        

        # fetch spectral information from mzml file
        spectral = data["mzml"]
        # fetch spectral data from mzxml file
        if spectral is not None:
            mzml_fields = ['scan number', 'MS level',
                           'peaks count', 'retention time', 'total ion current', 
                           'precursor charge', 'precursor intensity', 
                           'ion injection time', 'precursor scan number']
            
            print(c, spectral)
            c += 1
            
            try:
                mzml_df, mzml_metadata = mzml_to_df(open(spectral, 'rb'),
                                                    fields=mzml_fields)
            except:
                print("failed to read mzml file.")
                mzml_df = None
                
            if mzml_df is not None and mzml_df.shape[0] > 0: 
                sample_output['timestamp'] = mzml_metadata['timestamp']
                # make complete dataframe numeric
                mzml_df = mzml_df.apply(pd.to_numeric, axis=0)
                
                # obtain statistics from mzxml file
                # mslevel_idx = mzxml_fields.index()
                ms2_df = mzml_df[mzml_df['MS level'] == 2]
                ms1_df = mzml_df[mzml_df['MS level'] == 1]
                sample_output['ms2 count'] = ms2_df.shape[0]
                sample_output['ms1 count'] = ms1_df.shape[0]
                
                sample_output['total rt'] = mzml_df.at[mzml_df.index[-1],
                                                       'retention time']
                
                
                # compute ms2 intensities at median and top percentiles
                ms2_int_vals, ms2_int_names = scan_intensity_percentiles(ms2_df,
                                                                         'total ion current',
                                                                         options.intensity_percentiles)
                ms1_int_vals, ms1_int_names = scan_intensity_percentiles(ms1_df,
                                                                         'total ion current',
                                                                         options.intensity_percentiles)
                sample_output['ms2 intensity'] = {'percentiles': ms2_int_names,
                                                  'values': ms2_int_vals}
                sample_output['ms1 intensity'] = {'percentiles': ms1_int_names,
                                                  'values': ms1_int_vals}

                # compute transmission losses, only if ms2 spectra present
                if ms2_df.shape[0] == 0:
                    sample_output['transmission loss'] = {
                        'percentiles': ms2_int_names,
                        'values': [np.nan, np.nan, np.nan]}
                    sample_output['transmission loss ion injection time scaled'] = {
                        'percentiles': ms2_int_names,
                        'values': [np.nan, np.nan, np.nan]}
                else:
                    # fetch precursor ion injection time
                    prec_inj_time = fetch_precursor_ion_injection_time(
                        ms2_df, ms1_df
                    )

                    # compute transmission losses (ion injection time scaled and unscaled)
                    trm_loss_ion_inj = transmission_loss(ms2_df['precursor intensity'],
                                                         ms2_df['total ion current'],
                                                         prec_inj_time,
                                                         ms2_df['ion injection time'])
                    trm_loss = transmission_loss(ms2_df['precursor intensity'],
                                                ms2_df['total ion current'])

                    # compute percentiles for both
                    trm_loss_ion_inj_vals, trm_loss_ion_inj_names = array_to_percentiles(
                        trm_loss_ion_inj.to_numpy(),
                        options.transmission_loss_percentiles,
                        True)
                    trm_loss_vals, trm_loss_names = array_to_percentiles(
                        trm_loss.to_numpy(),
                        options.transmission_loss_percentiles,
                        True)
                    
                    # store both in output dataset
                    sample_output['transmission loss'] = {
                        'percentiles': trm_loss_names,
                        'values': trm_loss_vals}
                    sample_output['transmission loss ion injection time scaled'] = {
                        'percentiles': trm_loss_ion_inj_names,
                        'values': trm_loss_ion_inj_vals}
            
                # count occurrences of each charge
                charges, counts = np.unique(ms2_df.loc[:, 'precursor charge'].to_numpy(),
                                            return_counts=True)
                # grab indices of charges sorted, to sort both lists in output
                idx_sort = np.argsort(charges)
                sample_output['charge dist'] = {'charge': charges.astype(str)[idx_sort].tolist(),
                                                'counts': counts[idx_sort].tolist()}
        
        statistics_dict['samples'][name] = sample_output
    
    
    db_search_norm_ms2_mean, db_search_norm_ms2_std, db_search_norm_x_vals = score_rank_dist_norm(
        file_loc_dict,
        'db search',
        options.db_search_format,
        'Confidence',
        normalize_ms2=True,
        ref_dict=statistics_dict
    )
    de_novo_norm_ms2_mean, de_novo_norm_ms2_std, de_novo_norm_x_vals = score_rank_dist_norm(
        file_loc_dict,
        'de novo',
        options.de_novo_format,
        'Confidence',
        normalize_ms2=True,
        ref_dict=statistics_dict
    )
    
    statistics_dict['global']['db search confidence dist norm ms2'] = {
        'mean': db_search_norm_ms2_mean.tolist(),
        'std': db_search_norm_ms2_std.tolist(),
        'x vals': db_search_norm_x_vals.tolist()
    }
    statistics_dict['global']['de novo confidence dist norm ms2'] = {
        'mean': de_novo_norm_ms2_mean.tolist(),
        'std': de_novo_norm_ms2_std.tolist(),
        'x vals': de_novo_norm_x_vals.tolist()
    }
    
    return statistics_dict
    

def write_reference_file(ref_dict: Dict[str, Dict[str, Any]],
                          filename: str, write_loc: Path):
    """Write spectral reference dataset to json file.

    Args:
        ref_dict (Dict[str, Dict[str, Any]]): Spectral reference dataset.
        filename (str): Name of json file
        write_loc (Path): Location in filesystem to store json.
    """
    def unexp_val_manager(val):
        if isinstance(val, (int, np.integer)):
            return int(val)
        if isinstance(val, float) and val != val:
            return None
        print(val, type(val))
        raise TypeError
    
    data_str = json.dumps(ref_dict,
                          default=unexp_val_manager)
    
    with open(Path(write_loc, filename.removesuffix(".json") + ".json"), 'w') as jsonfile:
        jsonfile.write(data_str)


def build_ref_data(
    data_root_dir: Path,
    output_loc: Path,
    file_name: str,
    db_search_format: DbSearchSource,
    de_novo_format: DeNovoSource,
    db_search_thresholds: List[float] = [30, 50, 80],
    de_novo_thresholds: List[float] = [50, 80, 90],
    intensity_percentiles: List[int | float] = [50, 90, 99],
    transmission_loss_percentiles: List[int | float] = [50, 90, 99],
    db_search_file_pattern: str | None = None,
    de_novo_file_pattern: str | None = None
    ) -> None:
    """Build experimental reference statistics dataset from collection of
    raw files, spectral files and metaproteomics data output files. By supplying
    a single root folder that contains all files (internal path structure does
    not matter) and statistics file build options, a single compact dataset will
    be constructed that contains experimental quality and performance metrics of
    all experimental datasets. This allows large experiment datasets 
    (>100 experiments) to be compressed into a small metrics file for quick 
    benchmarking.

    Args:
        data_root_dir (Path): Root directory containing all data files to be 
            included in the dataset.
        output_loc (Path): Path location to store reference statistics dataset.
        file_name (str): Name of reference statistics dataset.
        db_search_format (DbSearchSource): Format of db search data to collect.
            Exported file names should not be altered as files are collected
            by their expected file name formats (for example, Peaks 10 exports
            db search match data in "DB search psm.csv").
        de_novo_format (DeNovoSource): Format of de novo peptide identification 
            data to collect. Exported file names should not be altered as files 
            are collected by their expected file name formats (for example, 
            Peaks 10 exports de novo identification data in "de novo peptides.csv").
        db_search_thresholds (List[int  |  float], optional): Specify threshold
            values for quantifying db search identification confidence 
            distributions across experiments. Defaults to [30, 50, 80].
        de_novo_thresholds (List[int  |  float], optional): Specify threshold
            values for quantifying de novo identification confidence 
            distributions across experiments. Defaults to [50, 80, 90].
        intensity_percentiles (List[int  |  float], optional): Specify percentile
            values to compute summarized intensities accross experiments: "50" 
            computes the mean intensity, "90" computes the 90th percentile, etc.
            Defaults to [50, 90, 99].
        transmission_loss_percentiles (List[int  |  float], optional): Specify 
            percentile values to compute summarized signal transmission loss 
            accross experiments: "50" computes the mean transmission loss, "90" 
            computes the 90th percentile, etc. Defaults to [50, 90, 99].
        db_search_file_pattern (str | None, optional): Define custom regex 
            pattern to fetch db search files in the root directory. Can be used
            if filenames were changed, invalidating the default file match pattern,
            or for using formats that do not have a unique default pattern. 
            Defaults to None.
        de_novo_file_pattern (str | None, optional): Define custom regex 
            pattern to fetch de novo files in the root directory. Can be used
            if filenames were changed, invalidating the default file match pattern,
            or for using formats that do not have a unique default pattern. 
            Defaults to None.
    """
    
    # set defaults for file pattern to search through
    if db_search_file_pattern is None:
        db_search_file_pattern = db_search_file_name[db_search_format]
    if de_novo_file_pattern is None:
        de_novo_file_pattern = de_novo_file_name[de_novo_format]
    
    ref_options = RefBuilderOptions(
        root_dir=data_root_dir,
        db_search_format=db_search_format,
        de_novo_format=de_novo_format,
        db_search_thresholds=db_search_thresholds,
        de_novo_thresholds=de_novo_thresholds,
        intensity_percentiles=intensity_percentiles,
        transmission_loss_percentiles=transmission_loss_percentiles,
        db_search_file_pattern=db_search_file_pattern,
        de_novo_file_pattern=de_novo_file_pattern,
    )
    
    print("parse root directory for data files...")
    file_loc_dict = fetch_file_locations(ref_options)
    
    mzml_count = len([x["mzml"] for x in file_loc_dict.values() if x is not None])
    db_search_count = len([x["db search"] for x in file_loc_dict.values() if x['db search'] is not None])
    de_novo_count = len([x["de novo"] for x in file_loc_dict.values() if x['de novo'] is not None])

    print(f"found {len(file_loc_dict.keys())} experiments:")
    print(f"{mzml_count} mzML files")
    print(f"{db_search_count} DB search files")
    print(f"{de_novo_count} de novo files\n")
    

    print("Build reference dataset...")
    ref_data_dict = build_spectral_reference_data(file_loc_dict, ref_options)
    print("write reference dataset to json...")
    write_reference_file(ref_data_dict, file_name, output_loc)
    