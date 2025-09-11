from typing import List, Dict, Literal, overload
from pathlib import Path
import pandas as pd
from scipy.interpolate import CubicSpline

from metapepview.backend.spectral_ref_builder.definitions import RefBuilderOptions
from metapepview.backend.type_operations import * # load_metapep_db_search, load_metapep_de_novo, db_search_importers, de_novo_importers
from metapepview.backend.utils import *


def ident_file_source(ident_file: Path,
                      filetype: str,
                      options: RefBuilderOptions) -> Sequence[str]:
    if filetype == "db search":
        source = db_search_importers[options.db_search_format]\
            .read_file(ident_file)\
            .get_source_files()
    elif filetype == "de novo":
        source = de_novo_importers[options.de_novo_format]\
            .read_file(ident_file)\
            .get_source_files()
    else:
        raise ValueError("Invlid filetype supplied")
    return source
    

def add_to_source_dict(output_dict: Dict[str, Dict[str, Path | None]],
                       sources: Sequence[str],
                       file: Path, 
                       filetype: str) -> Dict[str, Dict[str, Path | None]]:
    """Add data file location information to source files that correspond to 
    the data files. For example, a db search output file with three source files
    in the dataset will be added to each of the three source files in the output
    dictionary.

    Args:
        output_dict (Dict[str, Dict[str, Path  |  None]]): Source file to data 
            files mapping dictionary.
        sources (Sequence[str]): Source file keys to update with new file data.
        file (Path): Path to file to update source file keys with.
        filetype (str): Type of file data to add into the mapping dictionary.

    Returns:
        Dict[str, Dict[str, Path | None]]: Updated source file to data files
            mapping dictionary.
    """
    for source_name in sources:
        # locate raw file in dictionary
        raw_path = output_dict.get(source_name, None)

        if raw_path is None:
            output_dict[source_name] = {"raw": None,
                                        "db search": None,
                                        "de novo": None,
                                        "mzxml": None,
                                        "mzml": None}
        
        # add db search to output data
        output_dict[source_name][filetype] = file
        
    return output_dict


def count_threshold_values(
    df: pd.DataFrame,
    column_name: str,
    threshold_set: Sequence[int | float],
    cumulative: bool = True,
    div_factor: int | float | None = None) -> Tuple[List[int], List[int]]:
    """Generate lists of threshold counts for a given column name in the
    dataset, and group names. Thresholds are set by the user, as well as
    counting cumulative (Counting values towards all threshold groups it applies
    to) or not (Only counting values towards most stringent threshold group).

    Args:
        df (pd.DataFrame): Input dataset
        column_name (str): Column name from dataset to process
        threshold_set (Sequence[int]): Set of threshold values to create groups
            from. Any value higher than threshold is counted towards group.
        cumulative (bool, optional): Count value towards all groups for which
            threshold is met. If False, value is only counted towards most
            stringent . Defaults to True.
        div_factor (int | float | None, optional): Divide counts by a factor for
            normalization purposes. Defaults to False
    Raises:
        ValueError: Invalid column name

    Returns:
        Tuple[List[int], List[int]]: List of threshold counts, List of threshold
            group names.
    """
    
    # assert that column name is valid
    if not column_name in df.columns:
        raise ValueError("Invalid column name supplied")
    
    # return empty dataset if no threshold values given
    if len(threshold_set) == 0:
        return ([], [])
    
    # initialize output array
    counts_set = []
    count_names_set = []
    
    # sort given list of threshold values ascending
    threshold_set = sorted(threshold_set)
    
    # first, compute count below all threshold values, if not cumulative
    if cumulative is not True:
        counts_set.append((df[column_name] < threshold_set[0]).sum())
        count_names_set.append(f"<{threshold_set[0]}")
    
    # loop through threshold set to compute counts
    for i, val in enumerate(threshold_set):
        if cumulative is True or i == len(threshold_set) - 1:
            counts_set.append((df[column_name] > val).sum())
            count_names_set.append(f">{val}")
        else:
            counts_set.append(((df[column_name] > val) & 
                               (df[column_name] < threshold_set[i+1])).sum())
            count_names_set.append(f"{val} < x < {threshold_set[i+1]}")
    
    # divide counts by division factor if value given
    if div_factor is not None:
        counts_set = (np.array(counts_set) / div_factor).tolist()
    
    return (counts_set, count_names_set)


@overload
def score_rank_dist(data_dict: Dict[str, Dict[str, Path | None]],
                    file_type: Literal['db search'],
                    file_format: DbSearchSource,
                    score_col: str) -> Tuple[np.ndarray, np.ndarray]:
    ...

@overload
def score_rank_dist(data_dict: Dict[str, Dict[str, Path | None]],
                    file_type: Literal['de novo'],
                    file_format: DeNovoSource,
                    score_col: str) -> Tuple[np.ndarray, np.ndarray]:
    ...

def score_rank_dist(data_dict: Dict[str, Dict[str, Path | None]],
                    file_type: Literal['db search', 'de novo'],
                    file_format: DbSearchSource | DeNovoSource,
                    score_col: str) -> Tuple[np.ndarray, np.ndarray]:
    """Construct a distribution of confidence scores from a set of proteomics
    (db search, de novo) experiments by extracting the confidence column,
    sorting them descending and computing mean and standard deviation of the
    distribution. The resulting distribution plots the mean confidence of the
    n'th peptide (when sorted descending) across the set of experiments.

    Args:
        data_dict (Dict[str, Dict[str, Path | None]]): Dataset of file locations.
        file_type (Literal['db search', 'de novo']): File type to extract.
            Either db search or de novo data.
        file_format (DbSearchSource | DeNovoSource): Format of db search or 
            de novo data.
        score_col (str): Name of confidence column.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Score confidence distribution data,
            {mean values, std values}.
    """
    rank_mat = []
    max_len = 0

    for name, data in data_dict.items():
        df_path = data[file_type]
        if df_path is None:
            continue
        
        # extract confidence column from data
        if file_type == 'db search':
            df = load_metapep_db_search(df_path.open('r'), name, file_format) #type:ignore
        elif file_type == 'de novo':
            df = load_metapep_de_novo(df_path.open('r'), name, file_format) #type:ignore
        
        # if file contains data from other source files, omit them
        if len(df.source_files) > 1:
            df = df.filter_spectral_name(name)
        
        score_rank = fetch_sort_column(df.data, score_col)
        rank_mat.append(score_rank.values.tolist())
        max_len = max(max_len, score_rank.shape[0])

    # compute error bars
    rank_mat = np.array([x + (max_len - len(x))*[0] for x in rank_mat])
    mean_array = rank_mat.mean(axis=0)
    std_array = rank_mat.std(axis=0)
    
    return (mean_array, std_array)


def score_rank_dist_norm(data_dict: Dict[str, Dict[str, Path | None]],
                         file_type: Literal['db search', 'de novo'],
                         file_format: DbSearchSource | DeNovoSource,
                         score_col: str,
                         npoints: int = 1000,
                         normalize_ms2: bool = False,
                         ref_dict: dict | None = None,
                         fillna: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct a normalized distribution of confidence scores from a set of
    proteomics (db search, de novo) experiments by extracting the confidence
    column, sorting them descending and computing mean and standard deviation of
    the distribution. Normalization is performed over the number of peptides to
    equate the range from max confidence to threshold confidence between all
    experiments. The resulting figure therefore plots not the mean confidence
    of the n'th peptide when sorted, buth the mean n%'th fraction that has
    equal or higher confidence.

    Args:
        data_dict (Dict[str, Dict[str, Path | None]]): Dataset of file locations.
        file_type (Literal['db search', 'de novo']): File type to extract.
            Either db search or de novo data.
        file_format (DbSearchSource | DeNovoSource): Format of db search or 
            de novo data.
        score_col (str): Name of confidence column.
        npoints (int, optional): Size of the output array. Defaults to 1000.
        normalize_ms2 (bool, optional): If false, confidence distribution is
            normalized by total peptide matches. If true, the distribution is
            normalized by total ms2 scans. Defaults to False.
        ref_dict (dict | None, optional): Reference dataset being built. Is used
            to extract MS2 scan number from experiments if normalized by ms2.
            Defaults to None.
        fillna (bool, optional): If True, the fraction of scans/matches below
            the threshold is filled with nan values, else it is filled with
            0's. Defaults to False.
            
    Note:
        Not an exact method, it samples n evenly spaced points from
        the peptide match dataset of an experiment. Ideally, a non-linear curve
        fit is performed to interpolate confidence values at desired x values.

    Raises:
        ValueError: No reference dict provided while normalized by total MS2 scans.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: Score confidence distribution
            data, {mean values, std values, x values}.
    """
    
    
    # Create score distribution line plot normalized between 0-100 PSM / scans.
    
    # Note: Suboptimal method, it samples n evenly spaced points from
    # the peptide match datasets. Ideally, a non-linear curve fit function
    # is performed.
    rank_mat = []
    max_len = 0
    max_x_axis = []

    for name, data in data_dict.items():
        df_path = data[file_type]
        if df_path is None:
            continue
        
        # extract confidence column from data
        if file_type == 'db search':
            df = load_metapep_db_search(df_path.open('r'), name, file_format) #type:ignore
        elif file_type == 'de novo':
            df = load_metapep_de_novo(df_path.open('r'), name, file_format) #type:ignore
        
        # if file contains data from other source files, omit them
        if len(df.source_files) > 1:
            df = df.filter_spectral_name(name)
        
        # either normalize by total matches, or by ms2 count
        if normalize_ms2 is False:
            score_vals, x_axis = pept_match_dist_normalize(df.data,
                                                           score_col,
                                                           npoints,
                                                           None)
            max_x_axis = x_axis
        elif normalize_ms2 is True and ref_dict is not None:
            # get MS2 count, if no data, skip file
            total_scans = ref_dict["samples"][name]["ms2 count"]
            if total_scans != total_scans or total_scans is None:
                continue 
            score_vals, x_axis = pept_match_dist_normalize(df.data,
                                                           score_col,
                                                           npoints,
                                                           total_scans)
            max_len = max(max_len, score_vals.shape[0])
            if max_len == score_vals.shape[0]:
                max_x_axis = x_axis
        else:
            raise ValueError("Need reference dict to normalize ms2.")
            
        rank_mat.append(score_vals.tolist())

    # add nan or 0's to smaller lists to equate matrix dimensions
    if fillna is True:
        rank_mat = np.array([x + (max_len - len(x))*[0.0] for x in rank_mat])
    else:
        rank_mat = np.array([x + (max_len - len(x))*[np.nan] for x in rank_mat])
        
    rank_mat = np.array(rank_mat)
    # compute mean and error bars
    mean_array = np.nanmean(rank_mat, axis=0)
    std_array = np.nanstd(rank_mat, axis=0)
    
    return (mean_array, std_array, np.array(max_x_axis))


def pept_match_dist_normalize(df: pd.DataFrame,
                              score_col: str,
                              npoints: int=1000, 
                              total_scans: int | None=None) -> Tuple[np.ndarray, np.ndarray]:
    """Extract confidence scores from an identification dataset and transform 
    the data into a sorted, normalized and defined length data array.
    
    Confidence score range (x-axis) is normalized towards total MS2 scans. Thus,
    if fewer matches are present than scans, remaining values in the output
    are filled with zeros. However, if matches exceed the number of scans, extra
    datapoints are added to the output array.
    
    The length of the array is specified with 'npoints'. Here, the original
    array length extracted from the dataset (which has a length equal to the
    number of matches) is compressed into an array of 'npoints' number of values
    These values represent equally spaced confidence values in descending order 
    from the range of 0 to the total number of MS2 scans. An exeption is when the
    total number of matches exceed the total number of MS2 scans. In that case,
    the output array is extended in length such that the spacing of data points
    in the range remains equal.
    
    Note:
        Altering of the output array length is normally an elaborate and
        suboptimal way to normalize the distribution in the x-axis, as you could
        simply add x-values to the data in the original dimension size. However,
        in that case, multiple distribution arrays would have different x-axis
        datapoints. Determination of variances would therefore become a complex
        statistical task. By standardizing the range values towards equally
        spaced data points, variance computation can be done to a point-by-point
        basis.
        

    Args:
        df (pd.DataFrame): identification dataset.
        score_col (str): Column to extract confidence scores from.
        npoints (int, optional): Number of points in the output array
            (if less matches than total scans) between 0-1 (or 0% - 100%).
            If normalized by MS2, points are extended if valid confidence values
            extend beyond 100%. Defaults to 1000.
        total_scans (int | None, optional): Total number of MS2 scans.
            Defaults to None.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Confidence scores, x-axis range array
    """
    score_rank = fetch_sort_column(df, score_col)
    yvals = score_rank.values
    xvals = score_rank.index.to_numpy()
    
    # sample n samples evenly spaced over the complete dataset
    total_matches = score_rank.shape[0]
    
    # normalization range over all MS2 scans, points outside match range are zeros
    if total_scans is not None:
        # rescale x axis by dividing by total scans
        match_scan_ratio = total_matches / total_scans
        norm_xvals = xvals / total_scans
        
        spline_model = CubicSpline(norm_xvals, yvals)
        
        # define consistent spacing based on original defined npoints
        spacing = 1 / npoints
        # normalization towards MS2 does not cap at 100%. To compute mean + std,
        # rank steps are kept equal, but samples above 100 get more datapoints
        sample_x = np.arange(0, match_scan_ratio, spacing)
        sample_y = spline_model(sample_x)

        return (sample_y, sample_x * 100) # convert x to %

    else:
        # rescale x axis by dividing by total matches
        norm_xvals = xvals / (total_matches - 1)
        
        spline_model = CubicSpline(norm_xvals, yvals)
        
        # define consistent spacing based on original defined npoints
        spacing = 1 / npoints
        sample_x = np.arange(0, 1, spacing)
        sample_y = spline_model(sample_x)
        
        return (sample_y, sample_x * 100) # convert x to %


def transmission_loss(prec_int_array: pd.Series,
                      ms2_tic_array: pd.Series,
                      prec_injection_time_array: pd.Series | None = None,
                      ms2_injection_time_array: pd.Series | None = None,
                      invert_comp: bool = False) -> pd.Series:
    """Compute transmission loss for selected peptide signals. Transmission loss
    is computed by dividing total ion current of the MS2 scan with the MS2 precursor
    signal. Computation can also be inverted to compute transmission efficiency.

    Args:
        prec_int_array (pd.Series): Series of precursor intensities.
        ms2_tic_array (pd.Series): Series of corresponding MS2 total signal
            intensities.
        prec_injection_time_array (pd.Series | None, optional): Series of
            injection times for precursor spectra. Defaults to None.
        ms2_injection_time_array (pd.Series | None, optional): Series of
            injection times for ms2 spectra. Defaults to None.
        invert_comp (bool, optional): Invert division to compute transmission
            efficiency. Defaults to False.

    Returns:
        pd.Series: Series of transmission losses.
    """
    # divide signal intensities by injection times to take into account accumulation time
    if prec_injection_time_array is not None and ms2_injection_time_array is not None:
        prec_int_array = prec_int_array / prec_injection_time_array
        ms2_tic_array = ms2_tic_array / ms2_injection_time_array

    # divide precursor with ms2 signal, or vice versa if invert_comp
    if invert_comp is False:
        return prec_int_array / ms2_tic_array
    else:
        return ms2_tic_array / prec_int_array

def fetch_precursor_ion_injection_time(ms2_df: pd.DataFrame,
                                       ms1_df: pd.DataFrame) -> pd.Series:
    """Fetch series of precursor ion injection times based on the precursor 
    scan numbers present in the MS2 dataset. These scan numbers are used to
    match MS1 scans of the precursor to the MS2 rows and fetches the ion
    injection time.

    Args:
        ms2_df (pd.DataFrame): Dataframe of MS2 scans.
        ms1_df (pd.DataFrame): Dataframe of MS1 scans.

    Returns:
        pd.Series: Series of precursor injection times based on ms2 df rows.
    """
    # fetch precursor scan number column as separate df
    prec_scan_nums = ms2_df[['precursor scan number']]
    
    # merge ms1 data with precursor scan number df and fetch ion injection time
    prec_inj_time = prec_scan_nums.merge(ms1_df[['scan number', 'ion injection time']],
                                        left_on='precursor scan number',
                                        right_on='scan number',
                                        how='left')['ion injection time']
    return prec_inj_time

def array_to_percentiles(input_array: np.ndarray,
                         percentiles: List[int | float],
                         ignore_nan: bool = False) -> Tuple[List[float],
                                                            List[str]]:
    """Compute percentile distributions from array of input values and user
    supplied percentile values.

    Args:
        input_array (np.ndarray): Array of input data to process.
        percentiles (List[int  |  float]): List of percentiles to calculate.
        ignore_nan (bool, optional): If true, percentile function will ignore
            nan values for processing.

    Returns:
        Tuple[List[float], List[str]]: Tuple with list of percentile values, and 
            percentile names.
    """
    out_names, out_vals = [], []
    
    for pval in percentiles:
        if pval == 50:
            name = 'median'
        else:
            name = f'{pval}th'
        
        if len(input_array) == 0:
            scan_int_pt = np.nan
        elif ignore_nan is True:
            scan_int_pt = np.nanpercentile(input_array, pval)
        else:
            scan_int_pt = np.percentile(input_array, pval)

        out_names.append(name)
        out_vals.append(scan_int_pt)
        
    return (out_vals, out_names)
    

# function to compute user specified percentiles
def scan_intensity_percentiles(ms_df: pd.DataFrame,
                               scan_int_col: str,
                               percentiles: List[int | float]) -> Tuple[List[float],
                                                                         List[str]]:
    """Extract ion intensity distribution information from mzml data by
    calculating intensities at a specific number of percentile values.

    Args:
        ms_df (pd.DataFrame): mzml data in dataframe format, should include
            intensity column.
        scan_int_col (str): Name of intensity column
        percentiles (List[int  |  float]): List of percentiles to calculate.

    Returns:
        Tuple[List[str], List[float]]: List of percentile names for graph display,
        and list of percentile values.
    """
    int_array = ms_df.loc[:, scan_int_col].to_numpy()
    return array_to_percentiles(int_array, percentiles)


def calculate_miscleavages(db_search: pd.DataFrame) -> Tuple[List[str], List[int]]:
    """Group peptide matches by the number of miscleavages and calculate
    the groupsize for each miscleavage. Returns two lists: the list of 
    miscleavage groups, and the groups sizes.

    Args:
        db_search (pd.DataFrame): DB search dataset

    Returns:
        Tuple[List[str], List[int]]: Distribution of miscleavages. The first list are
        miscleavage groups, the second list are group sizes
    """
    categories = ["0", "1", "2", "3", ">3"]

    # count all occurences of cleave amino acids not at the end of the sequence
    # by removing last character. Take into account potential empty strings or
    # nan values
    miscleave_series: pd.Series = db_search['Sequence']\
        .apply(lambda x: x if x != x or x == "" else x[:-1])\
        .str.count("[KR]")\
        .value_counts()
    
    # convert categories to string
    miscleave_series.index = [str(x) for x in miscleave_series.index]

    # store all cleavage categories in counts array
    counts = list()
    for cat in categories[:-1]:
        counts.append(miscleave_series.get(cat, 0))

    # final counts are all values outside of defined categories
    counts.append(miscleave_series.drop(categories[:-1]).sum())

    return (categories, counts)
