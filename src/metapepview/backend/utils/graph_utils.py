import math
import numpy as np
import pandas as pd

from typing import Tuple


def set_tickvals(upper_limit: float, maxticks: int=5) -> np.ndarray:
    """Generate set of tick values to display on the y-axis of a graph.
    
    If standard generated tick values from graphing libraries do not give
    the desired result. This function can be used to create a set of values
    with defined upper limit and a maximum of tick values. The function
    generates tick values rounded out as much as possible.

    Note:
        It is designed for use in Plotly. This library accepts an array
        ov custom tick values for display in the figure. Other graphing
        libraries may not support custom tick values.

    Args:
        upper_limit (float): Upper limit of tick values.
        maxticks (int, optional): Maximum number of values generated.
            Defaults to 5.

    Returns:
        np.ndarray: Array of tick values
    """
    # specify tickvalues for second y axis (required to remove 0 tick)
    decimals = math.floor(math.log10(upper_limit))
    power_decimal = 10**decimals

    first_digit_ceil = math.ceil(upper_limit / power_decimal)
    
    # if upper limit is of form 1xxx..., round above for second digit
    if first_digit_ceil == 2 and decimals > 0:
        new_p_dec = 10**(decimals-1)
        second_digit_ceil = math.ceil(upper_limit / new_p_dec)
        round_maxval = second_digit_ceil * new_p_dec
    else:    
        round_maxval = first_digit_ceil * power_decimal

    # decide on spacing based on first digit value and number of ticks
    # we multiply 10**decimals by 5 as that is the average of the range 1-10
    # then we divide by maxticks as more ticks mean smaller spacing
    spacing = power_decimal * 5 / (10**math.floor(math.log10(maxticks)))
    
    
    nticks = math.floor((round_maxval+spacing) / spacing)
    
    # update spacing if needed
    new_spacing = spacing

    # decrease spacing as long as max ticks not reached
    div_gaps = iter([2, 5, 10])
    div_num = next(div_gaps)
    last_div_num = 1
    
    # check if there is potential to decrease spacing
    while nticks <= maxticks/ (div_num / last_div_num):
        new_spacing = spacing / div_num
        nticks = math.floor((round_maxval+new_spacing) / new_spacing)
        
        # update array variables
        last_div_num = div_num
        div_num=next(div_gaps)
    
    # when maxticks is exceeded, increase spacing by 2
    mult_val = 2
    while nticks > maxticks:
        new_spacing = spacing * mult_val
        mult_val += 1
        nticks = math.floor((round_maxval+new_spacing) / new_spacing)
    
    spacing = new_spacing

    return np.arange(0, round_maxval+spacing, spacing)


def standardize_array_length(sorted_scores: pd.Series,
                             npoints: int=1000, 
                             total_scans: int | None=None) -> Tuple[np.ndarray, np.ndarray]:
    """Extract confidence scores from an identification dataset and normalize 
    values within distribution in defined range and array length.
    
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
            (if less matches than total scans). Defaults to 1000.
        total_scans (int | None, optional): Total number of MS2 scans.
            Defaults to None.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Confidence scores, x-axis range array
    """
    # sample n samples evenly spaced over the complete dataset
    total_matches = sorted_scores.shape[0]
    
    # normalization range over all MS2 scans, points outside match range are zeros
    if total_scans is not None:
        npoints_in_dataset = int(npoints * (total_matches / total_scans))
        score_vals = sorted_scores.values[np.linspace(0, total_matches-1, npoints_in_dataset, dtype=int).tolist()]
        
        # normalization towards MS2 does not cap at 100%. To compute mean + std,
        # rank steps are kept equal, but samples above 100 get more datapoints
        return (score_vals, np.linspace(0,
                                        100 * (total_matches / total_scans),
                                        npoints_in_dataset))
        
    else:
        score_vals = sorted_scores.values[np.linspace(0, total_matches-1, npoints, dtype=int).tolist()]
    
        return (score_vals, np.linspace(0, 100, npoints))
 