import numpy as np
from typing import List, Tuple


def match_mz_rt_peaks(mz_array_a: np.ndarray,
                      mz_array_b: np.ndarray,
                      rt_array_a: np.ndarray,
                      rt_array_b: np.ndarray,
                      charge_array_a: np.ndarray,
                      charge_array_b: np.ndarray,
                      mz_threshold: float = 0.01,
                      rt_threshold: float = 0.5) -> List[Tuple[int, int]]:
    """
    Finds index pairs (i, j) where mz and RT arrays are both below user defined
    thresholds and charge values are identical. Used to match features to 
    db search or de novo data.

    Note:
        It is recommended that the largest array is assigned array a. This way,
        batch processing is most memory efficient.

    Args:
        mz_array_a (np.ndarray): First m/z array.
        mz_array_b (np.ndarray): Second m/z array.
        rt_array_a (np.ndarray): First RT array, same shape as mz_array_a.
        rt_array_b (np.ndarray): Second RT array, same shape as mz_array_b.
        charge_array_a (np.ndarray): First charge array, same shape as mz_array_a.
        charge_array_b (np.ndarray): Second charge array, same shape as mz_array_b.
        mz_threshold (float, optional): m/z threshold value. Defaults to 0.01.
        rt_threshold (float, optional): RT threshold value. Defaults to 0.5.
        batch_size (int, optional): Elements to process in parallel.Defaults to 100.

    Returns:
        List[Tuple[int, int]]: List of (i, j) index pairs
    """
    pairs = []

    mz_b_sort_idx = np.argsort(mz_array_b)
    mz_arr_b_sorted = mz_array_b[mz_b_sort_idx]

    rt_b_sort_idx = np.argsort(rt_array_b)
    rt_arr_b_sorted = rt_array_b[rt_b_sort_idx]

    # parse b array with binary search to find matching
    for idx, a_elem in enumerate(mz_array_a):
        # Find indices in sorted B where B is within (a - t, a + t)
        mz_low = np.searchsorted(mz_arr_b_sorted, a_elem - mz_threshold, side="left")
        mz_high = np.searchsorted(mz_arr_b_sorted, a_elem + mz_threshold, side="right")

        rt_low = np.searchsorted(rt_arr_b_sorted, rt_array_a[idx] - rt_threshold, side="left")
        rt_high = np.searchsorted(rt_arr_b_sorted, rt_array_a[idx] + rt_threshold, side="right")

        # Map back to original indices of B
        matching_b_idx_mz = mz_b_sort_idx[mz_low:mz_high]
        matching_b_idx_rt = rt_b_sort_idx[rt_low:rt_high]

        for j in np.intersect1d(matching_b_idx_mz, matching_b_idx_rt):
            # only add if charges are correct
            if charge_array_b[j] == charge_array_a[idx]:
                pairs.append((idx, j))
    return pairs
