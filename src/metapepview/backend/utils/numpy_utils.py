import numpy as np
from typing import List, Tuple, Sequence


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

# definition for type hint
NumericArray = Sequence[int | float] | np.ndarray

class CubicSpline:
    """
    Lightweight cubic spline interpolator for 1D data.

    Behaves like scipy.interpolate.CubicSpline(x, y) for 1d arrays and default
    settings
    """

    def __init__(self, 
                 x: NumericArray, 
                 y: NumericArray, 
                 *, 
                 extrapolate: bool=True):
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.extrapolate = extrapolate

        if self.x.ndim != 1 or self.y.ndim != 1:
            raise ValueError("x and y must be 1-dimensional")
        if self.x.size != self.y.size:
            raise ValueError("x and y must have the same length")
        if self.x.size < 2:
            raise ValueError("At least 2 points are required")
        if np.any(~np.isfinite(self.x)) or np.any(~np.isfinite(self.y)):
            raise ValueError("x and y must contain only finite values")
        if np.any(np.diff(self.x) <= 0):
            raise ValueError("x must be strictly increasing")

        self.n = self.x.size
        self.h = np.diff(self.x)
        self.m = self._solve_second_derivatives()


    def _solve_second_derivatives(self) -> np.ndarray:
        """Solve for the second derivatives at the knots usingnot-a-knot boundary
        conditions.
        """
        x, y, h, n = self.x, self.y, self.h, self.n

        # Two points -> linear interpolation
        if n == 2:
            return np.array([0.0, 0.0], dtype=float)

        # Three points -> quadratic through the points
        # This matches the not-a-knot behavior for n=3.
        if n == 3:
            A = np.vstack([x**2, x, np.ones_like(x)]).T
            a, b, c = np.linalg.solve(A, y)
            return np.array([2.0 * a, 2.0 * a, 2.0 * a], dtype=float)

        A = np.zeros((n, n), dtype=float)
        rhs = np.zeros(n, dtype=float)

        # Not-a-knot at x# (m2 - m1) / h1 = (m1 - m0) / h0
        # => h0*m2 - (h0+h1)*m1 + h1*m0 = 0
        A[0, 0] = h[1]
        A[0, 1] = -(h[0] + h[1])
        A[0, 2] = h[0]

        # Interior equations:
        # h[i-1] * m[i-1] + 2(h[i-1] + h[i]) * m[i] + h[i] * m[i+1]
        #   = 6 * ((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
        for i in range(1, n - 1):
            A[i, i - 1] = h[i - 1]
            A[i, i] = 2.0 * (h[i - 1] + h[i])
            A[i, i + 1] = h[i]
            rhs[i] = 6.0 * (
                (y[i + 1] - y[i]) / h[i]
                - (y[i] - y[i - 1]) / h[i - 1]
            )

        # Not-a-knot at x[n-2]:
        # (m[n-1] - m[n-2]) / h[n-2] = (m[n-2] - m[n-3]) / h[n-3]
        # => h[n-2]*m[n-3] - (h[n-3]+h[n-2])*m[n-2] + h[n-3]*m[n-1] = 0
        A[n - 1, n - 3] = h[-1]
        A[n - 1, n - 2] = -(h[-2] + h[-1])
        A[n - 1, n - 1] = h[-2]

        return np.linalg.solve(A, rhs)

    def __call__(self, x_new: int | float | NumericArray) -> float | np.ndarray:
        """Calculate y values using cubic spline interpolation model from sample
        of x values.

        Args:
            x_new (int | float | NumericArray): Sample x values.

        Returns:
            float | np.ndarray: Calculated y values from cubic spline
                interpolation.
        """
        x_new = np.asarray(x_new, dtype=float)
        scalar_input = (x_new.ndim == 0)
        xq = np.atleast_1d(x_new)

        if self.extrapolate:
            idx = np.searchsorted(self.x, xq, side="right") - 1
            idx = np.clip(idx, 0, self.n - 2)
            yq = self._eval_segments(xq, idx)
        else:
            yq = np.full_like(xq, np.nan, dtype=float)
            valid = (xq >= self.x[0]) & (xq <= self.x[-1])
            if np.any(valid):
                idx = np.searchsorted(self.x, xq[valid], side="right") - 1
                idx = np.clip(idx, 0, self.n - 2)
                yq[valid] = self._eval_segments(xq[valid], idx)

        return yq.item() if scalar_input else yq

    def _eval_segments(self, 
                       xq: np.ndarray, 
                       idx: np.ndarray) -> np.ndarray:
        """Evaluate the spline on the selected intervals.
        """
        x, y, m = self.x, self.y, self.m

        h = x[idx + 1] - x[idx]
        a = (x[idx + 1] - xq) / h
        b = (xq - x[idx]) / h

        return (
            a * y[idx]
            + b * y[idx + 1]
            + ((a**3 - a) * m[idx] + (b**3 - b) * m[idx + 1]) * (h**2) / 6.0
        )