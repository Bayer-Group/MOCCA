from __future__ import annotations
from typing import Tuple, Literal
from numpy.typing import NDArray
from copy import deepcopy

import numpy as np
from scipy.interpolate import interp1d


class Data2D:
    """2D chromatogram data"""

    time: NDArray
    """Time points at which data was sampled"""
    wavelength: NDArray
    """Wavelengths at which data was sampled"""
    data: NDArray
    """Absorbances at given wavelength and time `absorbance[wavelength, time]`"""

    def __init__(self, time: NDArray, wavelength: NDArray, data: NDArray):
        self.time = time
        self.wavelength = wavelength
        self.data = data

    def closest_time(self, time: float) -> Tuple[int, float]:
        """Returns index and value of time point that is closest to specified `time`"""
        return _closest(self.time, time)

    def closest_wavelength(self, wavelength: float) -> Tuple[int, float]:
        """Returns index and value of wavelength point that is closest to specified `wavelength`"""
        return _closest(self.wavelength, wavelength)

    def extract_time(self, min_time: float, max_time: float, inplace: bool = False) -> Data2D:
        """
        Extracts the data in the given time range

        Parameters
        ----------

        min_time: float
            Start time of the extracted segment

        max_time: float
            End time of the extracted segment

        inplace: bool
            If True, modifies the data in-place and returns self

        Returns 
        -------

        Data2D
            The data in the given time interval

        """

        start_idx = self.closest_time(min_time)[0]
        end_idx = self.closest_time(max_time)[0]

        if inplace:
            self.time=self.time[start_idx:end_idx+1]
            self.data=self.data[:, start_idx:end_idx+1]
            return self

        return Data2D(
            time=self.time[start_idx:end_idx+1],
            wavelength=self.wavelength,
            data=self.data[:, start_idx:end_idx+1]
        )
    
    def extract_wavelength(self, min_wavelength: float, max_wavelength: float, inplace: bool = False) -> Data2D:
        """
        Extracts the data in the given wavelength range

        Parameters
        ----------

        min_wavelength: float
            Start wavelength of the extracted segment

        max_wavelength: float
            End wavelength of the extracted segment

        inplace: bool
            If True, modifies the data in-place and returns self

        Returns 
        -------

        Data2D
            The data in the given wavelength interval

        """

        start_idx = self.closest_wavelength(min_wavelength)[0]
        end_idx = self.closest_wavelength(max_wavelength)[0]

        if inplace:
            self.wavelength=self.wavelength[start_idx:end_idx+1]
            self.data=self.data[start_idx:end_idx+1]
            return self

        return Data2D(
            wavelength=self.wavelength[start_idx:end_idx+1],
            time=self.time,
            data=self.data[start_idx:end_idx+1]
        )

    def check_same_sampling(self, *others: Data2D, tol: float = 1e-3, time: bool = True, wavelength: bool = True) -> bool:
        """Checks that the time sampling and wavelength sampling is same as in all provided data"""

        for o in others:
            if time and self.time.shape != o.time.shape:
                return False
            if wavelength and self.wavelength.shape != o.wavelength.shape:
                return False
            
            if time:
                time_tol = tol * np.sqrt(np.mean(self.time**2))
            if wavelength:
                wl_tol = tol * np.sqrt(np.mean(self.wavelength**2))

            if time and not np.allclose(self.time, o.time, atol=time_tol):
                return False
            if wavelength and not np.allclose(self.wavelength, o.wavelength, atol=wl_tol):
                return False
           
        return True
    
    def interpolate_time(self, time: NDArray, kind: str = 'linear', inplace:bool = False) -> Data2D:
        """Interpolates the data to the given time points using specified interpolation, see scipy.interpolate.interp1d for details"""
        new_data = np.empty((self.wavelength.shape[0], time.shape[0]))
        
        for wl_index in range(self.wavelength.shape[0]):
            new_data[wl_index] = interp1d(self.time, self.data[wl_index], kind, bounds_error=False, fill_value=(self.data[wl_index,0], self.data[wl_index, -1]))(time)

        if inplace:
            self.time = time
            self.data = new_data
            return self

        return Data2D(time, self.wavelength, new_data)

    def time_step(self) -> float:
        """Returns the sampling step of the time axis"""
        return self.time[1] - self.time[0]

    def wavelength_step(self) -> float:
        """Returns the sampling step of the wavelength axis"""
        return self.wavelength[1] - self.wavelength[0]

    def contract(
        self,
        method: Literal['mean', 'max', 'weighted_mean'] = 'mean',
        damping: float = 0.2
    ) -> NDArray:
        """
        Contracts the first dimension of 2D data to get 1D data for peak picking

        Parameters
        ----------
        method: Literal['mean', 'max', 'weighted_mean']
            The method that should be used for contraction. 'weighted_mean' weights the wavelenghts average std

        damping: float
            damping factor for 'weighted_mean'

        Returns
        -------
        Contracted 1D array
        """

        y = self.data

        if method not in {'mean', 'max', 'weighted_mean'}:
            raise ValueError("Invalid method for Data2D.contract()")

        if method == 'mean':
            return np.mean(y, axis=0)

        elif method == 'max':
            return np.max(y, axis=0)

        elif method == 'weighted_mean':
            overall_std = np.std(y)
            wl_std = np.std(y, axis=1)
            weights = 1. / ((1-damping) * wl_std + damping * overall_std)
            weights /= np.sum(weights)
            return y.T @ weights

    # Methods for adding and subtracting absorbace data

    def __add__(self, other: Data2D | NDArray) -> Data2D:
        if isinstance(other, Data2D):
            assert self.check_same_sampling(
                other), "The data has different sampling"
            return Data2D(self.time, self.wavelength, self.data + other.data)

        return Data2D(self.time, self.wavelength, self.data + other)

    def __sub__(self, other: Data2D | NDArray) -> Data2D:
        if isinstance(other, Data2D):
            assert self.check_same_sampling(
                other), "The data has different sampling"
            return Data2D(self.time, self.wavelength, self.data - other.data)

        return Data2D(self.time, self.wavelength, self.data - other)

    def to_json(self):
        json_dict = deepcopy(self.__dict__)
        for key in json_dict.keys():
            json_dict[key] = json_dict[key].tolist()
        return json_dict

    def from_json(json_dict):
        for key in json_dict.keys():
            json_dict[key] = np.array(json_dict[key])
        data2d = Data2D(**json_dict)
        return data2d

def _closest(data: NDArray, point: float) -> Tuple[int, float]:
    """
    Returns the index and value of entry in `data` closest to `point`

    Time complexity is O(n), but using numpy.
    """
    if point == np.inf:
        idx = np.argmax(data)
    elif point == -np.inf:
        idx = np.argmin(data)
    else:
        idx = np.argmin(np.abs(data-point))
    return int(idx), data[idx]
