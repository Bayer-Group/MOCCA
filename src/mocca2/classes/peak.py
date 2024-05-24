from __future__ import annotations
from typing import List, Any, Dict
from numpy.typing import NDArray

from mocca2.classes import Data2D
from mocca2.serializing import dict_encoder


class Peak:
    """Information about single peak"""

    left: int
    """Index of peak start"""
    right: int
    """Index of peak end"""
    maximum: int
    """Index of peak maximum"""
    height: float
    """Absolute height of the peak"""
    prominence: float
    """Height of the peak from the base"""
    all_maxima: List[int]
    """Indeces of all maxima of the merged peaks"""

    def __init__(
        self,
        left: int,
        right: int,
        maximum: int,
        height: float,
        prominence: float,
        all_maxima: List[int] | None = None,
    ):
        self.left = left
        self.right = right
        self.maximum = maximum
        self.height = height
        self.prominence = prominence
        self.all_maxima = all_maxima if all_maxima is not None else [maximum]

    def data(self, data: NDArray | Data2D) -> NDArray:
        """Returns part of the 2D data that contains this peak"""
        y = data.data if isinstance(data, Data2D) else data
        return y[:, self.left : self.right]

    def time(self, data: NDArray | Data2D) -> NDArray:
        """Returns part of the timescale that contains this peak"""
        y = data.time if isinstance(data, Data2D) else data
        return y[self.left : self.right]

    def to_dict(self) -> Dict[str, Any]:
        """Converts the data to a dictionary for serialization"""
        return dict_encoder(self.__dict__.copy() | {"__classname__": "Peak"})

    @staticmethod
    def from_dict(data: Dict[str, Any]) -> Peak:
        """Creates a Peak object from a dictionary"""
        assert data["__classname__"] == "Peak"

        return Peak(
            int(data["left"]),
            int(data["right"]),
            int(data["maximum"]),
            float(data["height"]),
            float(data["prominence"]),
            [int(d) for d in data["all_maxima"]],
        )
