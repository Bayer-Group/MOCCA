import os
import numpy as np

from mocca2.classes import Data2D

def parse_chemstation(path) -> Data2D:
    """
    Chemstation read and processing function.
    """

    if os.path.isfile(path):
        raw = np.genfromtxt(path, delimiter=',', encoding='utf-16')
    else:
        raw = np.genfromtxt(os.path.join(path, 'DAD1.CSV'), delimiter=',', encoding='utf-16')


    time = raw[1:,0]
    wavelength = raw[0,1:]
    data = raw[1:,1:].T
    return Data2D(time, wavelength, data)
