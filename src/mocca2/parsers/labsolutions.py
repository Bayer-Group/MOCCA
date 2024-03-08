import numpy as np

from mocca2.classes import Data2D

def parse_labsolutions(path) -> Data2D:
    """Reads the .txt Lab Solutions file"""
    with open(path, "r") as f:
        text = f.read()
        start_idx = text.find("R.Time (min)")
        skip_lines = sum([ 1 for ch in text[:start_idx+1] if ch == '\n']) + 1

    raw = np.genfromtxt(path, delimiter=',', skip_header=skip_lines)

    time = raw[1:,0]
    wavelength = raw[0,1:]/100.
    data = raw[1:,1:].T/1000.
    return Data2D(time, wavelength, data)