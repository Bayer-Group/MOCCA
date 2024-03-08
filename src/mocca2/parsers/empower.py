import numpy as np
import pandas as pd

from mocca2.classes import Data2D

def parse_empower(path) -> Data2D:
    """Reads the .arw empower file"""
    with open(path, "r") as f:
        skip_lines = 0
        for line in f:
            if "Wavelength" in line:
                wavelength = np.array(line.split('\t')[1:]).astype(np.float32)
            if "Time" in line:
                break
            skip_lines += 1

    skiprows = list(range(skip_lines + 1))

    raw = pd.read_csv(
        path,
        delimiter='\t',
        engine='c',
        header=None,
        skiprows=skiprows,
        dtype=np.float32
    ).to_numpy(np.float32)

    time = raw[:,0]
    data = raw[:,1:].T  * 1000.
    return Data2D(time, wavelength, data)
