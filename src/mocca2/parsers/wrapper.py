from typing import Literal

from mocca2.classes import Data2D
from mocca2.parsers.empower import parse_empower
from mocca2.parsers.chemstation import parse_chemstation
from mocca2.parsers.labsolutions import parse_labsolutions

def load_data2d(path: str, format: Literal['auto', 'empower', 'chemstation', 'labsolutions'] = 'auto') -> Data2D:
    """
    Loads empower/chemstation/labsolutions file, returns 2D data

    Parameters
    ----------
    path: str
        Path to the raw data file

    format: Literal['auto', 'empower', 'chemstation', 'labsolutions']
        Format of the raw data file. If 'auto', the format is inferred from file extension

    Returns
    -------
    Data2D
        The 2D chromatogram data

    """

    if format not in {'auto', 'empower', 'chemstation', 'labsolutions'}:
        raise ValueError(f"Incorrect format `{format}` specified for load_data2d()")

    if format == 'auto':
        if path.lower().endswith('.arw'):
            data = parse_empower(path)
        elif path.lower().endswith('.csv') or path.lower().endswith('.d'):
            data = parse_chemstation(path)
        elif path.lower().endswith('.txt'):
            data = parse_labsolutions(path)
        else:
            raise Exception("Unknown file format in load_data2D(), consider specifying the format instead of using `auto`")
    elif format == 'empower':
        data = parse_empower(path)
    elif format == 'chemstation':
        data = parse_chemstation(path)
    elif format == 'labsolutions':
        data = parse_labsolutions(path)
        
    assert data.data.shape[0] == data.wavelength.shape[0], 'Parsing raw data by load_data2D() yields inconsistent shapes'
    assert data.data.shape[1] == data.time.shape[0], 'Parsing raw data by load_data2D() yields inconsistent shapes'

    return data