import numpy as np
import pandas as pd


def dict_encoder(obj):
    """
    Custom JSON encoder for the MoccaDataset object.
    """
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, pd.DataFrame):
        return obj.to_dict(orient="list")
    elif np.issubdtype(type(obj), np.integer):
        return int(obj)
    elif np.issubdtype(type(obj), np.floating):
        return float(obj)
    elif np.issubdtype(type(obj), np.bool_):
        return bool(obj)
    elif isinstance(obj, dict):
        return {k: dict_encoder(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [dict_encoder(v) for v in obj]
    elif isinstance(obj, tuple):
        return tuple(dict_encoder(v) for v in obj)
    elif isinstance(obj, set):
        return set(dict_encoder(v) for v in obj)
    else:
        return obj
