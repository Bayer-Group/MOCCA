"""These routines can contract 2D chromatogram to create 1D data suitable for peak picking"""

from typing import Literal
from numpy.typing import NDArray, ArrayLike

import numpy as np
        
def cosine_similarity(x: NDArray, y: NDArray, axis: int = -1) -> float | NDArray:
    """Calculates cosine similarity over specified axis"""

    dot_prod = np.sum(x * y, axis=axis)
    norm_x = np.linalg.norm(x, axis=axis)
    norm_y = np.linalg.norm(y, axis=axis)

    norm_x = np.clip(norm_x, 1e-5, None)
    norm_y = np.clip(norm_y, 1e-5, None)

    return dot_prod/norm_x/norm_y
