import os

from rainbow.waters import masslynx

from mocca2.classes import Data2D


def parse_masslynx(path: str) -> Data2D:
    """Reads MassLynx _FUNC*.DAT from the .raw directory"""

    if path.lower().endswith(".raw"):
        path = os.path.join(path, "_FUNC001.DAT")

    data = masslynx.parse_function(path)

    # masslinx reports absorbance in uAU, so conversion to mAU is needed
    return Data2D(data.xlabels, data.ylabels, data.data.T / 1000)
