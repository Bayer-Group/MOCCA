from typing import Literal

from importlib import resources as impresources
import os
import re
import numpy as np
import pandas as pd

from mocca2.classes.chromatogram import Chromatogram, Data2D
from mocca2 import toy_data


def example_1() -> Chromatogram:
    """
    Loads example chromatogram.

    Key features are large peak at 2.4 min and three overlapping peaks around 1.5 min.
    """

    chrom = impresources.files(toy_data) / "data/examples/chrom1.arw"
    blank = impresources.files(toy_data) / "data/examples/blank1.arw"

    return Chromatogram(str(chrom), str(blank), name="example_chromatrogram_1")


def example_2() -> Chromatogram:
    """
    Loads example chromatogram.

    Rather complicated chromatogram with many peaks around 1.5 - 2.5 min, as well as strong peak at 0.2 min and broad peak at 3.5 min.
    """

    chrom = impresources.files(toy_data) / "data/examples/chrom2.arw"
    blank = impresources.files(toy_data) / "data/examples/blank2.arw"

    return Chromatogram(str(chrom), str(blank), name="example_chromatrogram_2")


def example_3() -> Chromatogram:
    """
    Loads example chromatogram.

    Simple chromatogram with 4 pure peaks in the first 1 minute.
    """

    chrom = impresources.files(toy_data) / "data/examples/chrom3.arw"
    blank = impresources.files(toy_data) / "data/examples/blank3.arw"

    return Chromatogram(str(chrom), str(blank), name="example_chromatrogram_3")


def knoevenagel_calibration(which: Literal["1", "2", "both"] = "both") -> pd.DataFrame:
    """
    Loads all chromatograms used for calibration curves in the Knoevenagel reaction, published in 10.1021/acscentsci.2c01042.

    The columns are:
    - compound: compound name (ba, ome, nme2)
    - conc: concentration of the compound (in mM)
    - grad_len: gradient length (in minutes)
    - chromatogram: the chromatogram object
    """
    if which == "both":
        return pd.concat(
            [knoevenagel_calibration(which="1"), knoevenagel_calibration(which="2")]
        )
    if which not in ["1", "2"]:
        raise ValueError("which must be one of '1', '2' or 'both'")
    directory = str(
        impresources.files(toy_data) / ("data/knoevenagel/calibration" + which)
    )

    compound_pattern = re.compile(
        r".*_(?P<grad_len>\d+)_(?P<compound>\w+)_(?P<conc>\d+)\.D"
    )

    blank_pattern = re.compile(r".*_(?P<grad_len>\d+)_grad\.D")

    # parse file names
    compound_data = []
    blank_data = []
    for file in os.listdir(directory):
        match = compound_pattern.match(file)
        if match:
            compound_data.append(match.groupdict() | {"file": file})
            continue
        match = blank_pattern.match(file)
        if match:
            blank_data.append(match.groupdict() | {"file": file})
            continue

    time_interpolations = {
        "250": np.linspace(0, 3, 1800),
        "150": np.linspace(0, 2, 1200),
        "100": np.linspace(0, 1.5, 900),
        "075": np.linspace(0, 1.5, 900),
        "050": np.linspace(0, 1.5, 900),
    }

    # load the chromatograms
    for idx, data in enumerate(compound_data):
        chrom = Chromatogram(os.path.join(directory, data["file"]))
        chrom.interpolate_time(
            time_interpolations[data["grad_len"]], kind="cubic", inplace=True
        )
        compound_data[idx]["chrom"] = chrom

    for idx, data in enumerate(blank_data):
        chrom = Chromatogram(os.path.join(directory, data["file"]))
        chrom.interpolate_time(
            time_interpolations[data["grad_len"]], kind="cubic", inplace=True
        )
        blank_data[idx]["chrom"] = chrom

    # average the blank chromatograms
    blanks = {}
    for grad_len in set([data["grad_len"] for data in blank_data]):
        chroms = [data["chrom"] for data in blank_data if data["grad_len"] == grad_len]
        blank = Data2D(
            chroms[0].time,
            chroms[0].wavelength,
            np.mean([chrom.data for chrom in chroms], axis=0),
        )
        blanks[int(grad_len)] = blank

    # substract blanks
    for idx, data in enumerate(compound_data):
        data["chrom"] -= blanks[int(data["grad_len"])]

    # create the dataframe
    df = pd.DataFrame(compound_data)
    df = df.rename(columns={"chrom": "chromatogram"})
    df["grad_len"] = df["grad_len"].astype(float).map(lambda x: f"{x/100:0.2f}")
    df.drop(columns=["file"], inplace=True)

    # correct the concentrations
    # calibration standard concentrations
    concs = {
        "ba": {"100": 0.9947, "075": 0.7618, "050": 0.4841, "025": 0.2486},
        "ome": {"100": 1.3162, "075": 0.9911, "050": 0.6379, "025": 0.2968},
        "nme2": {"100": 0.9676, "075": 0.7219, "050": 0.4753, "025": 0.2483},
    }
    for idx, row in df.iterrows():
        row["conc"] = concs[row["compound"]][row["conc"]]
    df["conc"] = df["conc"].astype(float)

    return df


def knoevenagel(which: Literal["ba_ome", "ba_ome_nme2"]) -> pd.DataFrame:
    """
    Loads all chromatograms from Knoevenagel reaction with benzaldehyde, 4-methoxybenzaldehyde and optionally 4-(N,N-dimethyl)benzaldehyde, published in 10.1021/acscentsci.2c01042.

    The columns are:
    - grad_len: gradient length (in minutes)
    - time: reaction time in minutes
    - chromatogram: the chromatogram object
    """

    if which not in ["ba_ome", "ba_ome_nme2"]:
        raise ValueError("which must be one of 'ba_ome', 'ba_ome_nme2'")

    directory = str(impresources.files(toy_data) / f"data/knoevenagel/reaction_{which}")

    time_interpolations = {
        "250": np.linspace(0, 3, 1800),
        "150": np.linspace(0, 2, 1200),
        "100": np.linspace(0, 1.5, 900),
        "075": np.linspace(0, 1.5, 900),
        "050": np.linspace(0, 1.5, 900),
    }

    chromatogram_pattern = re.compile(
        r".*(?P<hrs>\d\d)-(?P<min>\d\d)-(?P<sec>\d\d)_(?P<grad_len>\d+)_ba_ome.*\.D"
    )
    blank_pattern = re.compile(r".*_(?P<grad_len>\d+)_grad\.D")

    # parse file names
    chromatogram_data = []
    blank_data = []
    for file in os.listdir(directory):
        match = chromatogram_pattern.match(file)
        if match:
            chromatogram_data.append(match.groupdict() | {"file": file})
            continue
        match = blank_pattern.match(file)
        if match:
            blank_data.append(match.groupdict() | {"file": file})
            continue

    # load the chromatograms
    for idx, data in enumerate(chromatogram_data):
        chrom = Chromatogram(os.path.join(directory, data["file"]))
        chrom.interpolate_time(
            time_interpolations[data["grad_len"]], kind="cubic", inplace=True
        )
        chromatogram_data[idx]["chrom"] = chrom

    for idx, data in enumerate(blank_data):
        chrom = Chromatogram(os.path.join(directory, data["file"]))
        chrom.interpolate_time(
            time_interpolations[data["grad_len"]], kind="cubic", inplace=True
        )
        blank_data[idx]["chrom"] = chrom

    # average the blank chromatograms
    blanks = {}
    for grad_len in set([data["grad_len"] for data in blank_data]):
        chroms = [data["chrom"] for data in blank_data if data["grad_len"] == grad_len]
        blank = Data2D(
            chroms[0].time,
            chroms[0].wavelength,
            np.mean([chrom.data for chrom in chroms], axis=0),
        )
        blanks[int(grad_len)] = blank

    # substract blanks
    for idx, data in enumerate(chromatogram_data):
        data["chrom"] -= blanks[int(data["grad_len"])]

    # create the dataframe
    df = pd.DataFrame(chromatogram_data)
    df = df.rename(columns={"chrom": "chromatogram"})
    df["grad_len"] = df["grad_len"].astype(float).map(lambda x: f"{x/100:0.2f}")
    df.drop(columns=["file"], inplace=True)

    # correct the time
    df["time"] = (
        df["hrs"].astype(int) * 60 + df["min"].astype(int) + df["sec"].astype(int) / 60
    )
    df.drop(columns=["hrs", "min", "sec"], inplace=True)

    if which == "ba_ome":
        start_time = 13 * 60 + 35 + 20 / 60
    elif which == "ba_ome_nme2":
        start_time = 15 * 60 + 53 + 30 / 60
    df["time"] -= start_time

    return df
