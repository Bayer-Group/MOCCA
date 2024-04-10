from typing import Literal, Dict, List, Tuple

from importlib import resources as impresources
import os
import re
import numpy as np
import pandas as pd
import scipy.io

from mocca2.classes.chromatogram import Chromatogram, Data2D
from mocca2 import example_data


def check_data_needs_downloading():
    """
    Check if the example data needs to be downloaded.
    """
    if not os.path.exists(impresources.files(example_data) / "data"):
        raise FileNotFoundError(
            "Example data not found. Please run 'python -m mocca2 --download-data' to download the example data."
        )


def example_1(substract_blank: bool = True) -> Chromatogram:
    """
    Loads example chromatogram.

    Key features are large peak at 2.4 min and three overlapping peaks around 1.5 min.
    """
    check_data_needs_downloading()

    chrom = impresources.files(example_data) / "data/examples/chrom1.arw"
    if not substract_blank:
        return Chromatogram(str(chrom), name="example_chromatrogram_1")

    blank = impresources.files(example_data) / "data/examples/blank1.arw"
    return Chromatogram(str(chrom), str(blank), name="example_chromatrogram_1")


def example_2(substract_blank: bool = True) -> Chromatogram:
    """
    Loads example chromatogram.

    Rather complicated chromatogram with many peaks around 1.5 - 2.5 min, as well as strong peak at 0.2 min and broad peak at 3.5 min.
    """
    check_data_needs_downloading()

    chrom = impresources.files(example_data) / "data/examples/chrom2.arw"
    if not substract_blank:
        return Chromatogram(str(chrom), name="example_chromatrogram_2")

    blank = impresources.files(example_data) / "data/examples/blank2.arw"
    return Chromatogram(str(chrom), str(blank), name="example_chromatrogram_2")


def example_3(substract_blank: bool = True) -> Chromatogram:
    """
    Loads example chromatogram.

    Simple chromatogram with 4 pure peaks in the first 1 minute.
    """
    check_data_needs_downloading()

    chrom = impresources.files(example_data) / "data/examples/chrom3.arw"
    if not substract_blank:
        return Chromatogram(str(chrom), name="example_chromatrogram_3")

    blank = impresources.files(example_data) / "data/examples/blank3.arw"
    return Chromatogram(str(chrom), str(blank), name="example_chromatrogram_3")


def knoevenagel_calibration(
    which: Literal["1", "2", "both"] = "both", substract_blank: bool = True
) -> pd.DataFrame:
    """
    Loads all chromatograms used for calibration curves in the Knoevenagel reaction, published in 10.1021/acscentsci.2c01042.

    The columns are:
    - compound: compound name (ba, ome, nme2)
    - conc: concentration of the compound (in mM)
    - grad_len: gradient length (in minutes)
    - chromatogram: the chromatogram object
    """
    check_data_needs_downloading()

    if which == "both":
        return pd.concat(
            [knoevenagel_calibration(which="1"), knoevenagel_calibration(which="2")]
        )
    if which not in ["1", "2"]:
        raise ValueError("which must be one of '1', '2' or 'both'")
    directory = str(
        impresources.files(example_data) / ("data/knoevenagel/calibration" + which)
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
    if substract_blank:
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


def knoevenagel(
    which: Literal["ba_ome", "ba_ome_nme2"], substract_blank: bool = True
) -> pd.DataFrame:
    """
    Loads all chromatograms from Knoevenagel reaction with benzaldehyde, 4-methoxybenzaldehyde and optionally 4-(N,N-dimethyl)benzaldehyde, published in 10.1021/acscentsci.2c01042.

    The columns are:
    - grad_len: gradient length (in minutes)
    - time: reaction time in minutes
    - chromatogram: the chromatogram object
    """
    check_data_needs_downloading()

    if which not in ["ba_ome", "ba_ome_nme2"]:
        raise ValueError("which must be one of 'ba_ome', 'ba_ome_nme2'")

    directory = str(
        impresources.files(example_data) / f"data/knoevenagel/reaction_{which}"
    )

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
    if substract_blank:
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


def cyanation(
    substract_blank: bool = True,
) -> Dict[str, Chromatogram | List[Chromatogram]]:
    """
    Loads the chromatograms from the cyanation reaction, published in 10.1021/acscentsci.2c01042.

    The entries are:
    - istd: chromatogram with just the internal standard
    - educt_1 and educt_2: standards of the starting material
    - product_1 and product_2: standards of the product
    - cn_source_a and cn_source_d: standards of the cyanation source a nad d (standards for other sources are not included)
    - reactions: list of 96 all reactions
    """
    check_data_needs_downloading()

    directory = str(impresources.files(example_data) / "data/cyanation")

    def filename(name: str) -> str:
        return os.path.join(directory, f"09072021_{name}.txt")

    # load all the chromatograms
    chromatograms = {}

    if substract_blank:
        blank = Chromatogram(filename("gradient_97"))
    else:
        blank = None

    chromatograms["istd"] = Chromatogram(filename("istd_96"), blank, name="istd")
    chromatograms["educt_1"] = Chromatogram(filename("educt_88"), blank, name="educt_1")
    chromatograms["educt_2"] = Chromatogram(filename("educt_89"), blank, name="educt_2")
    chromatograms["product_1"] = Chromatogram(
        filename("product_92"), blank, name="product_1"
    )
    chromatograms["product_2"] = Chromatogram(
        filename("product_93"), blank, name="product_2"
    )
    chromatograms["cn_source_a"] = Chromatogram(
        filename("cnsource_a_98"), blank, name="cn_source_a"
    )
    chromatograms["cn_source_d"] = Chromatogram(
        filename("cnsource_d_99"), blank, name="cn_source_d"
    )

    chromatograms["reactions"] = []
    for i in range(84):
        chromatograms["reactions"].append(
            Chromatogram(
                filename(f"sample_{i+4:d}"),
                blank,
                name=f"reaction_{i+1:d}",
            )
        )

    return chromatograms


def benzaldehyde(substract_blank: bool = True) -> Tuple[Chromatogram, Chromatogram]:
    """Loads tutorial data published with the original MOCCA package, these contain 1mM and 0.5mM benzaldehyde respectively."""
    check_data_needs_downloading()

    directory = str(impresources.files(example_data) / "data/benzaldehyde")

    if substract_blank:
        blank = Chromatogram(os.path.join(directory, "blank.D"))
    else:
        blank = None

    chrom_1 = Chromatogram(
        os.path.join(directory, "ba_1.D"), blank, interpolate_blank=True
    )
    chrom_2 = Chromatogram(
        os.path.join(directory, "ba_05.D"), blank, interpolate_blank=True
    )

    return chrom_1, chrom_2


def diterpene_esters() -> pd.DataFrame:
    """
    Loads the chromatograms of diterpene esters from  data from 10.5281/zenodo.5412345.

    The data includes both calibration chromatograms with known concentrations and samples with coffee extracts.

    The columns are:
    - sample: sample name
    - KO, CO, KP, CP: known concentrations of the compounds
    - chromatogram: Data2D objects with chromatograms
    """
    check_data_needs_downloading()

    directory = str(impresources.files(example_data) / "data/diterpene_esters")

    mat = scipy.io.loadmat(os.path.join(directory, "data.mat"))

    # load chromatogram data
    time = mat["Data"][0, 0][0][:, 0].astype(float)
    wavelength = mat["Data"][0, 0][1][0].astype(float)

    # sample, wavelength, time
    calib = mat["Data"][0, 0][2].transpose([2, 1, 0]).astype(float)
    samples = mat["Data"][0, 0][4].transpose([2, 1, 0]).astype(float)

    # load sample names and concentrations
    calib_names = [
        "2-2_Rep1",
        "2-2_Rep2",
        "5-5_Rep1",
        "5-5_Rep1",
        "10-10_Rep1",
        "10-10_Rep2",
        "20-20_Rep1",
        "20-20_Rep2",
        "50-50_Rep1",
        "50-50_Rep2",
        "75-100_Rep1",
        "75-100_Rep2",
        "100-150_Rep1",
        "100-150_Rep2",
        "150-200_Rep1",
        "150-200_Rep2",
    ]
    sample_names = [
        "Boiled coffee - 1-Rep1",
        "Boiled coffee - 1-Rep2",
        "Boiled coffee - 2-Rep1",
        "Boiled coffee - 2-Rep2",
    ]
    calib_concs = pd.read_csv(os.path.join(directory, "calib.csv"))

    # prepare chromatograms
    calib = [Data2D(time, wavelength, c) for c in calib]
    samples = [Data2D(time, wavelength, c) for c in samples]

    # create the dataframe
    df = pd.DataFrame(
        [
            {
                "sample": sample_names[i],
                "KO": None,
                "CO": None,
                "KP": None,
                "CP": None,
                "chromatogram": sample,
            }
            for i, sample in enumerate(samples)
        ]
        + [
            {
                "sample": calib_names[i],
                "KO": calib_concs.loc[i, "KO"],
                "CO": calib_concs.loc[i, "CO"],
                "KP": calib_concs.loc[i, "KP"],
                "CP": calib_concs.loc[i, "CP"],
                "chromatogram": calib,
            }
            for i, calib in enumerate(calib)
        ]
    )

    return df
