import tarfile
import urllib.request
from importlib import resources as impresources
from pathlib import Path
from mocca2 import example_data

_FILES = ["examples", "benzaldehyde", "diterpene_esters", "knoevenagel", "cyanation"]


def _example_data_path(path: list[str] = []) -> Path:
    directory = Path(str(impresources.files(example_data))).resolve() / "data"
    return directory.joinpath(*path)


def _download_file(url, filename):
    urllib.request.urlretrieve(url, filename)


def _extract_bz2(filename, path):
    with tarfile.open(filename, "r:bz2") as tar:
        tar.extractall(path)


def download_data(verbose: bool = True):
    """Downloads .tar.bz2 archives containing example data"""
    directory = _example_data_path()

    directory.mkdir(exist_ok=True)

    for file in _FILES:
        url = f"https://github.com/bayer-group/MOCCA/raw/example-data/{file}.tar.bz2"
        filename = _example_data_path([f"{file}.tar.bz2"])

        if filename.is_file():
            if verbose:
                print(f"File {filename} already downloaded, skipping...")
            continue

        if verbose:
            print(f"Downloading {file} data to {filename}")
        _download_file(url, filename)

    if verbose:
        print("Done!")


def unpack_data(verbose: bool = True):
    """Unpacks .tar.bz2 archives containing example data"""

    directory = _example_data_path()

    if not directory.is_dir():
        raise ValueError("Data has not been downloaded yet.")

    for file in _FILES:
        filename = _example_data_path([f"{file}.tar.bz2"])
        if not filename.is_file():
            raise RuntimeError(
                f"Could not find archive {filename}, make sure it has been downloaded"
            )

        if verbose:
            print(f"Extracting {file} data")
        _extract_bz2(filename, directory)

    if verbose:
        print("Done!")
