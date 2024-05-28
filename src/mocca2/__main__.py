import argparse
import tarfile
import urllib.request
from importlib import resources as impresources
from mocca2 import example_data


def download_file(url, filename):
    urllib.request.urlretrieve(url, filename)


def extract_bz2(filename, path):
    with tarfile.open(filename, "r:bz2") as tar:
        tar.extractall(path)


def download_data():
    files = ["examples", "benzaldehyde", "diterpene_esters", "knoevenagel", "cyanation"]

    directory = impresources.files(example_data).joinpath("data")
    if directory.exists():
        raise ValueError("Data has been already downloaded.")
    directory.mkdir()

    for file in files:
        url = f"https://github.com/bayer-group/MOCCA/raw/example-data/{file}.tar.bz2"
        filename = directory.joinpath(f"{file}.tar.bz2")

        print(f"Downloading {file} data to {filename}")
        download_file(url, filename)
        print(f"Extracting {file} data")
        extract_bz2(filename, directory)

    print("Done!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MOCCA2 command line interface")
    parser.add_argument(
        "--download-data", help="Downloads example data", action="store_true"
    )
    args = parser.parse_args()
    if args.download_data:
        download_data()
    else:
        raise ValueError("No command provided")
