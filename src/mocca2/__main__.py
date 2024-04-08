import argparse
import tarfile
import urllib.request


def download_file(url, filename):
    urllib.request.urlretrieve(url, filename)


def extract_bz2(filename, path):
    with tarfile.open(filename, "r:bz2") as tar:
        tar.extractall(path)


def download_data():
    # Your code here
    print("This is my custom command!")


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
