import argparse
from mocca2.example_data.downloader import download_data, unpack_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MOCCA2 command line interface")
    parser.add_argument(
        "--download-data", help="Downloads example data", action="store_true"
    )
    args = parser.parse_args()
    if args.download_data:
        download_data()
        unpack_data()
    else:
        raise ValueError("No command provided")
