import argparse


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
