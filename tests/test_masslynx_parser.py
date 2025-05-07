from pathlib import Path

from mocca2.parsers import load_data2d
from mocca2.parsers.masslynx import parse_masslynx


def _test_data(path: list[str] = []) -> Path:
    return Path(__file__).resolve().parent.joinpath("test_data", *path)


def test_wrapper():
    path = str(_test_data(["masslynx_test.Raw"]))

    data = load_data2d(path)
    assert len(data.time) == 4801
    assert len(data.wavelength) == 191
    assert ((data.data > -500) & (data.data < 1500)).all()

    path = str(_test_data(["masslynx_test.Raw", "_FUNC001.DAT"]))

    data = load_data2d(path)
    assert len(data.time) == 4801
    assert len(data.wavelength) == 191
    assert ((data.data > -500) & (data.data < 1500)).all()


def test_loader():
    path = str(_test_data(["masslynx_test.Raw"]))

    data = parse_masslynx(path)
    assert len(data.time) == 4801
    assert len(data.wavelength) == 191
    assert ((data.data > -500) & (data.data < 1500)).all()

    path = str(_test_data(["masslynx_test.Raw", "_FUNC001.DAT"]))

    data = parse_masslynx(path)
    assert len(data.time) == 4801
    assert len(data.wavelength) == 191
    assert ((data.data > -500) & (data.data < 1500)).all()
