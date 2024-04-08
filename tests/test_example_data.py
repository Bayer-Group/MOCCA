from mocca2.example_data import loaders
from mocca2.classes.chromatogram import Chromatogram


def test_benzaldehyde():
    chroms = loaders.benzaldehyde()
    assert len(chroms) == 2


def test_examples():
    assert isinstance(loaders.example_1(), Chromatogram)
    assert isinstance(loaders.example_2(), Chromatogram)
    assert isinstance(loaders.example_3(), Chromatogram)


def test_knoevenagel():
    calib = loaders.knoevenagel_calibration()
    assert len(calib) == 120
    assert set("compound conc grad_len chromatogram".split()) == set(calib.columns)

    chroms = loaders.knoevenagel("ba_ome")
    assert len(chroms) == 25
    assert set("grad_len time chromatogram".split()) == set(chroms.columns)

    chroms = loaders.knoevenagel("ba_ome_nme2")
    assert len(chroms) == 50
    assert set("grad_len time chromatogram".split()) == set(chroms.columns)


def test_cyanation():
    data = loaders.cyanation()
    for (
        key
    ) in "istd educt_1 educt_2 product_1 product_2 cn_source_a cn_source_d".split():
        assert isinstance(data[key], Chromatogram)
    assert len(data["reactions"]) == 84


def test_diterpene():
    data = loaders.diterpene_esters()

    assert len(data) == 20
    assert set("sample chromatogram KO CO KP CP".split()) == set(data.columns)
