import mocca2, mocca2.example_data, mocca2.classes.chromatogram
import numpy


def deep_equal(a, b):
    if type(a) != type(b):
        return False
    if isinstance(a, dict):
        if len(a) != len(b):
            return False
        for key in a:
            if key not in b:
                return False
            if not deep_equal(a[key], b[key]):
                return False
        return True
    if isinstance(a, list):
        if len(a) != len(b):
            return False
        for i in range(len(a)):
            if not deep_equal(a[i], b[i]):
                return False
        return True
    if isinstance(a, tuple):
        if len(a) != len(b):
            return False
        for i in range(len(a)):
            if not deep_equal(a[i], b[i]):
                return False
        return True
    if (
        isinstance(a, float)
        or isinstance(a, int)
        or isinstance(a, str)
        or isinstance(a, bool)
        or isinstance(a, type(None))
    ):
        return a == b

    if isinstance(a, numpy.ndarray):
        return numpy.allclose(a, b)

    return deep_equal(a.__dict__, b.__dict__)


def test_dict_serialization():
    chrom = mocca2.example_data.example_1()
    chrom.time = chrom.time[::5]
    chrom.data = chrom.data[:, ::5]
    chrom.extract_wavelength(220, None, inplace=True)
    chrom.extract_time(1.4, 1.8, inplace=True)
    chrom.find_peaks(min_height=2)
    chrom.deconvolve_peaks(
        model="FraserSuzuki", min_r2=0.99, relaxe_concs=False, max_comps=4
    )
    chrom_dict = chrom.to_dict()
    chrom2 = mocca2.classes.chromatogram.Chromatogram.from_dict(chrom_dict)

    print(chrom_dict)
    print(chrom2.to_dict())

    assert deep_equal(chrom, chrom2)
    assert deep_equal(chrom_dict, chrom2.to_dict())
