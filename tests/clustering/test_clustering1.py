import numpy as np
from matplotlib import pyplot as plt # type: ignore

from mocca2 import Chromatogram
from mocca2.clustering.cluster_components import cluster_components
from mocca2.classes import Component, Compound, DeconvolvedPeak
from mocca2.math import cosine_similarity

from mocca2.deconvolution.fit_peak_model import fit_peak_model
from mocca2.deconvolution.peak_models import FraserSuzuki, BiGaussianTailing

print("Loading chromatograms")
chrom1 = Chromatogram('tests/test_data/09072021_istd_96.txt', 'tests/test_data/09072021_gradient_97.txt')
chrom2 = Chromatogram('tests/test_data/09072021_sample_4.txt', 'tests/test_data/09072021_gradient_97.txt')

assert chrom1.check_same_sampling(chrom2)

print("Correcting baseline")
chrom1.correct_baseline()
chrom2.correct_baseline()

print("Picking peaks")
chrom1.find_peaks('mean', 0.05)
chrom2.find_peaks('mean', 0.05)

print("Deconvolving peaks")
chrom1.deconvolve_peaks('BiGaussianTailing', 0.95, False, 3)
chrom2.deconvolve_peaks('BiGaussianTailing', 0.95, False, 3)

print("Clustering components")

components = chrom1.all_components() + chrom2.all_components()

print(f"There are {len(chrom1.peaks) + len(chrom2.peaks)} peaks with {len(components)} components.")

for p in chrom2.peaks:
    assert isinstance(p, DeconvolvedPeak)
    print(p.r2)

dt = chrom1.time_step()

def are_same(c1:Component, c2:Component)->bool:
    if np.abs(c1.elution_time - c2.elution_time) * dt > 0.1:
        return False
    if cosine_similarity(c1.spectrum, c2.spectrum) < 0.90:
        return False
    return True

compounds = cluster_components(components, are_same, lambda c: c.integral)

t = chrom1.time
print('Compounds:')
for idx, c in compounds.items():
    print(idx, t[c.elution_time], c.absorption_maxima())


plt.plot(chrom1.time, np.mean(chrom1.data, axis=0))
plt.plot(chrom2.time, np.mean(chrom2.data, axis=0))

plt.vlines([t[comp.elution_time] for comp in compounds.values()], 0, 50, colors='red', lw=0.5)

plt.show()