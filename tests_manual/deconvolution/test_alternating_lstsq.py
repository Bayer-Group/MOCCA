from mocca2 import Chromatogram
from mocca2.deconvolution.alternating_lstsq import alternating_lstsq

import numpy as np

from matplotlib import pyplot as plt

chrom = Chromatogram('tests/test_data/09072021_sample_4.txt', 'tests/test_data/09072021_gradient_97.txt')

chrom.correct_baseline('flatfit', 1)

peak = chrom.extract_time(3.4, 3.7)

conc, spect, mse = alternating_lstsq(peak, 2)
print(mse/np.mean(peak.data**2))

plt.plot(peak.time, peak.contract())
for c in conc:
    plt.plot(peak.time, c)
plt.show()
