from mocca2 import Chromatogram
from mocca2.deconvolution.fit_peak_model import fit_peak_model
from mocca2.deconvolution.nonnegative_lstsq import spectra_from_concentrations
from mocca2.deconvolution.peak_models import Bemg
from matplotlib import pyplot as plt
import numpy as np

chrom = Chromatogram('tests/test_data/raw_data1052.arw', 'tests/test_data/blank2.arw')

chrom.correct_baseline()
chrom.extract_wavelength(220, 400, inplace=True)

peak = chrom.extract_time(2.1, 2.3)

conc, mse, _ = fit_peak_model(peak.data, Bemg(), n_compounds=3)
spectra, _ = spectra_from_concentrations(peak.data, conc)
norms = np.mean(spectra, axis=1)
conc = (conc.T*norms).T
spectra = (spectra.T/norms).T

# for s in spectra:
#     plt.plot(peak.wavelength, s)
# plt.show()

for c in conc:
    plt.plot(peak.time, c)

plt.plot(peak.time, peak.contract())
plt.show()