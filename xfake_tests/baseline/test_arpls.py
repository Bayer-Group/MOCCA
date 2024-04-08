import numpy as np
from matplotlib import pyplot as plt # type: ignore

from mocca2.baseline.arpls import arpls


def gen_data(only_baseline=False):
    def gaussian(x, mu, sigma):
        return np.exp(-(x-mu)**2 / 2 / sigma**2)
    
    x = np.linspace(0, 10, 1000)
    
    peaks = 0.5 * gaussian(x, 2, 0.3) + 2.0 * gaussian(x, 3, 0.1) + 0.3 * gaussian(x, 3.5, 0.1) + 1.5 * gaussian(x, 6, 0.2) + 1.0 * gaussian(x, 7, 0.5)

    low_freq_noise = 0.3 * gaussian(x, 4, 8) + 0.5 * gaussian(x, 8, 6)

    high_freq_noise = np.random.randn(x.shape[0]) * 0.05

    if only_baseline:
        return x, low_freq_noise
    return x, peaks + low_freq_noise + high_freq_noise

x, signal = gen_data()
_, true_baseline = gen_data(only_baseline=True)
baseline = arpls(signal, 2., max_iter=1000)[0]
baseline4 = arpls(signal[::4], 2.)[0]
baseline4s = arpls(signal[::4]*4, 2.)[0]/4

plt.plot(x, signal)
plt.plot(x, true_baseline)
#plt.plot(x, baseline)
plt.plot(x[::4], baseline4)
plt.plot(x[::4], baseline4s - baseline4)
plt.show()