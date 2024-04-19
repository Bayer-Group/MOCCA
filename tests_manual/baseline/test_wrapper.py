import numpy as np
from matplotlib import pyplot as plt # type: ignore

from mocca2.baseline import estimate_baseline


def gen_data(only_baseline=False):
    def gaussian(x, mu, sigma):
        return np.exp(-(x-mu)**2 / 2 / sigma**2)
    
    x = np.linspace(0, 10, 1000)
    
    peaks = 0.5 * gaussian(x, 2, 0.3) + 2.0 * gaussian(x, 3, 0.1) + 0.3 * gaussian(x, 3.5, 0.1) + 1.5 * gaussian(x, 6, 0.2) + 1.0 * gaussian(x, 7, 0.5)

    low_freq_noise = 0.3 * gaussian(x, 4, 8) + 0.5 * gaussian(x, 8, 6)

    high_freq_noise = np.random.randn(x.shape[0]) * 0.02

    if only_baseline:
        return x, low_freq_noise
    return x, peaks + low_freq_noise + high_freq_noise

x, signal = gen_data()
_, true_baseline = gen_data(only_baseline=True)
b1 = estimate_baseline(signal, method='asls', smoothness=10., p=0.0001)
b2 = estimate_baseline(signal, method='arpls', smoothness=10.)
b3 = estimate_baseline(signal, method='flatfit', smoothness=10.)

lots_of_data = np.array([gen_data()[1] * (np.sin(i)+1.2) for i in np.linspace(0,20,300)])
lots_of_baselines = estimate_baseline(lots_of_data)

print(b1.shape, b2.shape, x.shape)

plt.plot(x, signal)
plt.plot(x, true_baseline)
plt.plot(x, b1)
plt.plot(x, b2)
plt.plot(x, b3)
# plt.plot(x, lots_of_baselines[10])
plt.show()