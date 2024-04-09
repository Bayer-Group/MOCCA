from mocca2 import load_data2d
from mocca2 import estimate_baseline

from time import time

from matplotlib import pyplot as plt # type: ignore

data = load_data2d('tests/test_data/raw_data3734.arw')
data -= load_data2d('tests/test_data/raw_data3725.arw')

start = time()
baseline = estimate_baseline(data, 'flatfit', 1, smooth_wl=10)
print("Elapsed", time()-start)

plt.imshow(baseline, aspect='auto', cmap='gist_ncar')

plt.figure()
plt.imshow(data.data - baseline, aspect='auto', cmap='gist_ncar')

plt.figure()
plt.plot(data.data[60])
plt.plot(baseline[60])
plt.show()