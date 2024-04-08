from mocca2.deconvolution.peak_models import Bemg
import numpy as np
from matplotlib import pyplot as plt

t = np.linspace(-10, 10, 500)
model = Bemg()

plt.plot(t, model(t, 1, 0, 1, 0.1, 0.1))
plt.plot(t, model(t, 1, 0, 1, 1., 1.))
plt.plot(t, model(t, 1, 0, 1, 4, 4))
plt.plot(t, model(t, 1, 0, 1, 4, 0.1))
plt.show()