from mocca2.deconvolution.peak_models import BiGaussian, BiGaussianTailing, FraserSuzuki, Bemg, test_gradients

import numpy as np

# Test gradients of model
assert test_gradients(BiGaussian(), np.linspace(15, 125), [2., 18., 0.5, 0.7])

assert test_gradients(BiGaussianTailing(), np.linspace(15, 125), [2., 18., 0.5, 0.7, 0.6, 0.3])

assert test_gradients(FraserSuzuki(), np.linspace(15, 125), [2., 18., 0.5, 0.7])

assert test_gradients(Bemg(), np.linspace(15, 125), [2., 18., 0.5, 0.7, 2.5])

print("Test finished")