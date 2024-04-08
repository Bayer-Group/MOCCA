import mocca2
from mocca2 import classes
from mocca2.toy_data import loaders

from matplotlib import pyplot as plt

chrom = loaders.knoevenagel("ba_ome_nme2")

plt.plot(chrom["time"], chrom["grad_len"].astype(float), "x")
plt.show()
