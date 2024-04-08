import mocca2
from mocca2 import classes
from mocca2 import toy_data

from matplotlib import pyplot as plt

chrom = toy_data.benzaldehyde()

for c in chrom:
    plt.plot(c.time, c.contract())
plt.show()
