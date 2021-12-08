import numpy as np
from matplotlib import pyplot as plt
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
size = np.array([10, 20, 30, 40, 50])
nopart = np.array([910, 2810, 5710, 9610, 14510])
withpart = np.array([510, 1210, 2110, 3210, 4510])

fig = plt.figure()
ax = fig.gca()

plt.plot(size, nopart, label="Disable partition", marker="x")
plt.plot(size, withpart, label="Enable partition", marker="x")
plt.xlabel("Size")
plt.ylabel("Time/ns")
plt.grid()
plt.legend()
plt.savefig("partition.svg", format="svg", bbox_inches='tight')
plt.show()
