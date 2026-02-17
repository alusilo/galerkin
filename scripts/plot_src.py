#!/usr/bin/env python3
"""Plot initial source distribution (run from project root: uv run python scripts/plot_src.py)."""
import json
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as it

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT_NAME = "erase"
OUTPUT_DIR = os.path.join(_PROJECT_ROOT, "resources", "output", PROJECT_NAME)
PARAM_PATH = os.path.join(OUTPUT_DIR, "model.param")

# Source order N: from model.param if present, else default 3
if os.path.exists(PARAM_PATH):
    with open(PARAM_PATH) as f:
        N = json.load(f)["source"]["order"]
else:
    N = 3

data_path = os.path.join(OUTPUT_DIR, "initial_source_n{}.npy".format(N))
data = np.load(data_path)

fig = plt.figure()
ax2 = fig.add_subplot(121)
ax1 = fig.add_subplot(122, projection="3d")

xmin, xmax = np.min(data[0]), np.max(data[0])
ymin, ymax = np.min(data[1]), np.max(data[1])

points = np.array([data[0], data[1]]).T
X, Y = np.mgrid[xmin : xmax : 100 * 1j, ymin : ymax : 100 * 1j]
Z = it.griddata(points, data[2], (X, Y), method="linear", fill_value=0)

ax1.plot_surface(X, Y, Z, cmap="jet", linewidth=0, antialiased=False)
ax2.imshow(Z.T, extent=[xmin, xmax, ymin, ymax], aspect="equal", cmap="jet", origin="upper")

plt.show()
