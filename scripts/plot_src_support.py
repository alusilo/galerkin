#!/usr/bin/env python3
"""Plot spatial source support (run from project root: uv run python scripts/plot_src_support.py)."""
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as it
from scipy import linalg as la

# Import from galerkin package (run from project root with uv run)
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_src = os.path.join(_PROJECT_ROOT, "src")
if _src not in sys.path:
    sys.path.insert(0, _src)
from galerkin.nodalDG import Nodes2D, Vandermonde2D, xytors


def barycentric(vertices, p):
    a = vertices[0] - vertices[1]
    b = vertices[2] - vertices[1]
    s = p - vertices[1]
    l1 = np.cross(s, a) / np.cross(b, a)
    l3 = np.cross(s, b) / np.cross(a, b)
    l2 = 1.0 - l1 - l3
    return l1 + 0, l2 + 0, l3 + 0


PROJECT_NAME = "erase"
PARAM_FILENAME = "model.param"
LOCAL_PROJECT_DIR = os.path.join("resources", "output", PROJECT_NAME)
PROJECT_DIR = os.path.join(_PROJECT_ROOT, LOCAL_PROJECT_DIR)
PARAM_FILE_PATH = os.path.join(PROJECT_DIR, PARAM_FILENAME)

param = json.load(open(PARAM_FILE_PATH))
N = param["source"]["order"]
data = np.load(os.path.join(PROJECT_DIR, "initial_source_n{}.npy".format(N)))

elements = param["source"]["elements"]
Nfp = N + 1
Np = int(Nfp * (Nfp + 1) / 2)

samples = 50
dr = 2.0 / samples

X, Y = Nodes2D(N)
R, S = xytors(X, Y)
V = Vandermonde2D(N, R, S)
invV = la.inv(V)

PR = []
PS = []
for i in range(samples):
    s = samples - i
    for j in range(s):
        PR.append(-1 + j * dr)
        PS.append(-1 + i * dr)

PR = np.array(PR)
PS = np.array(PS)

phi_n = Vandermonde2D(N, PR, PS)
interp = np.dot(phi_n, invV)

Xdata = np.array([])
Xdata.shape = (0, Np)
Ydata = np.array([])
Ydata.shape = (0, Np)
Zdata = np.array([])
Zdata.shape = (0, Np)

fig = plt.figure()
for k in range(elements):
    local_data = data[:, k * Np : (k + 1) * Np]
    mapV = [0, Nfp - 1, Np - 1]
    vertices = local_data[:2, mapV].T
    x, y, z = local_data
    PX, PY = 0.5 * (
        -np.dot(vertices[0].reshape(-1, 1), [PR + PS])
        + np.dot(vertices[1].reshape(-1, 1), [1 + PR])
        + np.dot(vertices[2].reshape(-1, 1), [1 + PS])
    )
    PZ = np.dot(interp, z)
    Xdata = np.append(Xdata, PX)
    Ydata = np.append(Ydata, PY)
    Zdata = np.append(Zdata, PZ)

global_data = np.array([Xdata, Ydata, Zdata])
xmin, xmax = np.min(global_data[0]), np.max(global_data[0])
ymin, ymax = np.min(global_data[1]), np.max(global_data[1])

points = np.array([global_data[0], global_data[1]]).T
X, Y = np.mgrid[xmin : xmax : 500 * 1j, ymin : ymax : 500 * 1j]
Z = it.griddata(points, global_data[2], (X, Y), method="linear", fill_value=0)

sx, sy = param["source"]["position"]
plt.xlabel(r"$x$ [$km$]")
plt.ylabel(r"$z$ [$km$]")
cmap = plt.cm.hot
plt.plot(
    sx / 1000,
    sy / 1000,
    markerfacecolor="yellow",
    fillstyle="full",
    markeredgecolor="black",
    marker="*",
    markersize=15,
)
plt.imshow(
    Z.T[::-1, :],
    extent=[xmin / 1000, xmax / 1000, ymin / 1000, ymax / 1000],
    aspect="equal",
    cmap=cmap,
    vmin=0,
)
bar = plt.colorbar()
plt.gca().invert_yaxis()
plt.show()
