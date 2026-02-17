import json
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as it
from scipy import linalg as la

from .nodalDG import Nodes2D, Vandermonde2D, xytors


def barycentric(vertices, p):
    # get vertices of the elements
    a = vertices[0] - vertices[1]
    b = vertices[2] - vertices[1]
    s = p - vertices[1]
    l1 = np.cross(s, a) / np.cross(b, a)
    l3 = np.cross(s, b) / np.cross(a, b)
    l2 = 1.0 - l1 - l3
    return l1 + 0, l2 + 0, l3 + 0


PROJECT_NAME = "gauss75_tests_s60"
PARAM_FILENAME = "model.param"
DG_ROOT = os.path.dirname(os.getcwd())
LOCAL_PROJECT_DIR = os.path.join("resources/output", PROJECT_NAME)
PROJECT_DIR = os.path.join(DG_ROOT, LOCAL_PROJECT_DIR)
PARAM_FILE_PATH = os.path.join(PROJECT_DIR, PARAM_FILENAME)

param = json.load(open(PARAM_FILE_PATH))
N = param["source"]["order"]
data = np.load("../resources/output/" + PROJECT_NAME + "/initial_source_n{}.npy".format(N))
# data = np.array([
# 	1.73443823E-05, 6.16954785E-05, 8.39711174E-06, 4.99628622E-05, 2.19158919E-05, 4.23596102E-05,
# 	1.35145492E-05, 3.62187461E-06, 4.81449661E-07, 1.41035634E-05, 3.82887083E-06, 2.06194431E-06,
# 	6.28286725E-05, 3.34721267E-06, 1.24896842E-05, 2.61333207E-05, 4.09113127E-05, 1.30340313E-05,
# 	8.87039914E-06, 3.29862758E-07, 1.83219654E-05, 2.36706546E-06, 2.31511403E-05, 5.53463769E-06,
# 	2.88116553E-07, 4.90284862E-08, 1.60032050E-05, 1.75086072E-07, 4.83419399E-06, 1.76904939E-06,
# 	3.36332027E-06, 6.31310249E-05, 9.07407320E-06, 2.62590802E-05, 9.87841941E-06, 3.65945489E-05,
# 	1.15718713E-05, 1.69176708E-08, 1.30413838E-07, 9.65411914E-07, 3.34046877E-06, 9.28157959E-08,
# 	1.20741088E-05, 1.36074007E-07, 6.17318074E-06, 3.48545041E-06, 2.00651357E-05, 1.62380672E-06,
# 	1.25628594E-05, 8.60144792E-06, 6.31967778E-05, 1.97166210E-05, 4.11510082E-05, 4.33903915E-05,
# 	3.32854347E-06, 8.98024700E-06, 1.02126158E-07, 9.77627678E-06, 1.00918407E-06, 1.52108976E-06,
# 	1.35803784E-05, 6.30445740E-08, 1.18290006E-06, 1.48211109E-06, 6.13032489E-06, 4.82913663E-07,
# 	1.88719059E-05, 9.64871924E-06, 6.71290109E-05, 3.13619312E-05, 5.43630995E-05, 3.89120250E-05,
# 	1.53929632E-05, 1.34078277E-06, 1.05391437E-05, 6.94854498E-06, 2.41582911E-05, 5.20341428E-06,
# 	1.40409502E-05, 5.00202475E-07, 6.51826966E-08, 3.97800795E-06, 1.53237613E-06, 3.08962541E-07])
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
# ax = fig.gca(projection='3d')
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
    # ax.plot_trisurf(PX,PY,PZ)

global_data = np.array([Xdata, Ydata, Zdata])
# Plot
xmin, xmax = np.min(global_data[0]), np.max(global_data[0])
ymin, ymax = np.min(global_data[1]), np.max(global_data[1])

points = np.array([global_data[0], global_data[1]]).T
X, Y = np.mgrid[xmin : xmax : 500 * 1j, ymin : ymax : 500 * 1j]
Z = it.griddata(points, global_data[2], (X, Y), method="linear", fill_value=0)

sx, sy = param["source"]["position"]
plt.xlabel(r"$x$ [$km$]")
plt.ylabel(r"$z$ [$km$]")
cmap = plt.cm.hot
# ax.plot_trisurf(Xdata,Ydata,Zdata,antialiased=True,linewidth=0,cmap=cmap)
plt.plot(
    sx / 1000,
    sy / 1000,
    markerfacecolor="yellow",
    fillstyle="full",
    markeredgecolor="black",
    marker="*",
    markersize=15,
)
# plt.scatter(Xdata/1000,Ydata/1000,color='black', marker='.', s=80)
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
