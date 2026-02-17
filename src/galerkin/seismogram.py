import json
import os

import numpy as np
import matplotlib.pyplot as plt

width = 3.487
height = width / 1.618

plt.rc("font", family="serif", serif="Times")
plt.rc("text", usetex=True)
plt.rc("xtick", labelsize=8)
plt.rc("ytick", labelsize=8)
plt.rc("axes", labelsize=8)
plt.rc("legend", fontsize=6)
plt.rc("axes", titlesize=6)

# yFormatter = FormatStrFormatter('%.2f')

PROJECT_NAME = "gauss75_tests_s30/test1"
NDG_BINARY_FILE = "tracesVx.npy"
PARAM_FILENAME = "model.param"
CDG_BINARY_PATH = "../../DGCrack2D/INPLANE/demos/Point_source_h300m/VX"
DG_ROOT = os.path.dirname(os.getcwd())
LOCAL_PROJECT_DIR = os.path.join("resources/output", PROJECT_NAME)
PROJECT_DIR = os.path.join(DG_ROOT, LOCAL_PROJECT_DIR)
PARAM_FILE_PATH = os.path.join(PROJECT_DIR, PARAM_FILENAME)

FILE_PATH = os.path.join(PROJECT_DIR, NDG_BINARY_FILE)

param = json.load(open(PARAM_FILE_PATH))

data = np.load(FILE_PATH)
nf = len(data) - 1
Ndata2 = 2147
Nsteps = param["time_steps"]
dt = param["dt"]
time = data[0]
traces = data[1:]
ntraces = len(data[1:])
with open(CDG_BINARY_PATH, "rb") as f:
    data2 = np.fromfile(f, count=ntraces * Ndata2, dtype=np.float32).reshape(-1, ntraces).T
dgtime = np.linspace(0, time[-1], Nsteps)

time2 = np.array([i * 0.000932 for i in range(Ndata2)])

for i in range(nf):
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.16, bottom=0.17, right=0.98, top=0.92)
    ax.plot(time2, data2[i], c="black", linewidth=0.5, label="DGCrack")
    ax.plot(
        time,
        traces[i],
        "r--",
        c="black",
        linewidth=0.5,
        marker="s",
        markevery=100,
        markerfacecolor="none",
        markersize=3,
        label="nodalDG",
    )
    ax.set_ylabel("Amplitud")
    ax.set_xlabel(r"tiempo ($s$)")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    # ax.yaxis.set_major_formatter(yFormatter)
    fig.set_size_inches(width, height)
    plt.legend(loc=0)
    plt.savefig("test1_30{}.pdf".format(i))

# plt.show()
