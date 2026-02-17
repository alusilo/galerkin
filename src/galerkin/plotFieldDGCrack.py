import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

from .Mesh2D import *  # noqa: F403

fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)

x_min = 0.0
x_max = 5600.0
y_min = 0.0
y_max = 3600.0
pix = 10.0
dt = 0.01
x_pix = int(np.ceil((x_max - x_min) / pix + 1))
y_pix = int(np.ceil((y_max - y_min) / pix + 1))
snap_pix = x_pix * y_pix
n_snaps = 100
data = []
PROJECT_NAME = "gauss75_tl_p1p2"
with open("../../DGCrack2D/INPLANE/demos/Point_source_h300m/VX.snap", "rb") as f:
    data = np.fromfile(f, count=snap_pix * n_snaps, dtype=np.float32).reshape(-1, snap_pix)

print("total frames: {}".format(n_snaps))

frame = data[0].reshape(-1, y_pix).T

vmin, vmax = data.min() / 100, data.max() / 100

im = plt.imshow(
    frame,
    extent=[x_min, x_max, y_max, y_min],
    aspect="equal",
    cmap=plt.get_cmap("jet"),
    origin="upper",
    animated=True,
)
im.set_clim([vmax, vmin])
title = ax.set_title("Wave propagation {}[s]".format(0.0))
ax.set_xlabel("Distance [m]")
ax.set_ylabel("Depth [m]")
bar = plt.colorbar()

pause = False


def simData():
    frame = data[0].reshape(-1, y_pix).T
    i = 0
    while i < n_snaps:
        if not pause:
            frame = data[i].reshape(-1, y_pix).T
            i = i + 1
        yield frame, i


def onClick(event):
    global pause
    pause ^= True


def updatefig(simData):
    frame, i = simData
    im.set_array(frame)
    ax.set_title("Wave propagation (t={0:.3f}[s])".format((i) * dt))
    return (im,)


fig.canvas.mpl_connect("button_press_event", onClick)
ani = animation.FuncAnimation(
    fig, updatefig, simData, interval=50, repeat_delay=400, blit=False, repeat=True
)
plt.show()
