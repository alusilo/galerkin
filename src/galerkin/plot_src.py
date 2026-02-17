import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import interpolate as it

data = np.load("../resources/output/initial_source_n3.npy")

fig = plt.figure()
ax2 = fig.add_subplot(121)
ax1 = fig.add_subplot(122, projection="3d")

xmin, xmax = np.min(data[0]), np.max(data[0])
ymin, ymax = np.min(data[1]), np.max(data[1])

points = np.array([data[0], data[1]]).T
X, Y = np.mgrid[xmin : xmax : 100 * 1j, ymin : ymax : 100 * 1j]
Z = it.griddata(points, data[2], (X, Y), method="linear", fill_value=0)

ax1.plot_trisurf(X.flatten(), Y.flatten(), Z.flatten(), cmap=cm.jet, linewidth=0, antialiased=False)
ax2.imshow(Z.T, extent=[xmin, xmax, ymin, ymax], aspect="equal", cmap=cm.jet, origin="upper")

plt.show()
