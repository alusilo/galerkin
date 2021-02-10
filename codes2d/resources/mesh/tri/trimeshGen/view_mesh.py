import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
import matplotlib.tri as tri

import os, json

from Mesh2D import *

fig, ax = plt.subplots(1, 1, sharex=False, sharey=False, figsize=(14,7))
src_position = (2138.1344775252105, 382.34938629690299)
MESH_FILE = '../tlr.ele'

obj = MeshReader(
	mesh_file 	 = MESH_FILE,
	src_position = (0, 0),
	gather 		 = [(0, 0)],
	src_smooth 	 = 0,
	plot = False
)

K, Nv, VX, VY, e2v = obj.get_values()

#title = ax.set_title("Malla Triangular no Estructurada")
ax.set_xlabel(r'$x$ [$km$]')
ax.set_ylabel(r'$z$ [$km$]')
sx,sy = src_position
plt.scatter(sx/1000,sy/1000,color='red', marker='*', s=80)
mesh = tri.Triangulation(VX/1000,VY/1000)
plt.gca().invert_yaxis()
plt.triplot(mesh, lw=0.5, color='black')
ax.set_aspect('equal')

plt.show()
