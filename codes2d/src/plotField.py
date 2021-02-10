import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
import matplotlib.tri as tri

import os, json

def load_vertex(path):
	if os.path.exists(path):
		with open(path, 'rt') as f:
			data = f.read().splitlines()
			# line identification for nodes
			nodeId = data.index('$Nodes')+1
			# line identification for elements
			elementId = data.index('$Elements')+1
			# Number of vertex
			Nv = int(data[nodeId])
			# Nodes array
			nodes 	 = np.array([list(map(float, a.split(' ')[1:3])) for a in data[nodeId+1:nodeId+Nv+1]])
			VX = nodes[:,0]
			VY = nodes[:,1]
			
		print("Mesh file {} loadded successfully.".format(path.split('/')[-1]))

		return VX, VY
	else:
		print("Error: Especified file '{}' does not exists.".format(path))
		exit(1)

def load_vertex2(path):
	if os.path.exists(path):
		with open(path, 'rt') as f:
			data = f.readlines()
			Nv,K  = np.array(data[6].split()[:2]).astype(int)
			VX,VY = np.array(list(map(lambda a: a.split()[1:], data[9:Nv+9])), dtype=float).T

		print("Mesh file {} loadded successfully.".format(path.split('/')[-1]))

		return VX, VY
	else:
		print("Error: Especified file '{}' does not exists.".format(path))
		exit(1)

fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)

param = json.load(open('../resources/output/model.param'))
duration = param['duration']
xmin,xmax,ymin,ymax = param['limits']
source = param['source']
gather = param['gather']
mesh_file = param['mesh']

data = np.load('../resources/output/movieUx.npy')
VX, VY = load_vertex(mesh_file)

frames = len(data)
print("total frames: {}".format(frames))
frame = data[0].T

vmin,vmax = data.min()/15,data.max()/15

dt = duration/frames

im = plt.imshow(frame, extent=[xmin,xmax,ymax,ymin], aspect='equal', cmap=plt.get_cmap('jet'), origin='upper', animated=True)
im.set_clim([vmax, vmin])
title = ax.set_title("Wave propagation {}[s]".format(0.0))
bar = plt.colorbar()

for gth in gather:
 	gx,gy = gth
 	ax.scatter(gx, gy, color='blue', marker='v', s=50)

sx,sy = source
ax.scatter(sx, sy, color='red', marker='*', s=50)

mesh = tri.Triangulation(VX,VY)
plt.triplot(mesh, lw=0.1, color='black')
def updatefig(*args):
	global frame
	i = args[0]
	frame = data[i].T
	im.set_array(frame)
	ax.set_title('Wave propagation (t={0:.3f}[s])'.format((i+1)*dt))
	return im,

ani = animation.FuncAnimation(fig, updatefig, frames=frames, interval=50, repeat_delay=400, blit=False, repeat=True)
plt.show()
