import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
import matplotlib.tri as tri

import os, json

from .Mesh2D import *

fig, ax = plt.subplots(1, 1, sharex=False, sharey=False, figsize=(9,5))

NDG_BINARY_FILE = 'movieVx.npy'
PARAM_FILENAME  = 'model.param'

PROJECT_NAME 		= 'erase'
DG_ROOT 			= os.path.dirname(os.getcwd())
LOCAL_PROJECT_DIR 	= os.path.join('resources/output', PROJECT_NAME)
PROJECT_DIR 		= os.path.join(DG_ROOT, LOCAL_PROJECT_DIR)

PARAM_FILE_PATH = os.path.join(PROJECT_DIR, PARAM_FILENAME)
MOVIE_FILE_PATH = os.path.join(PROJECT_DIR, NDG_BINARY_FILE)

param = json.load(open(PARAM_FILE_PATH))
duration = param['duration']
xmin,xmax,ymin,ymax = param['limits']
source = param['source']['position']
gather = param['gather']
print(gather)
mesh_file = param['mesh']
pml_layer = param['pml_layer']

pml_line1 = np.array([
	[xmin + pml_layer[0], xmin + pml_layer[0]],
	[ymin, ymax]
])

pml_line2 = np.array([
	[xmin, xmax],
	[ymax - pml_layer[3], ymax - pml_layer[3]]
])

pml_line3 = np.array([
	[xmax - pml_layer[1], xmax - pml_layer[1]],
	[ymax, ymin]
])

ax.plot(pml_line1[0]/1000, pml_line1[1]/1000, color='white', linewidth=0.5)
ax.plot(pml_line2[0]/1000, pml_line2[1]/1000, color='white', linewidth=0.5)
ax.plot(pml_line3[0]/1000, pml_line3[1], color='white', linewidth=0.5)

obj = MeshReader(
	mesh_file 	 = mesh_file,
	src_position = source,
	gather 		 = gather,
	src_smooth 	 = 0,
	plot = False
)

K, Nv, VX, VY, e2v = obj.get_values()

data = np.load(MOVIE_FILE_PATH)

frames = len(data)
print("total frames: {}".format(frames))
frame = data[0].T

vmin,vmax = data.min()/1.0,data.max()/1.0

dt = duration/frames

ax.set_xlim(xmin/1000,xmax/1000)
ax.set_ylim(ymax/1000,ymin/1000)
#ax.set_xlim(source[0]/1000-0.35,source[0]/1000+0.35)
#ax.set_ylim(source[1]/1000+0.35,source[1]/1000-0.35)
im = plt.imshow(frame, extent=[xmin/1000,xmax/1000,ymax/1000,ymin/1000], aspect='equal', cmap=plt.get_cmap('gray'), origin='upper', animated=True)
im.set_clim([vmax, vmin])
#title = ax.set_title("Wave propagation {}[s]".format(0.0))
ax.set_xlabel(r'$x$ [$km$]')
ax.set_ylabel(r'$z$ [$km$]')
bar = plt.colorbar()

bar.remove()
plt.draw()

for gth in gather:
 	gx,gy = gth
 	ax.scatter(gx/1000, gy/1000, color='white', marker='v', edgecolors='black', s=100)

sx,sy = source
ax.scatter(sx/1000, sy/1000, color='yellow', marker='*', edgecolors='black', s=100)

mesh = tri.Triangulation(VX/1000,VY/1000)
#plt.triplot(mesh, lw=0.1, color='black')

pause = False
def simData():
	global pause
	frame = data[0].T
	i = 0
	while i < frames:
		if not pause:
			frame = data[i].T
			#if i == 25:
			#	pause = True
			i = i+1
		yield frame, i

def onClick(event):
	global pause
	pause ^= True

def updatefig(simData):
	frame, i = simData
	im.set_array(frame)
	#ax.set_title('Wave propagation (t={0:.3f}[s])'.format((i)*dt))
	return im,

fig.canvas.mpl_connect('button_press_event', onClick)
ani = animation.FuncAnimation(fig, updatefig, simData, interval=50, repeat_delay=400, blit=False, repeat=True)
#ani.save('animation1.mp4')
plt.show()
