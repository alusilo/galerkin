import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation

# computational domain
xmin,xmax = -0,2.
ymin,ymax = -0,2.
# wave velocities
vs = 1.
vp = 2.
# material parameters
r = 1.
m = 1.
l = 2.

res = 400
# wave parameters
n = np.array([1.,1.])
k = (2*np.pi/2)*n
wl = 2*np.pi/np.linalg.norm(k)
fmax = vp/wl
w = 2*np.pi*fmax
T = 2.0#1/fmax
# total frames
frames = 100
dt = T/frames

print("====== VELOCITIES ======")
print("vs = {}[m/s]".format(vs))
print("vp = {}[m/s]".format(vp))
print("====== LAME PARAMETERS ======")
print("rho = {}[m/s]".format(r))
print("mu = {}[m/s]".format(m))
print("lambda = {}[m/s]".format(l))
print("====== WAVE PARAMETERS ======")
print("k = {}[rad/m]".format(k))
print("w = {}[rad/s]".format(w))
print("T = {}[s]".format(T))

x = np.linspace(xmin,xmax,res)
y = np.linspace(ymin,ymax,res)

X, Y = np.meshgrid(x,y)

def V(t):
	#Uo = -(vs*n[1]*np.sin(k[0]*X + k[1]*Y) + vp*n[0]*np.sin(k[0]*X + k[1]*Y))
	return self.x*self.y*np.sin(w*t)

fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)

frame = V(0.0)
vmin,vmax = frame.min(),frame.max()
im = plt.imshow(frame, extent=[xmin,xmax,ymax,ymin], aspect='equal', cmap=plt.get_cmap('jet'), origin='upper', animated=True)
im.set_clim([-2, 2])
title = ax.set_title("Wave propagation {}[s]".format(0.0))
bar = plt.colorbar()

def updatefig(*args):
	global frame
	t = (args[0]+1)*dt
	frame = V(t)
	im.set_array(frame)
	ax.set_title('Wave propagation (t={0:.3f}[s])'.format(t))
	#im.set_clim([frame.min(), frame.max()])
	return im,

ani = animation.FuncAnimation(fig, updatefig, frames=frames, interval=50, repeat_delay=400, blit=False, repeat=True)

plt.show()