import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
import matplotlib.tri as tri

data = np.load('../resources/output/tracesUz.npy')
nf = len(data)-1
fig, ax = plt.subplots(nf, 1, sharex=False, sharey=False)

Nsteps = 1300
dt = 0.001232

time = data[0]
traces = data[1:]
ntraces = len(data[1:])
with open('/media/jdsg/CODE/DGCrack2D/INPLANE/demos/Point_source_h300m/VY', 'rb') as f:
	data2 = np.fromfile(f, count=ntraces*Nsteps, dtype=np.float32).reshape(-1,ntraces).T

dgtime = np.linspace(0,time[-1],Nsteps)

if nf == 1:
	dt = np.zeros(traces[0].shape)
	dt[0:-1] = np.diff(traces[0])/np.diff(time)
	dt[-1] = (traces[0][-1] - traces[0][-2])/(time[-1] - time[-2])
	ax.plot(time,traces[0])
	ax.grid()
	ax.set_title('Seismic Traces')
	ax.set_xlabel(r'$t[s]$')
	ax.set_ylabel('Amplitude')
else:
	for i in range(nf):
		ax[i].plot(time,traces[i],c='b')
		ax[i].plot(dgtime,data2[i,:Nsteps],c='g')
		ax[i].grid()
		ax[i].set_title('Seismic Traces')
		ax[i].set_xlabel(r'$t[s]$')
		ax[i].set_ylabel('Amplitude')
		ax[i].legend(['nodalDG', 'DGCrack'], loc=0)

plt.show()
