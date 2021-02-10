import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
import matplotlib.tri as tri

import os, json
from scipy import interpolate

width = 3.487
height = width/1.618

plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
plt.rc('legend', fontsize=6)
plt.rc('axes', titlesize=6)

def slope(data,dt):
	derivative = np.zeros(data.shape)
	for i in range(1,len(data)):
		dx = (data[i]-data[i-1])
		derivative[i] = dx/dt
	
	return derivative

PROJECT_NAME 		= 'gauss75_lamb_s0'
NDG_BINARY_FILE 	= 'tracesVy.npy'
PARAM_FILENAME  	= 'model.param'
ANALYTIC_DATA_PATH  = '../../EX2DDIR_Berg/Uz_file_ascii'
DG_ROOT 			= os.path.dirname(os.getcwd())
LOCAL_PROJECT_DIR 	= os.path.join('resources/output', PROJECT_NAME)
PROJECT_DIR 		= os.path.join(DG_ROOT, LOCAL_PROJECT_DIR)
PARAM_FILE_PATH 	= os.path.join(PROJECT_DIR, PARAM_FILENAME)

FILE_PATH = os.path.join(PROJECT_DIR, NDG_BINARY_FILE)

param = json.load(open(PARAM_FILE_PATH))
data = np.load(FILE_PATH)
nf = len(data)-1

Nsteps = param['time_steps']
dt = param['dt']
time = data[0]
traces = data[1:]
ntraces = len(data[1:])
with open(ANALYTIC_DATA_PATH, 'r') as f:
	data2 = []
	for line in f.readlines():
		data2.append(float(line))
	data2 = np.array(data2).reshape(ntraces,-1)

fig, ax = plt.subplots(nf, 1, sharex=True, sharey=False, figsize=(10,4))
fig.subplots_adjust(left=.12, bottom=.17, right=.98, top=.92)

if nf == 1:
	dev = np.gradient(data2[0],time)
	ynew = -slope(data2[0],dt)
	xnew = time
	f = interpolate.interp1d(time, dev)
	xnew = np.arange(0, time[-1], 0.008)
	ynew = f(xnew)
	ymin, ymax = np.min(traces[0]), np.max(traces[0]) 
	ax.plot(time,traces[0],'r--',c='black', linewidth=1, marker='s', markevery=50,markerfacecolor='none',  markersize=3, label='nodalDG')
	ax.plot(xnew,ynew, c='black', linewidth=0.5, label=r"Anal\'itica")
	#ax.grid()
	#ax.set_title('Registro Sísmico')
	ax.set_xlabel(r'tiempo ($s$)')
	ax.set_ylabel('Amplitud')
	#ax.set_ylim(ymin-0.2*np.abs(ymax-ymin),ymax+0.2*np.abs(ymax-ymin))
	ax.legend(loc=0)
else:
	for i in range(nf):
		t,dev = slope(data2[i],time)
		ymin, ymax = np.min(traces[i]), np.max(traces[i]) 
		#ax[0].set_title('Registro Sísmico')
		ax[i].plot(t,-dev, c='black', linewidth=0.5, label=r"Anal\'itica")
		ax[i].plot(time,traces[i],'r--',c='black',linewidth=0.5,marker='s',markevery=50,markerfacecolor='none', markersize=1, label='nodalDG')
		#ax[i].grid()
		ax[i].set_ylabel('Amplitud')
		#ax[i].set_ylim(ymin-0.2*np.abs(ymax-ymin),ymax+0.2*np.abs(ymax-ymin))
		ax[i].legend(loc=0)
	ax[nf-1].set_xlabel(r'$t$ ($s$)')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

fig.set_size_inches(width, height)

plt.savefig('seismogram_uz.pdf')
#plt.show()
