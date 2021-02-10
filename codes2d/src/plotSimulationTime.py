from __future__ import unicode_literals
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

#sigma=0m (1 element), t=0.5s, order=2:  3m26.229s

#sigma=10m (2 elements), t=0.5s, order=2, order_src=2: 3m25.963s
#sigma=10m (2 elements), t=0.5s, order=2, order_src=3: 4m41.244s
#sigma=10m (2 elements), t=0.5s, order=2, order_src=4: 5m35.427s
#sigma=10m (2 elements), t=0.5s, order=2, order_src=5: 7m55.755s

#sigma=10m (2 elements), t=0.5s, order=2, order_src=2: 3m25.963s
#sigma=20m (4 elements), t=0.5s, order=2, order_src=2: 3m38.397s
#sigma=30m (9 elements), t=0.5s, order=2, order_src=2: 3m26.763s
#sigma=40m (17 elements), t=0.5s, order=2, order_src=2: 3m27.900s
#sigma=50m (17 elements), t=0.5s, order=2, order_src=2: 5m16.075s

#sigma=0m (2 elements), t=0.5s, order=2, order_src=3: 5m13.439s
#sigma=10m (2 elements), t=0.5s, order=2, order_src=3: 4m53.768s
#sigma=20m (4 elements), t=0.5s, order=2, order_src=3: 5m16.340s
#sigma=30m (9 elements), t=0.5s, order=2, order_src=3: 5m19.279s
#sigma=40m (17 elements), t=0.5s, order=2, order_src=3: 5m15.212s
#sigma=50m (17 elements), t=0.5s, order=2, order_src=3: 5m17.358s

plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
plt.rc('legend', fontsize=6)
plt.rc('axes', titlesize=6)

yFormatter = FormatStrFormatter('%.1f')
xIntFormatter = FormatStrFormatter('%d')

width = 3.487
height = width/1.618

# data for fixed sigma = 10m
data = {
	'order': [2, 3, 4, 5],
	'time': [205.963, 281.244, 335.427, 475.755]
}
#####################################################################
fig, ax = plt.subplots()
fig.subplots_adjust(left=.18, bottom=.17, right=.98, top=.92)

ax.plot(data['order'], data['time'], '-', markersize=1, linestyle='--', c='black', linewidth=0.25)

ax.set_xlabel(r"orden de interpolaci\'on")
ax.set_ylabel(u'tiempo ($s$)')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_major_formatter(yFormatter)
ax.xaxis.set_major_formatter(xIntFormatter)

fig.set_size_inches(width, height)
fig.savefig('order_vs_time.pdf')

# data for fixed order = 2
data1 = {
	'sigma': [10, 20, 30, 40],
	'time': [205.963, 218.397, 206.763, 207.900]
}
data2 = {
	'sigma': [10, 20, 30, 40],
	'time': [293.768, 316.34, 319.279, 315.212]
}
#####################################################################
fig, ax = plt.subplots()
fig.subplots_adjust(left=.18, bottom=.17, right=.98, top=.92)

ax.plot(data1['sigma'], data1['time'], '-', markersize=1, linestyle='--', c='black', linewidth=0.25, label="2do orden")
ax.plot(data2['sigma'], data2['time'], '-', markersize=1, linestyle='-.', c='black', linewidth=0.25, label="3er orden")

ax.set_xlabel(u'radio del soporte ($m$)')
ax.set_ylabel(u'tiempo ($s$)')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_major_formatter(yFormatter)

plt.legend(loc=0)

fig.set_size_inches(width, height)
fig.savefig('sigma_vs_time.pdf')

src_freq = 7.5
src_delay = 1/src_freq

#####################################################################
def createWavelet(time):
		dispv 	= 0
		half 	= 0.5
		one 	= 1.0
		two 	= 2.0
		tau 	= time-src_delay
		twopi 	= two*np.pi
		a1 		= -2.0
		a2   	=-half*half*(twopi*src_freq)**2
		temp 	= a2*tau*tau

		if dispv == 0:
			'''
			-------------------------------------------
			--- Choix source, sortie en deplacement ---
			-------------------------------------------
			'''
			# *** Gaussienne, sortie en deplacement
			wout = np.exp(temp)
			# *** derivee premiere de Gaussienne, sortie en deplacement
			# wout = a1*( tau )*np.exp(temp)
			# *** Ricker (derivee seconde de gaussienne), sortie en deplacement
			# wout = a1*( half + temp )*np.exp(temp)
		else:
			a1 = two*a1
			'''
			---------------------------------------
			--- Choix source, sortie en vitesse ---
			---------------------------------------
			'''
			# *** Gaussienne, sortie en vitesse
			wout = a1*( tau )*np.exp(temp)
			# *** derivee premiere de Gaussienne, sortie en vitesse
			# wout = a1*( half + temp )*np.exp(temp)
			# *** Ricker (derivee seconde de gaussienne), sortie en vitesse
			# wout = a1*a2*tau*np.exp(temp)*(3.0 + 2.0*temp)
		#wout = 4.513*np.exp(-(self.src_freq**2)*(time-2*self.src_delay)**2)
		
		fig, ax = plt.subplots()
		fig.subplots_adjust(left=.12, bottom=.17, right=.98, top=.92)
		#plt.grid()
		plt.xlabel(r"tiempo ($s$)")
		plt.ylabel("Amplitud")
		plt.plot(time,wout, color='black', linewidth=1)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)

		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_major_formatter(yFormatter)

		fig.set_size_inches(width, height)
		fig.savefig('wavelet.pdf')

time = np.linspace(0, 2.0, 1000)
createWavelet(time)