import numpy as np

# m/s = kg.m^-2.s^-2.(m^3.kg^-1.s)

def source(t):
	fc = 7.5 # Hz
	a1 = -2000.0 # kg.m^-2.s^-2
	a2 = -(np.pi*fc)**2 # Hz^2
	td = 1/fc # sec
	rho = 2200. # Kg.m^-3
	return a1*(0.5 + a2*(t-td)**2)*np.exp(a2*(t-td)**2)/rho

def Dsource(t):
	fc = 7.5 # Hz
	a1 = -2000.0 # kg.m^-2.s^-2
	a2 = -(np.pi*fc)**2 # Hz^2
	td = 1/fc # sec
	rho = 2200. # Kg.m^-3
	return 2*a1*a2*(t-td)*(1.5 + a2*(t-td)**2)*np.exp(a2*(t-td)**2)/rho

def Gsource(t):
	fc = 7.5 # Hz
	a1 = -2000.0 # kg.m^-2.s^-2
	a2 = -(np.pi*fc)**2 # Hz^2
	td = 1/fc # sec
	rho = 2200. # Kg.m^-3
	return 4.513*np.exp(-(fc**2)*(t-2*td)**2)

dt = 0.001
time = np.arange(0.0,1.0,dt)
with open('../gauss_7_5a.src', 'w') as f:
	for t in time:
		f.write(str(Gsource(t)))
		f.write('\n')
