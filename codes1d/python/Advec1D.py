from numpy import zeros, min, pi, ceil
from AdvecRHS1D import *

def Advec1D(Np, K, x, rx, nx, Dr, u, LIFT, Fscale, Nfp, Nfaces, vmapM, vmapP, vmapI, mapI, mapO, FinalTime):
	time = 0.
	rk4a = [0.0, -567301805773.0/1357537059087.0, -2404267990393.0/2016746695238.0, -3550918686646.0/2091501179385.0, -1275806237668.0/842570457699.0]
	rk4b = [1432997174477.0/9575080441755.0, 5161836677717.0/13612068292357.0, 1720146321549.0/2090206949498.0, 3134564353537.0/4481467310338.0, 2277821191437.0/14882151754819.0]
	rk4c = [0.0, 1432997174477.0/9575080441755.0, 2526269341429.0/6820363962896.0, 2006345519317.0/3224310063776.0, 2802321613138.0/2924317926251.0]
	resu = zeros((Np,K))
	xmin = min(abs(x[0,:]-x[1,:]))
	CFL=0.75
	dt   = CFL/(2*pi)*xmin
	dt = .5*dt
	Nsteps = int(ceil(FinalTime/dt))
	dt = float(FinalTime)/Nsteps

	a = 2*pi;
	data = []
	data.append(u)
	for tstep in range(0,Nsteps):
		for INTRK in range(0,5):
			timelocal = time + rk4c[INTRK]*dt
			rhsu = AdvecRHS1D(u, timelocal, a, Nfp, Nfaces, K, LIFT, Fscale, vmapM, vmapP, vmapI, mapI, mapO, nx, rx, Dr)
			resu = rk4a[INTRK]*resu + dt*rhsu
			u = u+rk4b[INTRK]*resu
		if tstep%60 == 0:
			data.append(u)
		time = time+dt

	return data
