from numpy import *
from numpy import transpose as Trans

def getm(v,w):
	return array([v[i] for i in w])

def AdvecRHS1D(u,time,a,Nfp,Nfaces,K,LIFT, Fscale, vmapM,vmapP, vmapI, mapI, mapO, nx, rx, Dr):
	alpha=0
	du = zeros(Nfp*Nfaces*K)
	du = (Trans(u).ravel()[vmapM]-Trans(u).ravel()[vmapP])*(a*Trans(nx).ravel()-(1-alpha)*abs(a*Trans(nx).ravel()))/2
	uin = -sin(a*time)
	du[mapI] = (Trans(u).ravel()[vmapI] - uin)*(a*Trans(nx).ravel()[mapI]-(1-alpha)*abs(a*Trans(nx).ravel()[mapI]))/2
	du[mapO] = 0
	du = reshape(du, (K,Nfp*Nfaces)).transpose()
	rhsu = -a*rx*dot(Dr,u) + dot(LIFT,(Fscale*(du)))
	
	return rhsu
