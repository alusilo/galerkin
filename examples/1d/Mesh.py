from numpy import *
from numpy.linalg import *
from numpy import transpose as Trans

# Returns the index of the matrix in a row
def sub2ind(array_shape, rows, cols):
	ind = rows*array_shape[1] + cols
	ind[ind < 0] = -1
	ind[ind >= array_shape[0]*array_shape[1]] = -1
	return ind
# Returns a sub array from v taking the index values in w
def getm(v,w):
	return array([v[i] for i in w])

# Return the index where the arrays a and b are equal
def myfind(a,b):
	return array([i for i in range(0,len(a)) if a[i] == b[i]])

# Generate a 1D mesh
def MeshGen1D(xmin,xmax,K):
	dx = (xmax-xmin)/K
	return array([xmin+i*dx for i in range(0,K+1)]), array([[i,i+1] for i in range(0,K)])

# Compute the metric elements for the local mappings of the 1D elements
def GeometricFactors1D(x,Dr):
	xr = dot(Dr,x)
	J = xr
	rx = 1/J
	return rx, J

# Build global connectivity arrays for 1D grid based on standard EToV input array from grid generator
def Connect1D(EToV, Nfaces, K, Nv):
	TotalFaces = Nfaces*K
	vn = [0,1]
	SpFToV = zeros((TotalFaces,Nv))#spalloc(TotalFaces, Nv, 2*TotalFaces);
	sk = 0
	for k in range(0,K):
		for face in range(0,Nfaces):
			SpFToV[sk, EToV[k, vn[face]]] = 1;
			sk = sk+1;

	SpFToF = dot(SpFToV,Trans(SpFToV)) - identity(TotalFaces);
	faces2, faces1 = nonzero(SpFToF==1.0)
	element1 = array(floor((faces1-2)/Nfaces),dtype=int) + 1
	face1    = array(mod((faces1-2), Nfaces),dtype=int)
	element2 = array(floor( (faces2-2)/Nfaces ), dtype=int)  + 1
	face2    = array(mod( (faces2-2), Nfaces ), dtype=int)
	i=element1.ravel()
	j=face1.ravel()
	EToE = dot(Trans(array([range(0,K)])),ones((1,Nfaces), dtype=int))
	EToF = dot(ones((K,1), dtype=int),array([range(0,Nfaces)]))
	EToE[i,j] = element2
	EToF[i,j] = face2
	return EToE, EToF

def BuildMaps1D(K, Np, Nfp, Nfaces, EToE, EToF, Fmask, x, NODETOL):
	nodeids = array(reshape(range(0,K*Np), (K, Np)).transpose())
	vmapM   = zeros((Nfp, Nfaces, K), dtype=int)
	vmapP   = zeros((Nfp, Nfaces, K), dtype=int)
	for k1 in range(0,K):
		for f1 in range(0,Nfaces):
			vmapM[:,f1,k1] = nodeids[Fmask[f1], k1]

	for k1 in range(0,K):
		for f1 in range(0,Nfaces):
			k2 = EToE[k1,f1]
			f2 = EToF[k1,f1]
			vidM = vmapM[:,f1,k1]
			vidP = vmapM[:,f2,k2]

			x1  = Trans(x).ravel()[vidM[0]]
			x2  = Trans(x).ravel()[vidP[0]]
			D = (x1 - x2)**2
			if D < NODETOL:
				vmapP[:,f1,k1] = vidP

	vmapP = Trans(vmapP[0]).ravel()
	vmapM = Trans(vmapM[0]).ravel()
	mapB = myfind(vmapP,vmapM)
	vmapB = getm(vmapM,mapB)

	mapI = 0
	mapO = K*Nfaces-1
	vmapI = 0
	vmapO = K*Np-1

	return vmapM, vmapP, vmapI, mapI, vmapB, mapB, mapO