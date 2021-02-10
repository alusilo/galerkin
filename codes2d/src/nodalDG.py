from Mesh2D import *

#import mkl

#mkl.set_num_threads(8)
#print("max number of threads: {}".format(mkl.get_max_threads()))

#np.set_printoptions(precision=4,suppress=True)

# Low storage Runge-Kutta coefficients
rk4a = [
	 0.0,
    -567301805773.0/1357537059087.0,
    -2404267990393.0/2016746695238.0,
    -3550918686646.0/2091501179385.0,
    -1275806237668.0/842570457699.0
]
rk4b = [
	1432997174477.0/9575080441755.0,
	5161836677717.0/13612068292357.0,
	1720146321549.0/2090206949498.0,
	3134564353537.0/4481467310338.0,
	2277821191437.0/14882151754819.0
]
rk4c = [
	0.0,
	1432997174477.0/9575080441755.0,
	2526269341429.0/6820363962896.0,
	2006345519317.0/3224310063776.0,
	2802321613138.0/2924317926251.0
]

Wall 		= 1
Source 		= 2
Free		= 3
Gather 		= 4

"""
	START MATH PART
"""
def Nodes2D(N):
	Np = int((N+1)*(N+2)/2)
	alpopt = [0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999, 1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258]
	alpha  = alpopt[N-1] if N < 16 else 5/3

	L1 = np.zeros(Np)
	L2 = np.zeros(Np)
	L3 = np.zeros(Np)
	i=0
	for n in range(N+1):
		for m in range(N-n+1):
			if N == 0: L1[i] = 0
			else: L1[i] = n/N
			if N == 0: L3[i] = 0
			else: L3[i] = m/N
			i=i+1

	L2 = 1.0 - L1 - L3
	x = -L2+L3
	y = (-L2-L3+2*L1)/np.sqrt(3)

	blend1 = 4*L2*L3
	blend2 = 4*L1*L3
	blend3 = 4*L1*L2

	warpf1 = Warpfactor(N,L3-L2)
	warpf2 = Warpfactor(N,L1-L3)
	warpf3 = Warpfactor(N,L2-L1)

	warp1 = blend1*warpf1*(1+(alpha*L1)**2)
	warp2 = blend2*warpf2*(1+(alpha*L2)**2)
	warp3 = blend3*warpf3*(1+(alpha*L3)**2)

	x = x + 1*warp1 + np.cos(2*np.pi/3)*warp2 + np.cos(4*np.pi/3)*warp3
	y = y + 0*warp1 + np.sin(2*np.pi/3)*warp2 + np.sin(4*np.pi/3)*warp3

	return x,y

def xytors(x,y):
	L1 = (np.sqrt(3)*y + 1)/3
	L2 = (-3*x - np.sqrt(3)*y + 2)/6
	L3 = ( 3*x - np.sqrt(3)*y + 2)/6

	r = -L2 + L3 - L1
	s = -L2 - L3 + L1

	return r,s

def rstoab(r,s):
	Np = len(r)
	a = np.zeros(Np)
	for i in range(Np):
		if s[i] != 1:
			a[i] = 2*(1 + r[i])/(1 - s[i]) - 1
		else:
			a[i] = -1

	b = s

	return a,b

def JGLPoints(alpha,beta,N):
	x = [] if (N==1) else sf.j_roots(N-1,alpha,beta)[0]
	x = np.append(np.append([-1.],x),[1.])
	return x

def JPNormalized(N,alpha,beta,r):
	gamma0 = (2.**(alpha+beta+1)/(2*N+alpha+beta+1))*(sf.gamma(N+alpha+1)*sf.gamma(N+beta+1)/(sf.gamma(N+alpha+beta+1)*sf.factorial(N)))
	norm = 1./np.sqrt(gamma0)
	x = norm*sf.eval_jacobi(N,alpha,beta,r)
	return x

def GradJacobiP(r,alpha,beta,N):
	dP = np.zeros(len(r))
	if N==0:
		return dP
	else:
		return np.sqrt(N*(N+alpha+beta+1))*JPNormalized(N-1,alpha+1,beta+1,r)

def Vandermonde1D(N,r):
	V = np.zeros((len(r),N+1))
	for i in range(N+1):
		V[:,i] = JPNormalized(i,0,0,r)
	return V

def Vandermonde2D(N,r,s):
	V2D = np.zeros((len(r),int((N+1)*(N+2)/2)))
	a,b = rstoab(r,s)
	k=0
	for i in range(N+1):
		for j in range(N-i+1):
			V2D[:,k] = Simplex2DP(a,b,i,j)
			k = k+1
	return V2D

def GradVandermonde2D(N,r,s):
	V2Dr = np.zeros((len(r),int((N+1)*(N+2)/2)))
	V2Ds = np.zeros((len(r),int((N+1)*(N+2)/2)))
	a,b = rstoab(r,s)
	k=0
	for i in range(N+1):
		for j in range(N-i+1):
			V2Dr[:,k],V2Ds[:,k] = GradSimplex2DP(a,b,i,j)
			k = k+1
	return V2Dr,V2Ds

def Simplex2DP(a,b,i,j):
	h1 = JPNormalized(i,0,    0,a)
	h2 = JPNormalized(j,2*i+1,0,b)
	P = np.sqrt(2.0)*h1*h2*(1-b)**i
	return P

def GradSimplex2DP(a,b,idx,jdx):
	fa 	= JPNormalized(idx,0,0,a)
	dfa = GradJacobiP(a,0,0,idx)
	gb 	= JPNormalized(jdx,2*idx+1,0,b)
	dgb = GradJacobiP(b,2*idx+1,0,jdx)
	dmodedr = dfa*gb
	dmodeds = dfa*(gb*(0.5*(1+a)))
	tmp = dgb*((0.5*(1-b))**idx)
	if idx>0:
		dmodedr = dmodedr*((0.5*(1-b))**(idx-1))
		dmodeds = dmodeds*((0.5*(1-b))**(idx-1))
		tmp = tmp-0.5*idx*gb*((0.5*(1-b))**(idx-1))

	dmodedr = 2**(idx+0.5)*dmodedr
	dmodeds = 2**(idx+0.5)*(dmodeds+fa*tmp)

	return dmodedr, dmodeds

def Warpfactor(N,r):
	LGLr = JGLPoints(1,1,N)
	req = np.linspace(-1,1,N+1)
	Veq = Vandermonde1D(N,req)
	Pmat = np.zeros((N+1,len(r)))
	for i in range(N+1):
		Pmat[i,:] = JPNormalized(i,0,0,r)

	Lmat = la.solve(Veq.T,Pmat)
	warp = np.dot(Lmat.T,(LGLr - req))
	zerof = 1*(abs(r) < 1.0 - 1.0e-10)
	sf = 1.0 - (zerof*r)**2
	warp = warp*(1/sf) + warp*(zerof-1)

	return warp

"""
	END MATH PART
"""

class NodalDG2D(MeshReader):
	"""docstring for NodalDG2D"""
	def __init__(self, **kwargs):
		super(NodalDG2D, self).__init__(**kwargs)
		self.order = kwargs['order']
		self.source_order = kwargs['source_order']

		self.N = 0
		self.Nfp = self.N + 1
		self.Np = int(self.Nfp*(self.Nfp+1)/2)

		self.Nfaces = 3
		self.NODETOL = 1e-12

		self.pml_layer = kwargs['pml_layer']

		self.vmapG = [[] for __ in range(len(self.gather))]
		self.vmapGO = np.zeros(len(self.gather), dtype=int)

		self.e2e, self.e2f = self.tiConnect2D()

		self.bc_wall = np.zeros((3,self.K))
		self.bc_free = np.zeros((3,self.K))
		self.SrcDom = np.zeros(self.K,dtype=int)

		e4f = np.dot(np.ones((self.Nfaces,1), dtype=int),[np.arange(self.K)])
		walls = np.nonzero(self.e2e == e4f)
		self.bc_wall[walls] = Wall

		self.map_walls = walls[1]

		vtx = self.VY[self.e2v[:,walls[1]]] - self.ymin <= self.NODETOL
		ids = np.nonzero(np.sum(vtx,axis=0) == 2)[0]

		free = (walls[0][ids],walls[1][ids])

		self.bc_free[free] = Free

		self.map_free = free[1]

		self.SrcDom[self.surrounding] = Source
		self.GthDom = np.zeros((len(self.gather),self.K))

		for i,k in enumerate(self.mapg):
			self.GthDom[i,k] = Gather

		self.pml_elements = self.GetPMLElements()

		self.Norder = self.order*np.ones(self.K,dtype=int)
		self.Norder[self.surrounding] = self.source_order

		self.Nmin = np.min(self.Norder)
		self.Nmax = np.max(self.Norder)

		self.param_file = kwargs['param_file']

		# read vs, vp and rho from file
		data = []
		if os.path.exists(self.param_file):
			with open(self.param_file, 'rt') as f:
				for line in f.readlines():
					if not line.startswith('#'):
						data.append(line)
				self.rho, self.vs, self.vp = np.array(list(map(lambda a: a.split()[:], data[:])), dtype=float).T
			print("Info: Physical parameters file {} loadded successfully.".format(self.param_file))
		else:
			print("Error: Especified file '{}' does not exists.".format(self.param_file))
			exit(-1)

		# lame parameters
		self.xmu = self.rho*self.vs**2
		self.xlambda = self.rho*(self.vp**2 - 2*self.vs**2)

		# Set up arbitrary order elements mesh
		self.pinfo = self.BuildPNonCon2D()

		self.vmapPML = []
		for obj in self.pinfo:
			if obj.K > 0:
				self.vmapPML += obj.vmapPML

		self.vmapPML = list(set(np.array(self.vmapPML)))

	def GetPMLElements(self):
		print("Info: Identifying elements belonging ABC zone...")
		cx = np.sum(self.VX[self.e2v[:,range(self.K)]],axis=0)/3
		cy = np.sum(self.VY[self.e2v[:,range(self.K)]],axis=0)/3
		ids1 = np.nonzero(cx < self.xmin + self.pml_layer[0])[0]
		ids2 = np.nonzero(cx > self.xmax - self.pml_layer[1])[0]
		ids3 = np.nonzero(cy < self.ymin + self.pml_layer[2])[0]
		ids4 = np.nonzero(cy > self.ymax - self.pml_layer[3])[0]

		return np.concatenate((ids1,ids2,ids3,ids4))

	def Dmatrices2D(self):
		print("Info: Calculating derivative matrices...")
		Vr,Vs = GradVandermonde2D(self.N,self.r,self.s)
		Dr = np.dot(Vr,la.inv(self.V))
		Ds = np.dot(Vs,la.inv(self.V))
		return Dr,Ds

	def Nodes2D(self):
		print("Info: Calculating elements nodes...")
		alpopt = [0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999, 1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258]
		alpha  = alpopt[self.N-1] if self.N < 16 else 5/3

		L1 = np.zeros(self.Np)
		L2 = np.zeros(self.Np)
		L3 = np.zeros(self.Np)
		i=0
		for n in range(self.N+1):
			for m in range(self.N-n+1):
				L1[i] = n/self.N
				L3[i] = m/self.N
				i=i+1

		L2 = 1.0 - L1 - L3
		x = -L2+L3
		y = (-L2-L3+2*L1)/np.sqrt(3)

		blend1 = 4*L2*L3
		blend2 = 4*L1*L3
		blend3 = 4*L1*L2

		warpf1 = Warpfactor(self.N,L3-L2)
		warpf2 = Warpfactor(self.N,L1-L3)
		warpf3 = Warpfactor(self.N,L2-L1)

		warp1 = blend1*warpf1*(1+(alpha*L1)**2)
		warp2 = blend2*warpf2*(1+(alpha*L2)**2)
		warp3 = blend3*warpf3*(1+(alpha*L3)**2)

		x = x + 1*warp1 + np.cos(2*np.pi/3)*warp2 + np.cos(4*np.pi/3)*warp3
		y = y + 0*warp1 + np.sin(2*np.pi/3)*warp2 + np.sin(4*np.pi/3)*warp3

		return x,y

	def Lift2D(self):
		print("Info: Calculating fluxes matrices...")
		Emat = np.zeros((self.Nfaces*self.Nfp,self.Np))
		# face 1
		faceR = self.r[self.Fmask[0]]
		V1D = Vandermonde1D(self.N,faceR)
		massEdge1 = la.inv(np.dot(V1D,V1D.T))
		Emat[:self.Nfp,self.Fmask[0]] = massEdge1

		# face 2
		faceR = self.r[self.Fmask[1]]
		V1D = Vandermonde1D(self.N,faceR)
		massEdge2 = la.inv(np.dot(V1D,V1D.T))
		Emat[self.Nfp:2*self.Nfp,self.Fmask[1]] = massEdge2

		# face 3
		faceS = self.s[self.Fmask[2]]
		V1D = Vandermonde1D(self.N,faceS)
		massEdge3 = la.inv(np.dot(V1D,V1D.T))
		Emat[2*self.Nfp:3*self.Nfp,self.Fmask[2]] = massEdge3

		LIFT = np.dot(np.dot(Emat,self.V),self.V.T)

		return LIFT

	def GeometricFactors2D(self):
		xr = np.dot(self.x,self.Dr.T)
		yr = np.dot(self.y,self.Dr.T)
		xs = np.dot(self.x,self.Ds.T)
		ys = np.dot(self.y,self.Ds.T)
		J  = -xs*yr + xr*ys
		rx =  ys/J
		sx = -yr/J
		ry = -xs/J
		sy =  xr/J
		
		return rx,sx,ry,sy,J

	def Normals2D(self):
		xr = np.dot(self.x,self.Dr.T)
		yr = np.dot(self.y,self.Dr.T)
		xs = np.dot(self.x,self.Ds.T)
		ys = np.dot(self.y,self.Ds.T)
		J = - xs*yr + xr*ys

		fxr = xr[:,self.Fmask.flat]
		fxs = xs[:,self.Fmask.flat]
		fyr = yr[:,self.Fmask.flat]
		fys = ys[:,self.Fmask.flat]

		nx = np.zeros((np.size(self.x,0),3*self.Nfp))
		ny = np.zeros((np.size(self.x,0),3*self.Nfp))

		fid1 = np.arange(self.Nfp)
		fid2 = np.arange(self.Nfp,2*self.Nfp)
		fid3 = np.arange(2*self.Nfp,3*self.Nfp)

		# face 1
		nx[:,fid1] =  fyr[:,fid1]
		ny[:,fid1] = -fxr[:,fid1]

		# face 2
		nx[:,fid2] =  fys[:,fid2] - fyr[:,fid2]
		ny[:,fid2] = -fxs[:,fid2] + fxr[:,fid2]

		# face 3
		nx[:,fid3] = -fys[:,fid3]
		ny[:,fid3] =  fxs[:,fid3]

		sJ = np.sqrt(nx**2 + ny**2)
		nx = nx/sJ
		ny = ny/sJ

		return nx,ny,sJ

	def tiConnect2D(self):
		print("Info: Calculating elements connections...")
		K = np.size(self.e2v,1)
		Nnodes = np.max(self.e2v)+1
		# create list of all faces 1, then 2, & 3
		fnodes = np.concatenate((self.e2v[[0,1],:],self.e2v[[1,2],:],self.e2v[[2,0],:]),axis=1)
		fnodes = np.sort(fnodes,axis=0)
		# set up default element to element and Element to faces connectivity
		e2e = np.dot(np.ones((self.Nfaces,1),dtype=int),[np.arange(0,K)])
		e2f = np.dot(np.arange(0,self.Nfaces).reshape(-1,1),np.ones((1,K),dtype=int))
		index = (fnodes[0,:]*Nnodes + fnodes[1,:])
		spNodeToNode = np.concatenate(([index], [np.arange(0,self.Nfaces*K)], [e2e.flatten()], [e2f.flatten()]),axis=0)
		sortedNodes = np.array(sorted(spNodeToNode.T, key=operator.itemgetter(0, 1))).T
		# find matches in the sorted face list
		indices = np.nonzero(sortedNodes[0,0:-1]==sortedNodes[0,1:])[0]
		# make links reflexive
		matchL = np.concatenate((sortedNodes[:,indices],sortedNodes[:,indices+1]),axis=1)
		matchR = np.concatenate((sortedNodes[:,indices+1],sortedNodes[:,indices]),axis=1)
		# insert matches
		e2e.flat[matchL[1,:]] = matchR[2,:]
		e2f.flat[matchL[1,:]] = matchR[3,:]

		return e2e,e2f

	def BuildMaps2D(self):
		print("Info: Mapping elements nodes...")
		nodeids = np.arange(self.K*self.Np).reshape(self.K,self.Np)
		vmapM = np.zeros((self.K,self.Nfaces,self.Nfp), dtype=int)
		vmapP = np.zeros((self.K,self.Nfaces,self.Nfp), dtype=int)
		mapM  = np.arange(self.K*self.Nfp*self.Nfaces)
		mapP  = np.arange(self.K*self.Nfp*self.Nfaces).reshape((self.K,self.Nfaces,self.Nfp))
		# find index of face nodes with respect to volume node ordering
		for k1 in range(self.K):
			for f1 in range(self.Nfaces):
				vmapM[k1,f1,:] = nodeids[k1,self.Fmask[f1]]

		one = np.ones((1,self.Nfp), dtype=int)
		for k1 in range(self.K):
			for f1 in range(self.Nfaces):
				# find neighbor
				k2 = self.e2e[f1,k1]
				f2 = self.e2f[f1,k1]
				# reference length of edge
				v1 = self.e2v[f1,k1]
				v2 = self.e2v[np.mod(f1+1,self.Nfaces),k1]
				refd = np.sqrt((self.VX[v1]-self.VX[v2])**2 + (self.VY[v1]-self.VY[v2])**2)
				# find find volume node numbers of left and right nodes
				vidM = vmapM[k1,f1,:]
				vidP = vmapM[k2,f2,:]
				x1 = self.x.flat[vidM].reshape(-1,1)
				y1 = self.y.flat[vidM].reshape(-1,1)
				x2 = self.x.flat[vidP].reshape(-1,1)
				y2 = self.y.flat[vidP].reshape(-1,1)

				x1 = np.dot(x1,one)
				y1 = np.dot(y1,one)
				x2 = np.dot(x2,one)
				y2 = np.dot(y2,one)

				# Compute distance matrix
				D = (x1 - x2.T)**2 + (y1-y2.T)**2;
				idM, idP = np.nonzero(np.sqrt(abs(D))<self.NODETOL*refd)
				vmapP[k1,f1,idM] = vidP[idP]
				mapP[k1,f1,idM] = idP + f2*self.Nfp + k2*self.Nfaces*self.Nfp

		vmapP = vmapP.flatten()
		vmapM = vmapM.flatten()
		mapP = mapP.flatten()
		mapB = np.nonzero(vmapP==vmapM)[0]
		vmapB = vmapM[mapB]

		return mapM, mapP, vmapM, vmapP, vmapB, mapB

	def InterpMatrix2D(N, invV, rout, sout):
		# compute Vandermonde at (rout,sout)
		Vout = Vandermonde2D(N, rout, sout)
		# build interpolation matrix
		IM = np.dot(Vout,invV)
		return IM

	def BuildWallBCMaps2D(self):
		bct = self.bc_wall.T
		bnodes = np.dot(np.ones((self.Nfp,1),dtype=int),[bct.flatten()]).T
		bnodes = bnodes.flatten()
		# find location of boundary nodes in face and volume node lists
		mapW 	= np.nonzero(bnodes==Wall)[0]
		vmapW 	= self.vmapM[mapW]

		return mapW,vmapW

	def BuildFreeBCMaps2D(self):
		bct = self.bc_free.T
		bnodes = np.dot(np.ones((self.Nfp,1),dtype=int),[bct.flatten()]).T
		bnodes = bnodes.flatten()
		# find location of boundary nodes in face and volume node lists
		mapF 	= np.nonzero(bnodes==Free)[0]
		vmapF 	= self.vmapM[mapF]

		return mapF,vmapF

	def BuildBCSource(self,ids):
		bnodes = np.dot(np.ones((self.Np,1),dtype=int),[self.SrcDom]).T
		anodes = np.dot(np.ones((self.Np,1),dtype=int),[self.area_glob]).T
		mapS = np.nonzero(bnodes==Source)
		vmapS = ids[mapS]
		srcArea = anodes[mapS]

		return vmapS, srcArea

	def BuildGthDom(self,ids):
		for i,row in enumerate(self.GthDom):
			bnodes = np.dot(np.ones((self.Np,1),dtype=int),[row]).T
			mapG = np.nonzero(bnodes==Gather)
			if len(ids[mapG]) > 0:
				self.vmapG[i]  = ids[mapG]
				self.vmapGO[i] = self.N-1

	def BuildBCPML(self,ids):
		mapPL = np.nonzero((self.x - self.xmin) < self.pml_layer[0])
		mapPR = np.nonzero((self.xmax - self.x) < self.pml_layer[1])
		mapPD = np.nonzero((self.y - self.ymin) < self.pml_layer[2])
		mapPU = np.nonzero((self.ymax - self.y) < self.pml_layer[3])

		vmapPL = ids[mapPL]
		vmapPR = ids[mapPR]
		vmapPD = ids[mapPD]
		vmapPU = ids[mapPU]

		return list(set(np.concatenate((vmapPL,vmapPR,vmapPD,vmapPU))))

	def BuildPNonCon2D(self):
		# Find maximum requested polynomial order
		Nmax = np.max(self.Norder)
		# Mesh details
		VX = self.VX
		VY = self.VY
		e2v = self.e2v
		bc_wall = self.bc_wall
		bc_free = self.bc_free
		SrcDom = self.SrcDom
		GthDom = self.GthDom
		area_glob = self.area_glob
		K = self.K
		xmu = self.xmu
		xlambda = self.xlambda
		rho = self.rho

		pinfo = [Element() for __ in range(0,Nmax)]

		kmap = np.zeros((2,self.K), dtype=int)
		sk = 0
		# Perform a mini StartUp2D for the elements of each order
		print("Info: Initializing DG elements per approximation order...")
		for n in range(Nmax):
			# load N'th order polynomial nodes
			self.N = n+1
			self.Nfp = self.N+1
			self.Np = int(self.Nfp*(self.Nfp+1)/2)
			r,s = self.Nodes2D()
			self.r,self.s = xytors(r,s)
			# Find list of N'th order nodes on each face of the reference element
			ids1 = np.nonzero(np.abs(1+self.s) < self.NODETOL)[0]
			ids = np.argsort(self.r[ids1], kind='quicksort')
			ids1 = ids1[ids]
			ids2 = np.nonzero(np.abs(self.r+self.s) < self.NODETOL)[0]
			ids = np.argsort(-self.r[ids2], kind='quicksort')
			ids2 = ids2[ids]
			ids3 = np.nonzero(np.abs(1+self.r) < self.NODETOL)[0]
			ids = np.argsort(-self.s[ids3], kind='quicksort')
			ids3 = ids3[ids]
			self.Fmask = np.array([ids1,ids2,ids3])

			# Build reference element matrices
			self.V = Vandermonde2D(self.N,self.r,self.s)
			self.invV = la.inv(self.V)
			self.Dr,self.Ds = self.Dmatrices2D()
			self.LIFT = self.Lift2D()

			# store information for N'th order elements
			pinfo[n].Np = self.Np
			pinfo[n].Nfp = self.Nfp
			pinfo[n].Fmask = self.Fmask
			pinfo[n].r	= self.r
			pinfo[n].s	= self.s
			pinfo[n].Dr	= self.Dr
			pinfo[n].Ds	= self.Ds
			pinfo[n].LIFT = self.LIFT
			pinfo[n].V = self.V
			pinfo[n].invV = self.invV
			pinfo[n].massM = np.dot(self.invV.T,self.invV)
			# Find elements of polynomial order N
			ksN = np.nonzero(self.Norder==self.N)[0]
			self.K = len(ksN)
			if self.K > 0:
				# Use the subset of elements of order N
				self.e2v = e2v[:,ksN]
				self.bc_wall	= bc_wall[:,ksN]
				self.bc_free	= bc_free[:,ksN]
				self.SrcDom		= SrcDom[ksN]
				self.GthDom		= GthDom[:,ksN]
				self.area_glob	= area_glob[ksN]
				self.xmu		= xmu[ksN]
				self.xlambda	= xlambda[ksN]
				self.rho		= rho[ksN]

				kmap[0,ksN] = np.arange(self.K)
				kmap[1,ksN] = n
				imap = np.arange(sk,sk+self.K*self.Np).reshape((self.K,self.Np))

				# Build coordinates of all the nodes
				va,vb,vc = self.e2v
				self.x = 0.5*(-np.dot(VX[va].reshape(-1,1),[(self.r+self.s)])+np.dot(VX[vb].reshape(-1,1),[(1+self.r)])+np.dot(VX[vc].reshape(-1,1),[(1+self.s)]))
				self.y = 0.5*(-np.dot(VY[va].reshape(-1,1),[(self.r+self.s)])+np.dot(VY[vb].reshape(-1,1),[(1+self.r)])+np.dot(VY[vc].reshape(-1,1),[(1+self.s)]))

				# Calculate geometric factors
				self.rx,self.sx,self.ry,self.sy,self.J = self.GeometricFactors2D()
				self.nx,self.ny,self.sJ = self.Normals2D()
				self.Fscale = self.sJ/self.J[:,self.Fmask.flat]
				# Calculate element connections on this mesh
				self.e2e, self.e2f = self.tiConnect2D()
				self.mapM, self.mapP, self.vmapM, self.vmapP, self.vmapB, self.mapB = self.BuildMaps2D()
				self.mapW, self.vmapW = self.BuildWallBCMaps2D()
				self.mapF, self.vmapF = self.BuildFreeBCMaps2D()
				pinfo[n].mapW  		= self.mapW
				pinfo[n].vmapW 		= self.vmapW
				pinfo[n].mapF  		= self.mapF
				pinfo[n].vmapF 		= self.vmapF

				params = np.ones((self.K,self.Np))
				pinfo[n].xmu   		= self.xmu.reshape(-1,1)*params
				pinfo[n].xlambda   	= self.xlambda.reshape(-1,1)*params
				pinfo[n].rho   		= self.rho.reshape(-1,1)*params

				mapSource, srcArea = self.BuildBCSource(imap)
				pinfo[n].mapSource = mapSource
				pinfo[n].srcArea = srcArea
				self.BuildGthDom(imap)

				self.vmapPML = self.BuildBCPML(imap)
				pinfo[n].vmapPML = self.vmapPML

				# Compute triangulation of N'th order nodes on mesh
				triN = sp.Delaunay(np.array([self.r,self.s]).T).simplices
				#, qhull_options='Qt Qb Qc'
				alltri = triN
				for k in range(1,self.K):
					alltri = np.vstack((alltri,triN+k*self.Np))

				# Store geometric inforatmion in pinfo struct
				pinfo[n].rx = self.rx
				pinfo[n].sx = self.sx
				pinfo[n].ry = self.ry
				pinfo[n].sy = self.sy
				pinfo[n].nx = self.nx
				pinfo[n].ny = self.ny
				pinfo[n].sJ = self.sJ
				pinfo[n].J  = self.J
				pinfo[n].x  = self.x
				pinfo[n].y  = self.y
				pinfo[n].Fscale = self.Fscale
				pinfo[n].K = self.K
				pinfo[n].ks = ksN
				pinfo[n].tri = alltri
				# Store location of the N'th order nodes in a global vector
				pinfo[n].ids = imap
				sk = sk+self.K*self.Np
				pinfo[n].vx = self.VX[self.e2v]
				pinfo[n].vy = self.VY[self.e2v]
				pinfo[n].calc_area()

		# For each possible order
		for n1 in range(Nmax):
			# generate face L2projection matrices (from order N2 to order N1 face space)
			pinfo[n1].interpP = np.zeros(Nmax,dtype=object)
			for n2 in range(Nmax):
				# Set up sufficient Gauss quadrature to exactly perform surface integrals
				gz, gw = sf.j_roots(max(n1+2,n2+2),0,0)

				# All edges have same distribution (note special Fmask)
				rM = pinfo[n1].r[pinfo[n1].Fmask[0]]
				rP = pinfo[n2].r[pinfo[n2].Fmask[0]]

				# Build N2 to N1 projection matrices for '+' trace data
				interpM = np.dot(Vandermonde1D(n1+1,gz),la.inv(Vandermonde1D(n1+1,rM)))
				interpP = np.dot(Vandermonde1D(n2+1,gz[::-1]),la.inv(Vandermonde1D(n2+1,rP)))

				# Face mass matrix used in projection
				mmM = np.dot(interpM.T,np.dot(np.diag(gw),interpM))
				pinfo[n1].interpP[n2] = la.lstsq(mmM,np.dot(interpM.T,np.dot(np.diag(gw),interpP)))[0]

		# Generate neighbor information for all faces
		self.e2v = e2v
		self.area_glob = area_glob
		self.xmu = xmu
		self.xlambda = xlambda
		self.rho = rho
		self.e2e, self.e2f = self.tiConnect2D()

		# For each possible polynomial order
		for n1 in range(Nmax):
			# Create a set of indexing arrays, one for each possible neighbor order
			pinfo[n1].fmapM = [np.array([],dtype=int) for __ in range(Nmax)]
			pinfo[n1].vmapP = [np.array([],dtype=int).reshape(0,n2+2) for n2 in range(Nmax)]

			# Loop through all elements of order N1
			for k1 in range(pinfo[n1].K):
				# Find element in original mesh
				k1orig = pinfo[n1].ks[k1]
				# Check all it's faces
				for f1 in range(self.Nfaces):
					# Find neighboring element (i.e. it's order and location in N2 mesh)
					k2orig = self.e2e[f1,k1orig]
					f2     = self.e2f[f1,k1orig]
					k2     = kmap[0,k2orig]
					n2     = kmap[1,k2orig]

					# Compute location of face nodes of '-' trace
					idsM = k1*pinfo[n1].Nfp*self.Nfaces + f1*pinfo[n1].Nfp
					idsM = np.arange(idsM,idsM+pinfo[n1].Nfp)

					# Find location of volume nodes on '+' trace of (k1orig,f1)
					idsP = pinfo[n2].ids[k2,[pinfo[n2].Fmask[f2]]]
					# Store node locations in cell arrays
					pinfo[n1].fmapM[n2] = np.append(pinfo[n1].fmapM[n2], idsM)
					pinfo[n1].vmapP[n2] = np.concatenate((pinfo[n1].vmapP[n2], idsP), axis=0)

		return pinfo
