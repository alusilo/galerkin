from .nodalDG import *
import sys

def Grad2D(info,u):
	ur = np.dot(u,info.Dr.T)
	us = np.dot(u,info.Ds.T)
	ux = info.rx*ur + info.sx*us
	uy = info.ry*ur + info.sy*us

	return ux,uy

class WaveDrive2D(NodalDG2D):
	"""docstring for WaveDrive2D"""
	def __init__(self, **kwargs):
		super(WaveDrive2D, self).__init__(**kwargs)
		self.source_path = kwargs['src_path']
		self.freq = kwargs['freq']
		self.vp = kwargs['p_velocity']
		self.vs = 1.0#3464.#self.vp/1.8#kwargs['s_velocity']
		self.rho = 1.0#2670.#310.0*self.vp**0.25#kwargs['density']
		
		# lame parameters
		self.param1 = self.rho*self.vs*self.vs
		self.param2 = self.rho*self.vp*self.vp - 2*self.param1
		
		# resolution
		self.pixel_size = kwargs['pixel_size']
		self.frame_time = kwargs['frame_time']
		
		# stress tensor
		self.stress = kwargs['stress'].reshape(-1,1)

		# data arrays and initial conditions
		self.x = np.array([])
		self.y = np.array([])
		for pinf in self.pinfo:
			self.x = np.append(self.x,pinf.x)
			self.y = np.append(self.y,pinf.y)

		# fig, ax = plt.subplots()
		# ax.set_aspect('equal')
		# mesh = tri.Triangulation(self.VX,self.VY)
		# plt.triplot(mesh, lw=0.5, color='blue')
		# for obj in self.pinfo:
		# 	if len(obj.vmapW) > 0:
		# 		plt.scatter(self.x[obj.vmapW],self.y[obj.vmapW], color='red', marker='o', s=50)
		
		# plt.show()
		nmode = mmode = 1
		self.Ux  = np.sin(nmode*np.pi*self.x/2)*np.sin(mmode*np.pi*self.y/2)
		self.Uy  = np.sin(nmode*np.pi*self.x/2)*np.sin(mmode*np.pi*self.y/2)
		self.Sxx = np.sin(nmode*np.pi*self.x/2)*np.sin(mmode*np.pi*self.y/2)
		self.Syy = np.sin(nmode*np.pi*self.x/2)*np.sin(mmode*np.pi*self.y/2)
		self.Sxy = np.sin(nmode*np.pi*self.x/2)*np.sin(mmode*np.pi*self.y/2)
		
		# Set simulation time
		self.FinalTime = kwargs['duration']
		self.source_dt = kwargs['source_dt']
		# Solve the problem
		self.WavePNonCon2D()

	def updateSource(self,t,dt):
		for obj in self.pinfo:
			sx,sy = self.src
			x = self.x[obj.mapSource]
			y = self.y[obj.mapSource]
			A = self.wavelet[t]
			Mxx,Myy,Mxy,Fx,Fy = self.stress
			dist = np.exp(-((x-sx)**2 + (y-sy)**2)/(2*self.sigma_smooth**2))
			self.Sxx[obj.mapSource] += Mxx*A*dt*dist/(2*np.pi*(self.sigma_smooth)**2)
			self.Syy[obj.mapSource] += Myy*A*dt*dist/(2*np.pi*(self.sigma_smooth)**2)
			self.Sxy[obj.mapSource] += Mxy*A*dt*dist/(2*np.pi*(self.sigma_smooth)**2)
			self.Ux[obj.mapSource]  += Fx*A*dt*dist/(2*np.pi*(self.sigma_smooth)**2)
			self.Uy[obj.mapSource]  += Fy*A*dt*dist/(2*np.pi*(self.sigma_smooth)**2)


	def abc(self,x,y):
		res = np.ones(x.shape)
		a = 0.001
		ids = np.nonzero(x <= self.xmin+self.pml_layer)
		res[ids] *= np.exp(-(a*(self.xmin+self.pml_layer-x[ids]))**2)
		ids = np.nonzero(x >= self.xmax-self.pml_layer)
		res[ids] *= np.exp(-(a*(self.xmax-self.pml_layer-x[ids]))**2)
		ids = np.nonzero(y <= self.ymin+self.pml_layer)
		res[ids] *= np.exp(-(a*(self.ymin+self.pml_layer-y[ids]))**2)
		ids = np.nonzero(y >= self.ymax-self.pml_layer)
		res[ids] *= np.exp(-(a*(self.ymax-self.pml_layer-y[ids]))**2)

		return res

	def WavePNonCon2D(self):
		# Runge-Kutta residual storage  
		resUx  = np.zeros(len(self.Ux),dtype=np.float32)
		resUy  = np.zeros(len(self.Ux),dtype=np.float32)
		resSxx = np.zeros(len(self.Ux),dtype=np.float32)
		resSyy = np.zeros(len(self.Ux),dtype=np.float32)
		resSxy = np.zeros(len(self.Ux),dtype=np.float32)

		# stability parameter
		self.dt = (self.vs/self.vp)*(1/(40*self.Nmax*self.Nmax*self.freq))
		
		print("dt = {}[sec]".format(self.dt))
		print("computational domain = [({},{}), ({},{})]".format(self.xmin,self.xmax, self.ymin,self.ymax))
		print("minimun polynomials order = {}".format(self.Nmin))
		print("maximun polynomials order = {}".format(self.Nmax))
		print("p-wave velocity = {}[m/s]".format(self.vp))
		print("s-wave velocity = {}[m/s]".format(self.vs))
		print("density = {}[Kg/m^3]".format(self.rho))
		print("maximun frequency = {}[Hz]".format(self.freq))
		print("source position = {}".format(self.src))
		print("source cells spreading = {}".format(self.src_cells))
		print("gather(s) position = {}".format(self.gather))

		if self.dt > self.FinalTime or self.dt < 0:
			print("Error: simulation time must be greater than sampling time, and this must be positive.")
			exit(1)

		# total simulation time steps
		nsteps = int(np.ceil(self.FinalTime/self.dt)+1)

		# simulation time array
		time = np.array([i*self.dt if i*self.dt < self.FinalTime else self.FinalTime for i in range(nsteps)])

		# Total simulation snaphots
		snapdt = int(np.ceil(self.FinalTime/self.frame_time)+1)
		snapdt = int(np.ceil(nsteps/snapdt)+1)
		
		# number of pixels in each snapshot
		h_pixels = int(np.ceil((self.xmax-self.xmin)/self.pixel_size))
		v_pixels = int(np.ceil((self.ymax-self.ymin)/self.pixel_size))
		
		# interval to print out simulation progress
		pstep = int(np.ceil(nsteps/100))

		# read source information
		self.wavelet = self.ReadSource(time)

		# array to storage the information of perturbation in several points
		traces = np.zeros((len(self.gather),nsteps))
		# movie array
		movie = []

		# init figure
		fig, ax = plt.subplots(1, 2, sharex=False, sharey=False)
		ax[0].set_aspect('equal')
		ax[1].set_aspect('auto')
		
		# initial time
		t0 = 0.0
		# outer time step loop
		for n,t in enumerate(time):
			for INTRK in range(5):
				# compute right hand side of TM-mode Maxwell's equations
				rhsSxx, rhsSyy, rhsSxy, rhsUx, rhsUy = self.WavePNonConRHS2D()
				# initiate and increment Runge-Kutta residuals
				resSxx = rk4a[INTRK]*resSxx + (t-t0)*rhsSxx
				resSyy = rk4a[INTRK]*resSyy + (t-t0)*rhsSyy
				resSxy = rk4a[INTRK]*resSxy + (t-t0)*rhsSxy
				resUx  = rk4a[INTRK]*resUx  + (t-t0)*rhsUx
				resUy  = rk4a[INTRK]*resUy  + (t-t0)*rhsUy
				# update fields
				self.Sxx = self.Sxx + rk4b[INTRK]*resSxx
				self.Syy = self.Syy + rk4b[INTRK]*resSyy
				self.Sxy = self.Sxy + rk4b[INTRK]*resSxy
				self.Ux  = self.Ux  + rk4b[INTRK]*resUx
				self.Uy  = self.Uy  + rk4b[INTRK]*resUy
			
			# ABC conditions
			# self.Sxx[self.vmapPML] *= self.abc(self.x[self.vmapPML],self.y[self.vmapPML])
			# self.Syy[self.vmapPML] *= self.abc(self.x[self.vmapPML],self.y[self.vmapPML])
			# self.Sxy[self.vmapPML] *= self.abc(self.x[self.vmapPML],self.y[self.vmapPML])
			# self.Ux[self.vmapPML]  *= self.abc(self.x[self.vmapPML],self.y[self.vmapPML])
			# self.Uy[self.vmapPML]  *= self.abc(self.x[self.vmapPML],self.y[self.vmapPML])

			# print progress
			if n%pstep == 0 or n == nsteps-1:
				sys.stdout.write('\r')
				sys.stdout.write("[%-30s] %d%%" % ('='*int(30*n/(nsteps-1)), 100*(n/(nsteps-1))))
				sys.stdout.flush()

			# store simulation snapshot
			if n%snapdt == 0 or n == nsteps-1:
				# performe a linear interpolation funtion of the whole computational domain
				points = np.array([self.x,self.y]).T
				X, Y = np.mgrid[self.xmin:self.xmax:h_pixels*1j, self.ymin:self.ymax:v_pixels*1j]
				Z0 = it.griddata(points, self.Sxy,(X,Y), method='cubic', fill_value=0)
				movie.append(Z0)

			# store gather
			for i,gth in enumerate(self.gather):
				#l1,l2,l3 = self.barycentric(self.mapg[i],gth)
				#r,s = l2*np.array([[-1.],[-1.]])+l3*np.array([[1.],[-1.]])+l1*np.array([[-1.],[1.]])
				#V = Vandermonde2D(self.vmapGO[i]+1,r,s)
				#inter_matrix = np.dot(V,self.pinfo[self.vmapGO[i]].invV)
				points = np.array([self.x[self.vmapG[i]],self.y[self.vmapG[i]]]).T
				traces[i,n] = it.griddata(points, self.Ux[self.vmapG[i]], (gth[0],gth[1]), method='cubic', fill_value=0)#np.dot(inter_matrix, self.Ux[self.vmapG[i]])[0]
			
			# update perturbation
			#self.updateSource(n,t-t0)

			t0 = t
		
		movie = np.array(movie)
		np.save('../resources/output/movie.npy', movie)
		np.save('../resources/output/traces.npy', np.concatenate(([time], traces)))

		return

	def ReadSource(self,time):
		# reading source values from file
		if os.path.exists(self.source_path):
			with open(self.source_path) as f:
				data = f.readlines()
				source = np.array([np.array(line).astype(float) for line in data])
				ts = np.array([i*self.source_dt for i in range(len(source))])
		else:
			print("Error: Especified file '{}' does not exists.".format(self.source_path))
			exit(1)
		
		wout = np.zeros(len(time))
		# time iteration
		for n,t in enumerate(time):
			# check time sub domain associated to points in file
			for p in range(1,np.size(source,0)):
				if(t >= ts[p-1] and t <= ts[p]):
					t0 = ts[p-1]
					t1 = ts[p]
					w0 = source[p-1]
					w1 = source[p]
					break
				else:
					t0 = 0.0
					t1 = 0.0
			# Linear interpolation of the source values
			if(t1-t0 == 0.0):
				wout[n] = 0.0
			else:
				wout[n] = w0+(w1-w0)*(t-t0)/(t1-t0)

		plt.grid()
		plt.title("Wavelet")
		plt.xlabel("t[s]")
		plt.ylabel("Amplitude")
		plt.plot(time,wout)
		plt.show()
		
		return wout

	def WavePNonConRHS2D(self):
		# Initialize storage for right hand side residuals
		rhsUx  = np.zeros(len(self.Ux),dtype=np.float32)
		rhsUy  = np.zeros(len(self.Ux),dtype=np.float32)
		rhsSxx = np.zeros(len(self.Ux),dtype=np.float32)
		rhsSyy = np.zeros(len(self.Ux),dtype=np.float32)
		rhsSxy = np.zeros(len(self.Ux),dtype=np.float32)

		# For each possible polynomial order
		for N in range(len(self.pinfo)):
			# Extract information for this polynomial order
			pinf = self.pinfo[N]
			K = pinf.K
			# Check to see if any elements of this order exist
			if K > 0:
				# Find location of N'th order nodes
				ids = pinf.ids
				Fmask = pinf.Fmask
				# Extract N'th order nodes
				SxxN = self.Sxx[ids].reshape(ids.shape)
				SyyN = self.Syy[ids].reshape(ids.shape)
				SxyN = self.Sxy[ids].reshape(ids.shape)
				UxN  = self.Ux[ids].reshape(ids.shape)
				UyN  = self.Uy[ids].reshape(ids.shape)
				
				# Extract '-' traces of N'th order nodal data
				SxxM = SxxN[:,Fmask.flat]
				SyyM = SyyN[:,Fmask.flat]
				SxyM = SxyN[:,Fmask.flat]
				UxM  = UxN[:,Fmask.flat]
				UyM  = UyN[:,Fmask.flat]

				# Storage for '+' traces
				SxxP = np.zeros(SxxM.shape)
				SyyP = np.zeros(SyyM.shape)
				SxyP = np.zeros(SxyM.shape)
				UxP  = np.zeros(UxM.shape)
				UyP  = np.zeros(UyM.shape)

				# For each possible order
				for N2 in range(len(self.pinfo)):
					# Check to see if any neighbor nodes of this order were located
					if len(pinf.fmapM[N2]) > 0:
						# L2 project N2'th order neighbor data onto N'th order trace space
						interp = pinf.interpP[N2]
						fmapM  = pinf.fmapM[N2]
						vmapP  = pinf.vmapP[N2]
						SxxP.flat[fmapM] = np.dot(self.Sxx[vmapP].reshape(vmapP.shape),interp.T)
						SyyP.flat[fmapM] = np.dot(self.Syy[vmapP].reshape(vmapP.shape),interp.T)
						SxyP.flat[fmapM] = np.dot(self.Sxy[vmapP].reshape(vmapP.shape),interp.T)
						UxP.flat[fmapM]  = np.dot(self.Ux[vmapP].reshape(vmapP.shape),interp.T)
						UyP.flat[fmapM]  = np.dot(self.Uy[vmapP].reshape(vmapP.shape),interp.T)

				# Compute jumps of trace data at faces
				dSxx = SxxM - SxxP
				dSyy = SyyM - SyyP
				dSxy = SxyM - SxyP
				dUx  = UxM - UxP
				dUy  = UyM - UyP

				# Apply PEC boundary condition at wall boundary faces
				dSxx.flat[pinf.mapW] = 0.
				dSyy.flat[pinf.mapW] = 0.
				dSxy.flat[pinf.mapW] = 0.
				dUx.flat[pinf.mapW]  = 0.
				dUy.flat[pinf.mapW]  = 0.

				# evaluate upwind fluxes
				fluxSxx = -(pinf.ny*self.param2*dUy + pinf.nx*(self.param2+2*self.param1)*dUx)
				fluxSyy = -(pinf.nx*self.param2*dUx + pinf.ny*(self.param2+2*self.param1)*dUy)
				fluxSxy = -self.param1*(pinf.ny*dUx + pinf.nx*dUy)
				fluxUx  = -(1/self.rho)*(pinf.nx*dSxx + pinf.ny*dSxy)
				fluxUy  = -(1/self.rho)*(pinf.ny*dSyy + pinf.nx*dSxy)
				
				# local derivatives of fields
				Sxx_x, Sxx_y = Grad2D(pinf,SxxN)
				Syy_x, Syy_y = Grad2D(pinf,SyyN)
				Sxy_x, Sxy_y = Grad2D(pinf,SxyN)
				Ux_x,  Ux_y  = Grad2D(pinf,UxN)
				Uy_x,  Uy_y  = Grad2D(pinf,UyN)
				
				# compute right hand sides of the PDE's
				rhsSxx[ids] = (self.param2+2*self.param1)*Ux_x + self.param2*Uy_y + np.dot(pinf.Fscale*fluxSxx,pinf.LIFT)/2
				rhsSyy[ids] = (self.param2+2*self.param1)*Uy_y + self.param2*Ux_x + np.dot(pinf.Fscale*fluxSyy,pinf.LIFT)/2
				rhsSxy[ids] = self.param1*(Uy_x + Ux_y) + np.dot(pinf.Fscale*fluxSxy,pinf.LIFT)/2
				rhsUx[ids]  = (1/self.rho)*(Sxx_x + Sxy_y) + np.dot(pinf.Fscale*fluxUx,pinf.LIFT)/2
				rhsUy[ids]  = (1/self.rho)*(Sxy_x + Syy_y) + np.dot(pinf.Fscale*fluxUy,pinf.LIFT)/2

		return rhsSxx,rhsSyy,rhsSxy,rhsUx,rhsUy

M = np.array([0., 0., 0., 0., 1.])
# Galerkin Object
obj = WaveDrive2D(
	mesh_path	 = '../resources/mesh/small.msh', # Mesh file path
	mesh_type	 = 'msh', # type of mesh [neu, msh]
	src_path 	 = '../resources/source/source', # Source path
	source 		 = (1.,1.),#(2800.,200.), # Source position (x,y)
	sigma_smooth = 0.01, # Sigma spatial support
	gather		 = [(0.5,0.5)],#[(3000.,0.)], # array of positions (tuples) to record information
	max_order	 = 8, # maximun order
	min_order	 = 8, # minimun order
	freq		 = np.sqrt(2),#7.5, # frequency
	p_velocity	 = 2.,#6000., # p-wave velocity
	duration	 = 1.0, # simulation time
	pml_layer	 = 0.1,#300., # PML length from border
	pixel_size	 = 0.01, # movie resolution
	frame_time	 = 0.01, # sampling time between frames
	source_dt	 = 0.0001, # Source file sampling
	stress		 = M, # initial stress tensor: [Mxx, Myy, Mxy]
)
'''
	relations
	p-wave and dendity
	rho = a*vp^m, a = 0.31 and m = 0.25 rho in [g/cm^3] and vp in [m/s]
'''
#1580000