from scipy import linalg as la, special as sf, spatial as sp, interpolate as it
import numpy as np
import os
import operator

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
import matplotlib.tri as tri

class Element(object):
	"""docstring for Element"""
	def __init__(self):
		super(Element, self).__init__()
		self.elem 		= np.array([])
		self.vx			= np.array([])
		self.vy			= np.array([])
		self.mapW 		= np.array([])
		self.vmapW 		= np.array([])
		self.mapF  		= np.array([])
		self.vmapF 		= np.array([])
		self.mapSource 	= np.array([],dtype=int)
		self.rx 		= np.array([])
		self.sx 		= np.array([])
		self.ry 		= np.array([])
		self.sy 		= np.array([])
		self.nx 		= np.array([])
		self.ny 		= np.array([])
		self.sJ 		= np.array([])
		self.J  		= np.array([])
		self.x  		= np.array([])
		self.y  		= np.array([])
		self.Fscale 	= np.array([])
		self.K 			= 0
		self.ks 		= np.array([])
		self.tri 		= np.array([])
		self.ids 		= np.array([])
		self.area 		= np.array([])

	def calc_area(self):
		len1 = np.sqrt((self.vx[0]-self.vx[1])**2+(self.vy[0]-self.vy[1])**2)
		len2 = np.sqrt((self.vx[1]-self.vx[2])**2+(self.vy[1]-self.vy[2])**2)
		len3 = np.sqrt((self.vx[2]-self.vx[0])**2+(self.vy[2]-self.vy[0])**2)
		sper = (len1 + len2 + len3)/2.0
		self.area = np.sqrt(sper*(sper-len1)*(sper-len2)*(sper-len3))

class MeshReader(object):
	"""docstring for MeshReader"""
	def __init__(self, **kwargs):
		super(MeshReader, self).__init__()
		self.path = kwargs['mesh_path']
		self.dtype = kwargs['mesh_type']
		self.sigma_smooth = kwargs['sigma_smooth']
		# find source position
		self.src = kwargs['source']
		# find gather position
		self.gather = kwargs['gather']

		if self.dtype == 'neu':
			self.Nv, self.VX, self.VY, self.K, self.EToV = self.__neu()
		elif self.dtype == 'msh':
			self.Nv, self.VX, self.VY, self.K, self.EToV = self.__msh()
		else:
			print("Error: can't find this mesh format.")
			exit(-1)

		self.xmin = np.min(self.VX)
		self.xmax = np.max(self.VX)
		self.ymin = np.min(self.VY)
		self.ymax = np.max(self.VY)

		self.pml_layer = kwargs['pml_layer']

		self.maps = self.FindSource()

		self.mapg = self.FindGather()

		# Build source rounding elements
		self.surrounding = self.getAroundSource()

		self.src_cells = len(self.surrounding)

		self.__plot()

	def __plot(self):
		fig, ax = plt.subplots()
		mesh = tri.Triangulation(self.VX,self.VY)
		plt.triplot(mesh, lw=0.5, color='blue')
		plt.scatter(self.src[0],self.src[1],color='red', marker='*', s=50)
		for gth in self.gather:
			gx,gy = gth
			plt.scatter(gx,gy, color='blue', marker='v', s=50)

		for k in self.surrounding:
			cx = np.sum(self.VX[self.EToV[:,k]])/3
			cy = np.sum(self.VY[self.EToV[:,k]])/3
			ax.text(cx, cy, ('{}').format(k), fontsize=9, bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

		plt.gca().invert_yaxis()
		ax.set_aspect('equal')
		plt.show()

	def getAroundSource(self):
		mat = np.array(list(map(lambda a: np.sqrt((np.sum(self.VX[a])/3 - self.src[0])**2 + (np.sum(self.VY[a])/3 - self.src[1])**2), self.EToV.T)))
		ids = np.nonzero(mat <= 3*self.sigma_smooth)[0]
		if len(ids) == 0:
			ids = np.array([self.maps])

		return ids

	def barycentric(self,k,p):
		# get vertices of the elements
		v1 = np.array([self.VX[self.EToV[0,k]],self.VY[self.EToV[0,k]]])
		v2 = np.array([self.VX[self.EToV[1,k]],self.VY[self.EToV[1,k]]])
		v3 = np.array([self.VX[self.EToV[2,k]],self.VY[self.EToV[2,k]]])
		a = v1-v2
		b = v3-v2
		s =  p-v2
		l1 = np.cross(s,a)/np.cross(b,a)
		l3 = np.cross(s,b)/np.cross(a,b)
		l2 = 1.0 - l1 - l3
		return l1+0,l2+0,l3+0

	def FindSource(self):
		mapS = None
		k=0
		while k < self.K:
			l1,l2,l3 = self.barycentric(k,self.src)
			# condition to calculate the element that contains the source
			if((l1 >= 0.0 and l1 <= 1.0) and (l2 >= 0.0 and l2 <= 1.0) and (l3 >= 0.0 and l3 <= 1.0)):
				mapS = k
				break
			k=k+1

		if mapS is None:
			print('Error: Source outside the computacional domain.')
			exit(1)
		else:
			print('Info: Source found in triangle {}.'.format(mapS))

		return mapS

	def FindGather(self):
		mapG = [[],[]]
		k=0
		while k < self.K:
			for i,gth in enumerate(self.gather):
				l1,l2,l3 = self.barycentric(k,gth)
				# condition to calculate the element that contains the source
				if((l1 >= 0.0 and l1 <= 1.0) and (l2 >= 0.0 and l2 <= 1.0) and (l3 >= 0.0 and l3 <= 1.0)):
					mapG[0].append(i)
					mapG[1].append(k)
			
			if len(mapG[1]) == len(self.gather):
				break

			k=k+1

		result = np.array(sorted(np.array(mapG).T, key=operator.itemgetter(0), reverse=False)).T[1]


		if len(mapG[1]) == 0:
			print('Error: All gathers are outside the computacional domain.')
			exit(1)
		else:
			print('Info: Gather found in triangle {}.'.format(result))

		return result

	# Function to read the mesh
	def __neu(self):
		if os.path.exists(self.path):
			with open(self.path, 'rt') as f:
				data = f.readlines()
				Nv,K  = np.array(data[6].split()[:2]).astype(int)
				VX,VY = np.array(list(map(lambda a: a.split()[1:], data[9:Nv+9])), dtype=float).T
				EToV  = np.array(list(map(lambda a: a.split()[3:], data[Nv+11:K+Nv+11])), dtype=int).T-1

			print("Mesh file {} loadded successfully. ({} E, {} V)".format(self.path.split('/')[-1], K, Nv))

			return Nv, VX, VY, K, EToV
		else:
			print("Error: Especified file '{}' does not exists.".format(self.path))
			exit(1)

	def __msh(self):
		if os.path.exists(self.path):
			with open(self.path, 'rt') as f:
				data = f.read().splitlines()
				# line identification for nodes
				nodeId = data.index('$Nodes')+1
				# line identification for elements
				elementId = data.index('$Elements')+1
				# Number of vertex
				Nv = int(data[nodeId])
				# number of elements
				K = int(data[elementId])

				# Nodes array
				nodes 	 = np.array([list(map(float, a.split(' ')[1:3])) for a in data[nodeId+1:nodeId+Nv+1]])
				VX = nodes[:,0]
				VY = nodes[:,1]
				# Elements array (check if elements are of type 2)
				EToV = np.array([list(map(int, a.split(' ')[1:]))[-3:] for a in data[elementId+1:elementId+K+1] if list(map(int,   a.split(' ')[1:]))[0] == 2]).T-1
				# Update number of elements
				K = np.size(EToV,1)

			print("Mesh file {} loadded successfully ({} E, {} V).".format(self.path.split('/')[-1], K, Nv))

			return Nv, VX, VY, K, EToV
		else:
			print("Error: Especified file '{}' does not exists.".format(self.path))
			exit(1)