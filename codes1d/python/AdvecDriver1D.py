from Mesh import *
from Jacobi import *
from Advec1D import *

import matplotlib.pyplot as plt
from matplotlib import animation

# Parameters definition
NODETOL = 1e-10
# domain extension
xmin = 0.0
xmax = 6.28
# Subdomains
K = 10
# Vertices
Nv = K+1
'''
return number of vertex, equi-spaced grid points, elements to vertex
VX are the points in the mesh and EToV is a Kx2 matrix, the first row
is the vertex index to the left of each element K and the second row
is the vertex index to the right of each element K.
'''
VX, EToV = MeshGen1D(xmin,xmax,K)
# Aproximation order
N=1
# number of points
Np = N+1
# Number of points per face
Nfp = 1
# Number of faces
Nfaces=2
# Jacobi Gauss Lobatto points
r = jacobi_gauss_lobatto_points(1,1,N)
# Vandermonde Matrix
V = vandermonde_matrix_1d(N, r)
# Matrix derivatives
Dr = gradient_matrix_1d(N, r, V)
# M = inv(dot(V,Trans(V)))
# print dot(M,Dr)
# LIFT = Lift1D()
Emat = zeros((Np,Nfaces*Nfp))
Emat[0,0] = 1.0
Emat[Np-1,1] = 1.0
#print Emat
LIFT = dot(V,(dot(Trans(V),Emat)))
va = Trans(EToV)[0]
vb = Trans(EToV)[1]
x1 = dot(ones((N+1,1)),array([VX[:len(va)]]))
x2 = 0.5*dot(Trans(array([r+1])),array([vec(VX,vb)-vec(VX,va)]))
x = x1+x2
rx,J = GeometricFactors1D(x,Dr)
fmask1 = nonzero(abs(r+1) < NODETOL)
fmask2 = nonzero(abs(r-1) < NODETOL)
Fmask  = array([fmask1[0][0],fmask2[0][0]])

Fx = getm(x,Fmask)
#nx = Normals1D()

nx = zeros((Nfp*Nfaces, K))
nx[0] = -1.0
nx[Nfp*Nfaces-1] = 1.0
Fscale = 1./getm(J,Fmask);
EToE, EToF = Connect1D(EToV,Nfaces,K,Nv)

vmapM, vmapP, vmapI, mapI, vmapB, mapB, mapO = BuildMaps1D(K,Np,Nfp,Nfaces,EToE,EToF,Fmask,x,NODETOL)

u = sin(x)

FinalTime = 30.
data = Advec1D(Np,K,x,rx,nx,Dr,u,LIFT,Fscale,Nfp,Nfaces,vmapM,vmapP,vmapI,mapI,mapO,FinalTime)

# Create animation
fig = plt.figure()
ax = plt.axes(xlim=(xmin, xmax), ylim=(-2, 2))
line, = ax.plot([], [], 'bo')

def animate(i):
    t = x.ravel()
    y = data[i].ravel()
    line.set_data(t, y)
    return line,

anim = animation.FuncAnimation(fig, animate,
                               frames=len(data), interval=100, blit=True)

anim.save('basic_animation.mp4', fps=5, extra_args=['-vcodec', 'libx264'])

plt.show()