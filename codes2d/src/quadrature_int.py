import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from nodalDG import *
from cubatureData2D import cubatureData2D

#eps = 1e-12

def f(x,y):
	return np.exp(-2j*np.pi*(x+y)/5).real

def area(vertices):
	(x1,y1),(x2,y2),(x3,y3) = vertices
	return np.abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2)

# triangle vertices
vertices = np.array([
	[1, 2],
	[5, 4],
	[2, 6]
])

# plot triangle
triangle = np.append(vertices,[vertices[0]], axis=0)
plt.plot(triangle[:,0],triangle[:,1])

# plot text on vertex of the element
for (vx,vy) in vertices:
	plt.scatter(vx,vy, c='b')
	plt.axes().text(vx,vy,'P({},{})'.format(vx,vy))

plt.axes().set_title(r'Elemento sobre sistema global ($x$, $y$)')
plt.axes().set_xlabel(r'$x$')
plt.axes().set_ylabel(r'$y$')
plt.axes().set_aspect('equal')
plt.axes().grid()
plt.show()

# Element approximation order
N = 4
# Number of support points
Np = int((N+1)*(N+2)/2)
# Points distribution on element
x,y = Nodes2D(N)

# equilateral triangle vertices
v_equi = np.array([
	[-1, np.min(y)],
	[ 1, np.min(y)],
	[ 0, np.max(y)]
])
# plot equilateral triangle with quadrature points
equi_tri = np.append(v_equi,[v_equi[0]], axis=0)
plt.plot(equi_tri[:,0],equi_tri[:,1])
plt.scatter(x,y)
plt.axes().set_title(r"Puntos de soporte deformados sobre triangulo isoseles, sistema ($x'$, $y'$)")
plt.axes().set_xlabel(r"$x'$")
plt.axes().set_ylabel(r"$y'$")
plt.axes().set_aspect('equal')
plt.axes().grid()
plt.show()

# change to reference coordinates
r,s = xytors(x,y)

# fmask = np.array([
# 	np.nonzero(abs(1+s) < eps)[0],
# 	np.nonzero(abs(r+s) < eps)[0],
# 	np.nonzero(abs(1+r) < eps)[0]])

# vertices reference triangle
v_ref = np.array([
	[-1,-1],
	[ 1,-1],
	[-1, 1]
])

# calculating Jacobian
(x1,y1),(x2,y2),(x3,y3) = vertices
J = np.abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/4)

# plot reference coordinates
ref_tri = np.append(v_ref,[v_ref[0]], axis=0)
plt.plot(ref_tri[:,0],ref_tri[:,1])
plt.scatter(r,s)
plt.axes().set_title(r"Elemento y puntos de soporte sobre sistema estandar ($r$, $s$)")
plt.axes().set_xlabel(r"$r$")
plt.axes().set_ylabel(r"$s$")
plt.axes().set_aspect('equal')
plt.axes().grid()
plt.show()

# plot global coordinates
x,y = 0.5*(-np.dot(vertices[0].reshape(-1,1),[r+s])+np.dot(vertices[1].reshape(-1,1),[1+r])+np.dot(vertices[2].reshape(-1,1),[1+s]))
plt.plot(triangle[:,0],triangle[:,1])
plt.scatter(x,y)
plt.axes().set_title(r'Puntos de soporte sobre sistema global ($x$, $y$)')
plt.axes().set_xlabel(r'$x$')
plt.axes().set_ylabel(r'$y$')
plt.axes().set_aspect('equal')
plt.axes().grid()
plt.show()

# Vandermonde matrix
V = Vandermonde2D(N,r,s)

# Vandermonde derivatives
Vr,Vs = GradVandermonde2D(N,r,s)

# Polynomial by grade
grade = 0
function = f(x,y)
P_n = V[:,grade]
fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x,y,function, s=25, c='r')
ax.plot_trisurf(x,y,function,antialiased=True)
#ax.plot_trisurf(r,s,dwr,antialiased=True,shade=True)
ax.view_init(50, -210)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('amplitude')

# calculate Vandermonde inverse
invV = la.inv(V)
# which points to be interpolated?
samples = 100
pr = -1.0*np.ones(samples)
ps = np.linspace(-1,1,samples)
px,py = 0.5*(-np.dot(vertices[0].reshape(-1,1),[pr+ps])+np.dot(vertices[1].reshape(-1,1),[1+pr])+np.dot(vertices[2].reshape(-1,1),[1+ps]))
# calculate polynomial basis on points to be interpolated
phi_n = Vandermonde2D(N,pr,ps)
# calculate interpolating Lagrange polynomials
interp = np.dot(phi_n,invV)
# interpolate points
pz = np.dot(interp,function)

#ax.scatter(pr,ps,pz,s=10,c='y')
#ax.plot(px,py,pz,linewidth=4,c='y')

plt.show()

# calculate area
A = area(vertices)
# the triangle area is 2 times the Jacobian of the transformation
print('area={}, 2J = {}'.format(A,2*J))

cub = cubatureData2D(N=N)

cubx,cuby = 0.5*(-np.dot(vertices[0].reshape(-1,1),[cub.r+cub.s])+np.dot(vertices[1].reshape(-1,1),[1+cub.r])+np.dot(vertices[2].reshape(-1,1),[1+cub.s]))
# calculate numerical integration
integral = J*np.dot(cub.w,f(cubx,cuby))
print(integral)
