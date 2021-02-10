from numpy import *
from numpy import sum as Sum
from numpy import transpose as Trans
from numpy import diagflat as diag
from numpy.linalg import *
from scipy.special import *

set_printoptions(precision=4,suppress=True)

def vec(v,w):
	return array(v[w[0]:w[len(w)-1]+1])

def getm(v,w):
	return array([v[i] for i in w])

def jacobi_gauss_lobatto_points(alpha,beta,N):
	if (N==1):
		x = array([-1, 1])
		return x
	norm = 1./sqrt(2./(2.*N))
	xint, w = j_roots(N-1,alpha,beta)
	x = append(append([-1.],xint),[1.])
	return x

def jacobi_normalized_polynomials(N,alpha,beta,x):
	gamma0 = (2.**(alpha+beta+1)/(2*N+alpha+beta+1))*(gamma(N+alpha+1)*gamma(N+beta+1)/(gamma(N+alpha+beta+1)*factorial(N)))
	norm = 1./sqrt(gamma0)
	return norm*eval_jacobi(N,alpha,beta,x)

def gradient_jacobi_polynomials(r,alpha,beta,N):
	dP = zeros((len(r),1))
	if N==0:
		return dP
	else:
		return sqrt(N*(N+alpha+beta+1))*jacobi_normalized_polynomials(N-1,alpha+1,beta+1,r)

def vandermonde_matrix_1d(N,r):
	vm = zeros((len(r),N+1))
	for j in range(0,N+1):
		vm[j] = jacobi_normalized_polynomials(j,0,0,r)

	return Trans(vm)

def gradient_vandermonde_matrix_1d(N,r):
	dvr = zeros((len(r),N+1))
	for i in range(0,N+1):
		dvr[i] = Trans(gradient_jacobi_polynomials(r,0,0,i))
	
	return Trans(dvr)

def gradient_matrix_1d(N,r,V):
	Vr = gradient_vandermonde_matrix_1d(N,r)
	return dot(Vr,inv(V))
