from galerkin.wave2D import *

# Optional: Use Intel MKL for optimized BLAS if available
try:
    import mkl
    mkl.set_num_threads(4)
except ImportError:
    pass  # Continue without MKL optimizations

PROJECT_NAME 		= 'erase'

M = (1,1,0)
F = (0,0)
# Galerkin Object
obj = WaveDrive2D(
	project 	 = PROJECT_NAME,
	mesh_file	 = '../../resources/mesh/tri/tlr.ele', # Mesh file path
	param_file	 = '../../resources/mesh/tri/tlr.param', # Mesh file path
	src_position = (2146.76, 383.27),#(2144.,400.),# # Source position (x,y)
	src_smooth 	 = 50., # Sigma spatial support
	gather		 = [(3200.,0.01), (2500.,1100.), (3400.,500.)],#(3350.,0.01),(3600.,0.01)], # array of positions (tuples) to record information
	source_order = 3, # maximun order
	order	 	 = 2, # minimun order
	src_freq	 = 7.5, # frequency
	src_delay	 = 1/7.5, # frequency
	duration	 = 0.5, # simulation time
	pml_layer	 = (300., 300., 0., 300.), # PML length from border (xmin, xmax, ymin, ymax)
	pml_coef 	 = 0.0015,
	pixel_size	 = 10., # movie resolution
	tpf	 	 	 = 0.01, # sampling time between frames
	stress		 = M, # initial stress tensor: [Mxx, Myy, Mxy]
	displacement = F, # Velocity excitation [Fx, Fy]
)
'''
	relations
	p-wave and density
	rho = a*vp^m, a = 0.31 and m = 0.25 rho in [g/cm^3] and vp in [m/s]
'''
#1580000
