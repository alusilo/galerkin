from wave2D import *

M = np.array([1,1,0,0,0])
# Galerkin Object
obj = WaveDrive2D(
	mesh_path	 = '../resources/mesh/new-mesh2.msh', # Mesh file path
	mesh_type	 = 'msh', # type of mesh [neu, msh]
	src_path 	 = '../resources/source/gauss_rtp5.src', # Source path
	source 		 = (2800.,1800.), # Source position (x,y)
	sigma_smooth = 60., # Sigma spatial support
	gather		 = [(1790.,0.01), (2506.,0.01)], # array of positions (tuples) to record information
	source_order = 2, # maximun order
	order	 	 = 1, # minimun order
	freq		 = 7.5, # frequency
	p_velocity	 = 4000., # p-wave velocity
	s_velocity	 = 2310., # p-wave velocity
	density		 = 2000.,
	duration	 = 1.6, # simulation time
	pml_layer	 = 300., # PML length from border
	pml_coef 	 = 0.00001,
	pixel_size	 = 10., # movie resolution
	tpf	 		 = 0.01, # sampling time between frames
	source_dt	 = 0.001, # Source file sampling
	stress		 = M, # initial stress tensor: [Mxx, Myy, Mxy]
)
'''
	relations
	p-wave and density
	rho = a*vp^m, a = 0.31 and m = 0.25 rho in [g/cm^3] and vp in [m/s]
'''
#1580000
