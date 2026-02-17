import os

from galerkin.wave2D import WaveDrive2D

# Optional: Use Intel MKL for optimized BLAS if available
try:
    import mkl

    mkl.set_num_threads(4)
except ImportError:
    pass  # Continue without MKL optimizations

# Paths relative to project root (so they work regardless of cwd)
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_RESOURCES = os.path.join(_PROJECT_ROOT, "resources", "mesh")

PROJECT_NAME = "erase"

M = (1, 1, 0)
F = (0, 0)

# Mesh and param files: you must generate the mesh first (see README).
# Option A: Triangle format (.ele + .node) — generate with resources/mesh/tri/trimeshGen or Triangle.
# Option B: Use existing GMSH mesh, e.g. small.msh, and a matching .param file.
# Uncomment for Triangle mesh (tlr.ele / tlr.node); generate mesh first.
# MESH_FILE = os.path.join(_RESOURCES, "tri", "tlr.ele")
# PARAM_FILE = os.path.join(_RESOURCES, "tri", "tlr.param")

# small.msh: domain is [0, 2] x [0, 2]; generate param with scripts/generate_param_from_mesh.py
MESH_FILE = os.path.join(_RESOURCES, "mesh.msh")
PARAM_FILE = os.path.join(_RESOURCES, "mesh.param")

# Source and gather must lie inside the mesh domain.
# For small.msh (0 ≤ x,y ≤ 2) use coordinates in that range.
obj = WaveDrive2D(
    project=PROJECT_NAME,
    mesh_file=MESH_FILE,
    param_file=PARAM_FILE,
    src_position=(2000, 1000),  # inside [0,2] x [0,2] for small.msh
    src_smooth=10,  # spatial support (fraction of domain for small mesh)
    gather=[
        (2500, 100),
        (3000, 100),
        (3500, 100),
    ],  # receiver positions inside domain
    source_order=6,
    order=4,
    src_freq=7.5,
    src_delay=1 / 7.5,
    duration=1.5,
    pml_layer=(10, 10, 10, 10),  # PML thickness for [0,2] domain
    pml_coef=0.0015,
    pixel_size=0.5,  # resolution for 2x2 domain
    tpf=0.01,
    stress=M,
    displacement=F,
)
"""
	relations
	p-wave and density
	rho = a*vp^m, a = 0.31 and m = 0.25 rho in [g/cm^3] and vp in [m/s]
"""
# 1580000
