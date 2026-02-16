# Galerkin - Discontinuous Galerkin Method for Wave Propagation

A Python implementation of the **Discontinuous Galerkin (DG) Method** for solving wave propagation problems in 2D using the nodal DG approach. This code is based on the MATLAB implementation by **Jan S. Hesthaven** and **Tim Warburton** from their book *Nodal Discontinuous Galerkin Methods*.

## Features

- **2D Wave Propagation**: Full-featured 2D implementation with support for:
  - Triangular meshes
  - PML (Perfectly Matched Layer) boundary conditions
  - Custom source functions
  - Seismogram recording
  - Field visualization

- **1D Examples**: Simple 1D examples demonstrating basic DG concepts

## Project Structure

```
galerkin/
├── src/
│   └── galerkin/          # Main 2D DG implementation
│       ├── __init__.py
│       ├── wave2D.py      # Main wave propagation driver
│       ├── nodalDG.py     # Nodal DG base class
│       ├── Mesh2D.py      # 2D mesh handling
│       ├── polynomials.py # Polynomial basis functions
│       └── ...            # Additional utilities
├── examples/
│   └── 1d/                # 1D examples (educational)
│       ├── Advec1D.py
│       ├── AdvecDriver1D.py
│       └── ...
├── resources/             # Mesh files, sources, and output
│   ├── mesh/              # Mesh files (.ele, .neu, .msh)
│   ├── source/            # Source function definitions
│   └── output/            # Simulation outputs (gitignored)
├── pyproject.toml         # UV package configuration
└── README.md
```

## Installation

This project uses [uv](https://github.com/astral-sh/uv) for fast Python package management.

### Prerequisites

- Python 3.8 or higher
- [uv](https://github.com/astral-sh/uv) package manager

Install uv:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Setup

1. Clone the repository:
```bash
git clone <repository-url>
cd galerkin
```

2. Install dependencies using uv:
```bash
uv sync
```

This will create a virtual environment and install all dependencies.

3. Activate the virtual environment:
```bash
source .venv/bin/activate  # On Linux/Mac
# or
.venv\Scripts\activate     # On Windows
```

## Usage

### Running 2D Simulations

The main 2D code is located in `src/galerkin/`. You can run simulations using:

```python
from galerkin import WaveDrive2D

# Create a wave propagation simulation
obj = WaveDrive2D(
    project='my_simulation',
    mesh_file='resources/mesh/tri/your_mesh.ele',
    param_file='resources/mesh/tri/your_mesh.param',
    src_position=(2146.76, 383.27),
    src_smooth=50.0,
    gather=[(3200., 0.01), (2500., 1100.)],
    source_order=3,
    order=2,
    src_freq=7.5,
    src_delay=1/7.5,
    duration=0.5,
    pml_layer=(300., 300., 0., 300.),
    pml_coef=0.0015,
    pixel_size=10.0,
    tpf=0.01,
    stress=(1, 1, 0),
    displacement=(0, 0),
)

# Run the simulation
obj.run()
```

Or use the example script:
```bash
python src/galerkin/main.py
```

### Running 1D Examples

The 1D examples are located in `examples/1d/` and demonstrate basic DG concepts:

```bash
cd examples/1d
python AdvecDriver1D.py
```

## Dependencies

- **numpy** >= 1.20.0 - Numerical computations
- **scipy** >= 1.7.0 - Scientific computing utilities
- **matplotlib** >= 3.4.0 - Visualization

### Optional Dependencies

- **Intel MKL** - For optimized BLAS operations (automatically used if available)
  - Can significantly improve performance for large simulations
  - Install via conda: `conda install mkl`

### Development Dependencies

- **pytest** - Testing framework
- **black** - Code formatting
- **ruff** - Fast Python linter

## Configuration

### Mesh Files

The code supports triangular meshes in various formats:
- `.ele` files (Triangle format)
- `.neu` files (Gambit format)
- `.msh` files (Gmsh format)

Place mesh files in `resources/mesh/`.

### Source Functions

Source functions are defined in `resources/source/` and can be customized for different wave types (Gaussian, Ricker, etc.).

### Output

Simulation outputs are saved to `resources/output/<project_name>/` and include:
- Field snapshots
- Seismograms
- Parameter files

## Development

### Code Formatting

Format code using black:
```bash
uv run black src/
```

### Linting

Check code quality with ruff:
```bash
uv run ruff check src/
```

### Running Tests

```bash
uv run pytest
```

## Mathematical Background

The Discontinuous Galerkin Method is a numerical technique for solving partial differential equations. This implementation solves the 2D elastic wave equation:

\[
\rho \frac{\partial^2 \mathbf{u}}{\partial t^2} = \nabla \cdot \boldsymbol{\sigma} + \mathbf{f}
\]

where:
- \(\rho\) is density
- \(\mathbf{u}\) is displacement
- \(\boldsymbol{\sigma}\) is stress tensor
- \(\mathbf{f}\) is source term

The method uses:
- Nodal basis functions on triangular elements
- Upwind flux for numerical stability
- Runge-Kutta time integration
- PML for absorbing boundary conditions

## References

- Hesthaven, J. S., & Warburton, T. (2008). *Nodal Discontinuous Galerkin Methods: Algorithms, Analysis, and Applications*. Springer.

## License

[Add your license information here]

## Contributing

[Add contribution guidelines here]

## Contact

[Add contact information here]
