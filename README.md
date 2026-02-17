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
│       └── ...            # Additional utilities
├── scripts/               # Complement scripts (plotting, demos)
│   ├── visualize.py      # Seismograms and wave-field animation
│   ├── generate_param_from_mesh.py
│   ├── plot_src.py       # Initial source distribution
│   ├── plot_src_support.py
│   ├── plot_field.py     # Wave field animation
│   ├── wavelet_generator.py
│   ├── polynomials.py    # Polynomial basis demo
│   └── quadrature_int.py
├── examples/
│   └── 1d/                # 1D examples (educational)
├── resources/             # Mesh files, sources, and output
│   ├── mesh/              # Mesh files (.ele, .neu, .msh)
│   ├── source/            # Source function definitions
│   └── output/            # Simulation outputs (gitignored)
├── pyproject.toml         # UV package configuration
├── run.sh                 # Run 2D simulation (from project root)
└── README.md
```

## Installation

This project uses [uv](https://github.com/astral-sh/uv) for fast Python package management.

### Prerequisites

- Python 3.10 or higher
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

Or run the simulation from the project root using the provided script:
```bash
./run.sh
```
(This runs `uv run python src/galerkin/main.py`.)

Alternatively:
```bash
uv run python src/galerkin/main.py
```

### Visualizing results

After a 2D run, wave field and seismograms are under `resources/output/<project_name>/`. To plot them:

```bash
uv run python scripts/visualize.py
uv run python scripts/visualize.py --project erase
```

This opens:
- **Seismograms**: Vx and Vy at each receiver.
- **Wave field**: animated snapshot of Vx (click to run).

Options:
- `--no-movie`: only plot seismograms.
- `--snapshot N`: show a single wave-field frame instead of animation.

Other complement scripts (run from project root):
- `uv run python scripts/plot_field.py` — wave field animation (reads `model.param` and mesh from the project output folder).
- `uv run python scripts/plot_src.py` — plot initial source distribution.
- `uv run python scripts/plot_src_support.py` — plot spatial source support.
- `uv run python scripts/polynomials.py` — polynomial basis demo.
- `uv run python scripts/quadrature_int.py` — quadrature/interpolation demo.
- `uv run python scripts/wavelet_generator.py` — generate wavelet and timing PDFs.

### Running 1D Examples

The 1D examples are located in `examples/1d/` and demonstrate basic DG concepts:

```bash
cd examples/1d
python AdvecDriver1D.py
```

## Dependencies

- **numpy** >= 1.24.0 - Numerical computations
- **scipy** >= 1.11.0 - Scientific computing utilities
- **matplotlib** >= 3.7.0 - Visualization

### Optional Dependencies

- **Intel MKL** - For optimized BLAS operations (automatically used if available)
  - Can significantly improve performance for large simulations
  - Install via conda: `conda install mkl`

### Development Dependencies

- **pytest** - Testing framework
- **black** - Code formatting
- **ruff** - Fast Python linter

## Configuration

### Resources and mesh (you need a mesh first)

The 2D solver expects **mesh and parameter files to already exist**. Paths in the code are relative to the **project root**, so run from the repo root (e.g. `uv run python src/galerkin/main.py`).

1. **Mesh files**  
   Place meshes under `resources/mesh/`. Supported formats:
   - **Triangle** (`.ele` + `.node`): typically generated first (see below).
   - **GMSH** (`.msh`): e.g. from [Gmsh](https://gmsh.info/) or the included `resources/mesh/small.msh`.
   - **Gambit** (`.neu`): e.g. `resources/mesh/neu/*.neu`.

2. **Triangle format (`.ele` / `.node`)**  
   The example in `main.py` uses Triangle output (e.g. `resources/mesh/tri/tlr.ele` and `tlr.node`). These are **not** included in the repo; you must generate them:
   - Use the **Triangle** program to produce `.ele` and `.node` from a `.poly` description, or
   - Use the mesh generator under `resources/mesh/tri/trimeshGen/` (SCons + C++/RSF), which writes `.ele`/`.node` (and related files) from the JSON layer models there.

3. **Parameter file (`.param`)**  
   The solver needs a **parameter file** (e.g. `tlr.param`) with one line per **element**:  
   `rho vs vp` (density, S-wave speed, P-wave speed), space-separated. The number of lines must equal the number of triangular elements in the mesh.

4. **If you see “file does not exist”**  
   - Ensure the mesh (and for Triangle, both `.ele` and `.node`) and the `.param` file exist at the paths used in `main.py` (under `resources/mesh/`).
   - Or point `mesh_file` and `param_file` in `main.py` to an existing mesh (e.g. a `.msh` and a matching `.param` with the correct number of lines).

5. **Quick run with the included GMSH mesh**  
   To run without generating a Triangle mesh, use the provided `small.msh` and generate a param file:
   ```bash
   uv run python scripts/generate_param_from_mesh.py resources/mesh/small.msh -o resources/mesh/small.param
   ```
   Then in `main.py` set `MESH_FILE` to `os.path.join(_RESOURCES, "small.msh")` and `PARAM_FILE` to `os.path.join(_RESOURCES, "small.param")`.

### Source Functions

Source functions are defined in `resources/source/` and can be customized for different wave types (Gaussian, Ricker, etc.).

### Output

Simulation outputs are saved to `resources/output/<project_name>/` and include:
- Field snapshots
- Seismograms
- Parameter files

## Development

**Use `uv run`, not `uvx`**, for project tools. From the project root (after `uv sync`):

- `uv run <script>` — runs using this project’s virtualenv and dependencies.
- `uvx <tool>` — runs a tool in a separate, temporary environment (for one-off tools not in this project). If `uvx` fails, use `uv run` instead for tools listed below.

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

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

The MIT License is permissive and allows:
- Academic and commercial use
- Modification and distribution
- Private use

The only requirement is that the original copyright notice and license are included in any copies or substantial portions of the software.

## Contributing

Anyone is free to fork this repository and create a new repository from it. You may use, modify, and distribute the code in accordance with the MIT License.

## Contact

For questions or feedback, you can reach the author at **silva.l.a.l@gmail.com**.
