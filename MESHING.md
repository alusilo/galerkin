# Meshing Guide

This document describes how meshes are used in the 2D solver and how to create or obtain them.

## Overview

The 2D Discontinuous Galerkin solver needs:

1. **A mesh file** — defines the triangulation (vertices and elements).
2. **A parameter file** (`.param`) — one line per element with physical properties: `rho vs vp`.

Both must exist before running the solver. Paths are relative to the project root (e.g. `resources/mesh/...`).

---

## Supported mesh formats

The code detects the format from the file extension and reads it accordingly.

| Format   | Extension(s)     | Description                          |
|----------|------------------|--------------------------------------|
| Triangle | `.ele` + `.node` | Two files: elements and nodes.       |
| GMSH     | `.msh`           | Single file (version 2.2 ASCII).     |
| Gambit   | `.neu`           | Gambit Neutral File (e.g. from Fluent). |

- **Triangle**: you provide the path to the `.ele` file; the reader looks for a `.node` file in the same directory with the same base name (e.g. `tlr.ele` → `tlr.node`).
- **GMSH**: only **triangular** elements (type 2) are used; other element types are ignored.
- **Gambit**: read via `Mesh2D.__neu()`.

---

## Where to put mesh files

Place mesh (and param) files under:

```
resources/mesh/
├── tri/           # Triangle format: tlr.ele, tlr.node, tlr.param (generated)
├── msh/           # GMSH: e.g. small.msh, mesh.msh
├── neu/           # Gambit: e.g. Maxwell05.neu
├── small.msh      # example included
├── small.param    # generated for small.msh (see below)
└── ...
```

In `main.py` you set `MESH_FILE` and `PARAM_FILE` to point to these paths (relative to project root or built from `_RESOURCES`).

---

## Important: the param script does not create meshes

**`scripts/generate_param_from_mesh.py`** only creates the **`.param` file** from a mesh you already have. It does **not** generate the mesh itself. To get a bigger (or any) mesh you must create it first with GMSH, Triangle, or trimeshGen; then run the script on the resulting `.msh` or `.ele` to get the param file.

---

## Parameter file (`.param`)

The solver expects a **text file** with **one line per triangular element**. Each line has three numbers (space-separated):

```
rho  vs  vp
```

- **rho** — density (e.g. kg/m³)
- **vs**  — S-wave velocity (e.g. m/s)
- **vp**  — P-wave velocity (e.g. m/s)

The number of lines must equal the number of triangular elements in the mesh. Order of lines must match the order of elements in the mesh file.

### Generating a param file from a GMSH mesh

For a **constant** material (same rho, vs, vp everywhere), you can generate the param file from a `.msh`:

```bash
uv run python scripts/generate_param_from_mesh.py resources/mesh/small.msh -o resources/mesh/small.param
```

Optional arguments:

- `--rho`, `--vs`, `--vp` — override default values (defaults: 2670, 1540, 3464).

Example with custom values:

```bash
uv run python scripts/generate_param_from_mesh.py resources/mesh/small.msh -o resources/mesh/small.param --rho 2500 --vp 3000 --vs 1500
```

For **heterogeneous** models you must produce the param file by other means (e.g. one value per element from a velocity model).

The script accepts **GMSH (`.msh`)** and **Triangle (`.ele`)** input:

```bash
uv run python scripts/generate_param_from_mesh.py resources/mesh/bigger.msh -o resources/mesh/bigger.param
uv run python scripts/generate_param_from_mesh.py resources/mesh/tri/mymesh.ele -o resources/mesh/tri/mymesh.param
```

---

## Generating a bigger mesh

To get a **bigger** (or finer) mesh you need to **create the mesh** first, then run `generate_param_from_mesh.py` on it.

### Bigger mesh with GMSH

1. Install [Gmsh](https://gmsh.info/) and open or create a `.geo` file.
2. **Increase the domain size** by changing the geometry (e.g. in `small.geo`, set `x_max`, `y_max` to larger values than 2).
3. **Control element size** with the characteristic length `len` (smaller `len` → more elements). For example in `small.geo`, `len = vs/(3*freq)`; decrease `len` for a finer mesh.
4. In GMSH: **Mesh → 2D** to generate the triangulation, then **File → Export → .msh** (ASCII, version 2).
5. Save the mesh as e.g. `resources/mesh/bigger.msh`.
6. Generate the param file:
   ```bash
   uv run python scripts/generate_param_from_mesh.py resources/mesh/bigger.msh -o resources/mesh/bigger.param
   ```
7. In `main.py` set `MESH_FILE` and `PARAM_FILE` to `bigger.msh` and `bigger.param`, and set **source and receivers inside the new domain** (e.g. if the domain is 0–10, use coordinates in that range).

### Bigger mesh with Triangle

1. Install [Triangle](https://www.cs.cmu.edu/~quake/triangle.html).
2. Create a `.poly` file describing a **larger** domain (vertices and segments). You can use a text editor or a script; the format is documented on the Triangle site.
3. Run Triangle, e.g.:
   ```bash
   triangle -pDq30 mydomain.poly
   ```
   This produces `mydomain.1.ele` and `mydomain.1.node`. Copy them to e.g. `resources/mesh/tri/` and rename if desired.
4. Generate the param file:
   ```bash
   uv run python scripts/generate_param_from_mesh.py resources/mesh/tri/mydomain.1.ele -o resources/mesh/tri/mydomain.param
   ```
5. In `main.py` point `MESH_FILE` to the `.ele` file and `PARAM_FILE` to the generated `.param`.

---

## Creating meshes

### Option 1: GMSH (recommended for getting started)

1. Install [Gmsh](https://gmsh.info/).
2. Create a geometry (`.geo`) and generate a 2D triangular mesh; export as `.msh` (ASCII, version 2).
3. Put the `.msh` in `resources/mesh/` (or `resources/mesh/msh/`).
4. Generate a param file with `scripts/generate_param_from_mesh.py` (see above).

The repository includes `resources/mesh/small.msh` (and `small.geo`) as an example.

### Option 2: Triangle program (`.ele` / `.node`)

1. Install [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) by Jonathan Shewchuk.
2. Create a `.poly` file describing the domain (vertices, segments, holes if needed).
3. Run Triangle to generate `.ele` and `.node`:

   ```bash
   triangle -pDq30 your_domain.poly
   ```

   This produces `your_domain.1.ele` and `your_domain.1.node` (or similar, depending on options). Rename or copy to the names expected by the solver (e.g. `tlr.ele`, `tlr.node`).
4. Create a `.param` file with one line per triangle (same format as above).

The solver expects the **path to the `.ele` file**; it will look for the `.node` file in the same directory with the same base name.

### Option 3: trimeshGen (project-specific, layered models)

Under `resources/mesh/tri/trimeshGen/` there is a mesh generator that:

- Reads a **JSON** description of layered models (e.g. `twolayers.json`, `onelayer.json`).
- Uses **RSF** (CWP) and **Triangle** to build a refined triangular mesh and velocity models.
- Outputs `.poly`, then `.ele` and `.node` (and related files) after several refinement iterations.

Requirements:

- [SCons](https://scons.org/)
- [CWP/RSF](https://wiki.seismic-unix.org/) (e.g. `CWPROOT` set).
- Triangle installed.

Typical workflow (from `trimeshGen/`):

1. Edit the JSON (e.g. `twolayers.json`) to define layers and velocities.
2. Run SCons to generate the mesh (and optional vp/vs/rho outputs).
3. Copy the resulting `.ele` and `.node` (e.g. `geo_model.30.ele`, `geo_model.30.node` if `geometry='geo_model.poly'` and `mesh_iter=30`) to `resources/mesh/tri/` and rename if desired (e.g. to `tlr.ele`, `tlr.node`).
4. Create or generate a `.param` file with one line per element (trimeshGen may produce velocity files you can convert to param).

Clean build artifacts (from `trimeshGen/`):

```bash
scons -c
# or
python SConstruct -c
```

(Exact command may depend on how SCons and the project are set up.)

---

## Summary checklist

1. **Obtain or generate a mesh** (GMSH `.msh`, Triangle `.ele`/`.node`, or Gambit `.neu`).
2. **Place it** under `resources/mesh/` (or a subfolder).
3. **Create a `.param` file** with one line per triangle: `rho vs vp`. Use `scripts/generate_param_from_mesh.py` for constant material from a `.msh`.
4. **Point the solver** to the mesh and param in `main.py` (`MESH_FILE`, `PARAM_FILE`).
5. **Set source and receivers** inside the mesh domain (see README).

For a minimal test, use `resources/mesh/small.msh`, generate `small.param` with the script, and set `MESH_FILE` and `PARAM_FILE` to those paths in `main.py`.
