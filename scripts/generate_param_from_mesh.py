#!/usr/bin/env python3
"""
Generate a default .param file from an existing mesh so the 2D solver can run.

This script does NOT create the mesh. It only writes a .param file (one line per
triangle: "rho vs vp") from a mesh you already have. To get a bigger mesh, create
it first with GMSH or Triangle (see MESHING.md), then run this script on the
resulting .msh or .ele file.

Supported input: .msh (GMSH), .ele (Triangle).

Usage (from project root):
  uv run python scripts/generate_param_from_mesh.py resources/mesh/small.msh
  uv run python scripts/generate_param_from_mesh.py resources/mesh/small.msh -o resources/mesh/small.param
  uv run python scripts/generate_param_from_mesh.py resources/mesh/tri/mymesh.ele
"""
import argparse
import logging
import os

logger = logging.getLogger(__name__)


def count_triangles_msh(path: str) -> int:
    with open(path, "rt") as f:
        data = f.read().splitlines()
    element_id = data.index("$Elements") + 1
    K = int(data[element_id])
    count = 0
    for a in data[element_id + 1 : element_id + K + 1]:
        parts = list(map(int, a.split(" ")[1:]))
        if parts[0] == 2:  # type 2 = 3-node triangle
            count += 1
    return count


def count_triangles_ele(path: str) -> int:
    """Triangle format: first line is 'K 3 0' or similar; then K lines of elements."""
    with open(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            parts = line.split()
            if parts:
                return int(parts[0])
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate default .param from an existing mesh (.msh or .ele). Does not create meshes."
    )
    parser.add_argument("mesh_file", help="Path to .msh (GMSH) or .ele (Triangle) mesh file")
    parser.add_argument("-o", "--output", help="Output .param path (default: mesh basename + .param)")
    parser.add_argument("--rho", type=float, default=2670.0, help="Density (default 2670)")
    parser.add_argument("--vs", type=float, default=1540.0, help="S-wave velocity (default 1540)")
    parser.add_argument("--vp", type=float, default=3464.0, help="P-wave velocity (default 3464)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    if not os.path.exists(args.mesh_file):
        parser.error(f"Mesh file not found: {args.mesh_file}")

    ext = os.path.splitext(args.mesh_file)[1].lower()
    if ext == ".msh":
        ntri = count_triangles_msh(args.mesh_file)
    elif ext == ".ele":
        ntri = count_triangles_ele(args.mesh_file)
    else:
        parser.error("Mesh file must be .msh (GMSH) or .ele (Triangle). See MESHING.md to create a mesh.")

    out = args.output or os.path.splitext(args.mesh_file)[0] + ".param"
    with open(out, "w") as f:
        for _ in range(ntri):
            f.write(f"{args.rho} {args.vs} {args.vp}\n")
    logger.info("Wrote %s lines to %s", ntri, out)


if __name__ == "__main__":
    main()
