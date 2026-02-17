#!/usr/bin/env python3
"""Quadrature and interpolation demo (run from project root: uv run python scripts/quadrature_int.py)."""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as la

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_src = os.path.join(_PROJECT_ROOT, "src")
if _src not in sys.path:
    sys.path.insert(0, _src)
from galerkin.cubatureData2D import cubatureData2D
from galerkin.nodalDG import GradVandermonde2D, Nodes2D, Vandermonde2D, xytors


def f(x, y):
    return np.exp(-2j * np.pi * (x + y) / 5).real


def area(vertices):
    (x1, y1), (x2, y2), (x3, y3) = vertices
    return np.abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2)


vertices = np.array([[1, 2], [5, 4], [2, 6]])
triangle = np.append(vertices, [vertices[0]], axis=0)

fig, ax = plt.subplots()
ax.plot(triangle[:, 0], triangle[:, 1], "b-")
for vx, vy in vertices:
    ax.scatter(vx, vy, c="b")
    ax.text(vx, vy, " P({},{})".format(vx, vy), fontsize=9)
ax.set_title(r"Elemento sobre sistema global ($x$, $y$)")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_aspect("equal")
ax.grid(True)
plt.tight_layout()
plt.show()

N = 4
Np = int((N + 1) * (N + 2) / 2)
x, y = Nodes2D(N)
v_equi = np.array([[-1, np.min(y)], [1, np.min(y)], [0, np.max(y)]])
equi_tri = np.append(v_equi, [v_equi[0]], axis=0)

fig, ax = plt.subplots()
ax.plot(equi_tri[:, 0], equi_tri[:, 1], "b-")
ax.scatter(x, y, c="b", s=20)
ax.set_title(r"Puntos de soporte deformados sobre triangulo isoseles, sistema ($x'$, $y'$)")
ax.set_xlabel(r"$x'$")
ax.set_ylabel(r"$y'$")
ax.set_aspect("equal")
ax.grid(True)
plt.tight_layout()
plt.show()

r, s = xytors(x, y)
v_ref = np.array([[-1, -1], [1, -1], [-1, 1]])
(x1, y1), (x2, y2), (x3, y3) = vertices
J = np.abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 4)
ref_tri = np.append(v_ref, [v_ref[0]], axis=0)

fig, ax = plt.subplots()
ax.plot(ref_tri[:, 0], ref_tri[:, 1], "b-")
ax.scatter(r, s, c="b", s=20)
ax.set_title(r"Elemento y puntos de soporte sobre sistema estandar ($r$, $s$)")
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$s$")
ax.set_aspect("equal")
ax.grid(True)
plt.tight_layout()
plt.show()

# Map (r,s) reference to global (x,y); same as nodalDG: 0.5*(-va*(r+s) + vb*(1+r) + vc*(1+s))
v0, v1, v2 = vertices[0], vertices[1], vertices[2]
x = 0.5 * (-v0[0] * (r + s) + v1[0] * (1 + r) + v2[0] * (1 + s))
y = 0.5 * (-v0[1] * (r + s) + v1[1] * (1 + r) + v2[1] * (1 + s))

fig, ax = plt.subplots()
ax.plot(triangle[:, 0], triangle[:, 1], "b-")
ax.scatter(x, y, c="b", s=20)
ax.set_title(r"Puntos de soporte sobre sistema global ($x$, $y$)")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_aspect("equal")
ax.grid(True)
plt.tight_layout()
plt.show()

V = Vandermonde2D(N, r, s)
Vr, Vs = GradVandermonde2D(N, r, s)
grade = 0
function = f(x, y)
P_n = V[:, grade]
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection="3d")
ax.scatter(x, y, function, s=25, c="r")
ax.plot_trisurf(x, y, function, antialiased=True)
ax.view_init(50, -210)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("amplitude")

invV = la.inv(V)
samples = 100
pr = -1.0 * np.ones(samples)
ps = np.linspace(-1, 1, samples)
px = 0.5 * (-v0[0] * (pr + ps) + v1[0] * (1 + pr) + v2[0] * (1 + ps))
py = 0.5 * (-v0[1] * (pr + ps) + v1[1] * (1 + pr) + v2[1] * (1 + ps))
phi_n = Vandermonde2D(N, pr, ps)
interp = np.dot(phi_n, invV)
pz = np.dot(interp, function)
plt.show()

A = area(vertices)
print("area={}, 2J = {}".format(A, 2 * J))

cub = cubatureData2D(N=N)
cubx = 0.5 * (-v0[0] * (cub.r + cub.s) + v1[0] * (1 + cub.r) + v2[0] * (1 + cub.s))
cuby = 0.5 * (-v0[1] * (cub.r + cub.s) + v1[1] * (1 + cub.r) + v2[1] * (1 + cub.s))
integral = J * np.dot(cub.w, f(cubx, cuby))
print(integral)
