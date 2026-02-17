#!/usr/bin/env python3
"""Polynomial basis demo and plots (run from project root: uv run python scripts/polynomials.py)."""
import numpy as np
from scipy import linalg as la, special as sf
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

np.set_printoptions(precision=4, suppress=True)


def xytors(x, y):
    L1 = (np.sqrt(3) * y + 1) / 3
    L2 = (-3 * x - np.sqrt(3) * y + 2) / 6
    L3 = (3 * x - np.sqrt(3) * y + 2) / 6
    r = -L2 + L3 - L1
    s = -L2 - L3 + L1
    return r, s


def rstoab(r, s):
    Np = len(r)
    a = np.zeros(Np)
    for i in range(Np):
        if s[i] != 1:
            a[i] = 2 * (1 + r[i]) / (1 - s[i]) - 1
        else:
            a[i] = -1
    b = s
    return a, b


def JGLPoints(alpha, beta, N):
    x = [] if (N == 1) else sf.j_roots(N - 1, alpha, beta)[0]
    x = np.append(np.append([-1.0], x), [1.0])
    return x


def JPNormalized(N, alpha, beta, r):
    gamma0 = (2.0 ** (alpha + beta + 1) / (2 * N + alpha + beta + 1)) * (
        sf.gamma(N + alpha + 1)
        * sf.gamma(N + beta + 1)
        / (sf.gamma(N + alpha + beta + 1) * sf.factorial(N))
    )
    norm = 1.0 / np.sqrt(gamma0)
    x = norm * sf.eval_jacobi(N, alpha, beta, r)
    return x


def GradJacobiP(r, alpha, beta, N):
    dP = np.zeros(len(r))
    if N == 0:
        return dP
    return np.sqrt(N * (N + alpha + beta + 1)) * JPNormalized(N - 1, alpha + 1, beta + 1, r)


def Vandermonde1D(N, r):
    V = np.zeros((len(r), N + 1))
    for i in range(N + 1):
        V[:, i] = JPNormalized(i, 0, 0, r)
    return V


def Vandermonde1D2(N, r):
    V = np.zeros((len(r), N + 1))
    for i in range(N + 1):
        V[:, i] = r**i
    return V


def Vandermonde2D(N, r, s):
    V2D = np.zeros((len(r), int((N + 1) * (N + 2) / 2)))
    a, b = rstoab(r, s)
    k = 0
    for i in range(N + 1):
        for j in range(N - i + 1):
            V2D[:, k] = Simplex2DP(a, b, i, j)
            k = k + 1
    return V2D


def Vandermonde2D2(N, r, s):
    V2D = np.zeros((len(r), int((N + 1) * (N + 2) / 2)))
    k = 0
    for i in range(N + 1):
        for j in range(N - i + 1):
            V2D[:, k] = (r**i) * (s**j)
            k = k + 1
    return V2D


def GradVandermonde2D(N, r, s):
    V2Dr = np.zeros((len(r), int((N + 1) * (N + 2) / 2)))
    V2Ds = np.zeros((len(r), int((N + 1) * (N + 2) / 2)))
    a, b = rstoab(r, s)
    k = 0
    for i in range(N + 1):
        for j in range(N - i + 1):
            V2Dr[:, k], V2Ds[:, k] = GradSimplex2DP(a, b, i, j)
            k = k + 1
    return V2Dr, V2Ds


def Simplex2DP(a, b, i, j):
    h1 = JPNormalized(i, 0, 0, a)
    h2 = JPNormalized(j, 2 * i + 1, 0, b)
    return np.sqrt(2.0) * h1 * h2 * (1 - b) ** i


def GradSimplex2DP(a, b, idx, jdx):
    fa = JPNormalized(idx, 0, 0, a)
    dfa = GradJacobiP(a, 0, 0, idx)
    gb = JPNormalized(jdx, 2 * idx + 1, 0, b)
    dgb = GradJacobiP(b, 2 * idx + 1, 0, jdx)
    dmodedr = dfa * gb
    dmodeds = dfa * (gb * (0.5 * (1 + a)))
    tmp = dgb * ((0.5 * (1 - b)) ** idx)
    if idx > 0:
        dmodedr = dmodedr * ((0.5 * (1 - b)) ** (idx - 1))
        dmodeds = dmodeds * ((0.5 * (1 - b)) ** (idx - 1))
        tmp = tmp - 0.5 * idx * gb * ((0.5 * (1 - b)) ** (idx - 1))
    dmodedr = 2 ** (idx + 0.5) * dmodedr
    dmodeds = 2 ** (idx + 0.5) * (dmodeds + fa * tmp)
    return dmodedr, dmodeds


def Warpfactor(N, r):
    LGLr = JGLPoints(1, 1, N)
    req = np.linspace(-1, 1, N + 1)
    Veq = Vandermonde1D(N, req)
    Pmat = np.zeros((N + 1, len(r)))
    for i in range(N + 1):
        Pmat[i, :] = JPNormalized(i, 0, 0, r)
    Lmat = la.solve(Veq.T, Pmat)
    warp = np.dot(Lmat.T, (LGLr - req))
    zerof = 1 * (abs(r) < 1.0 - 1.0e-10)
    sf_val = 1.0 - (zerof * r) ** 2
    warp = warp * (1 / sf_val) + warp * (zerof - 1)
    return warp


class Poly(object):
    def __init__(self, **kwargs):
        super(Poly, self).__init__()
        self.N = kwargs["order"]
        self.Nfp = self.N + 1
        self.Np = int(self.Nfp * (self.Nfp + 1) / 2)
        self.Nfaces = 3
        self.NODETOL = 1e-12
        i = kwargs["i"]
        j = kwargs["j"]
        res = N
        self.x, self.y = self.Nodes2D(res)
        self.r, self.s = xytors(self.x, self.y)
        a, b = rstoab(self.r, self.s)
        self.m = int(j + (N + 1) * i + 1 - i * (i - 1) / 2)
        self.z = Simplex2DP(a, b, i, j)

    def Nodes2D(self, N):
        alpopt = [0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999, 1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258]
        alpha = alpopt[N - 1] if N < 16 else 5 / 3
        Nfp = N + 1
        Np = int(Nfp * (Nfp + 1) / 2)
        L1 = np.zeros(Np)
        L2 = np.zeros(Np)
        L3 = np.zeros(Np)
        i = 0
        for n in range(N + 1):
            for m in range(N - n + 1):
                L1[i] = n / N
                L3[i] = m / N
                i = i + 1
        L2 = 1.0 - L1 - L3
        x = -L2 + L3
        y = (-L2 - L3 + 2 * L1) / np.sqrt(3)
        blend1 = 4 * L2 * L3
        blend2 = 4 * L1 * L3
        blend3 = 4 * L1 * L2
        warpf1 = Warpfactor(N, L3 - L2)
        warpf2 = Warpfactor(N, L1 - L3)
        warpf3 = Warpfactor(N, L2 - L1)
        warp1 = blend1 * warpf1 * (1 + (alpha * L1) ** 2)
        warp2 = blend2 * warpf2 * (1 + (alpha * L2) ** 2)
        warp3 = blend3 * warpf3 * (1 + (alpha * L3) ** 2)
        x = x + 1 * warp1 + np.cos(2 * np.pi / 3) * warp2 + np.cos(4 * np.pi / 3) * warp3
        y = y + 0 * warp1 + np.sin(2 * np.pi / 3) * warp2 + np.sin(4 * np.pi / 3) * warp3
        return x, y


N = 2
i = 2
j = 0

font = {"family": "serif", "style": "italic", "weight": "normal", "size": 24}
plt.rc("font", **font)
plt.rc("text", usetex=True)

x = np.linspace(-1, 1, N)
y = np.linspace(-1, 1, N)
obj = Poly(order=N, i=i, j=j)
print(obj.r)
V1 = Vandermonde2D2(N, obj.r, obj.s)
V2 = Vandermonde2D(N, obj.r, obj.s)
M1 = la.inv(np.dot(V1, V1.T))
M2 = la.inv(np.dot(V2, V2.T))
print(M1)
print(M2)

if i >= 0 and j >= 0 and i + j <= N:
    obj = Poly(order=N, i=i, j=j)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_trisurf(obj.r, obj.s, obj.z, color="#555555", linewidth=0, antialiased=False)
    ax._axis3don = False
    ax.set_aspect("equal")
    ax.set_title(r"$\displaystyle m = ${}".format(obj.m))
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter("%.02f"))
    ax.view_init(20, -130)
    plt.show()
