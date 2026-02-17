import numpy as np

# m/s = kg.m^-2.s^-2.(m^3.kg^-1.s)


def source(t):
    fc = 14.5  # Hz
    a1 = -2000.0  # kg.m^-2.s^-2
    a2 = -((np.pi * fc) ** 2)  # Hz^2
    td = 0.08  # sec
    rho = 2200.0  # Kg.m^-3
    return a1 * (0.5 + a2 * (t - td) ** 2) * np.exp(a2 * (t - td) ** 2) / rho


def Dsource(t):
    fc = 14.5  # Hz
    a1 = -2000.0  # kg.m^-2.s^-2
    a2 = -((np.pi * fc) ** 2)  # Hz^2
    td = 0.08  # sec
    rho = 2200.0  # Kg.m^-3
    return 2 * a1 * a2 * (t - td) * (1.5 + a2 * (t - td) ** 2) * np.exp(a2 * (t - td) ** 2) / rho


dt = 0.0001
time = np.arange(0.0, 1.0, dt)
with open("../resources/source/wavelet_lamb.src", "w") as f:
    for t in time:
        f.write(str(Dsource(t)))
        f.write("\n")
