#!/usr/bin/env python3
"""Generate wavelet and timing plots (run from project root: uv run python scripts/wavelet_generator.py)."""
from __future__ import unicode_literals

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

plt.rc("font", family="serif", serif="Times")
plt.rc("text", usetex=True)
plt.rc("xtick", labelsize=8)
plt.rc("ytick", labelsize=8)
plt.rc("axes", labelsize=8)
plt.rc("legend", fontsize=6)
plt.rc("axes", titlesize=6)

yFormatter = FormatStrFormatter("%.1f")
xIntFormatter = FormatStrFormatter("%d")

width = 3.487
height = width / 1.618

# data for fixed sigma = 10m
data = {"order": [2, 3, 4, 5], "time": [205.963, 281.244, 335.427, 475.755]}
fig, ax = plt.subplots()
fig.subplots_adjust(left=0.18, bottom=0.17, right=0.98, top=0.92)
ax.plot(data["order"], data["time"], "-", markersize=1, linestyle="--", c="black", linewidth=0.25)
ax.set_xlabel(r"orden de interpolaci\'on")
ax.set_ylabel("tiempo ($s$)")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_major_formatter(yFormatter)
ax.xaxis.set_major_formatter(xIntFormatter)
fig.set_size_inches(width, height)
fig.savefig("order_vs_time.pdf")

# data for fixed order = 2
data1 = {"sigma": [10, 20, 30, 40], "time": [205.963, 218.397, 206.763, 207.900]}
data2 = {"sigma": [10, 20, 30, 40], "time": [293.768, 316.34, 319.279, 315.212]}
fig, ax = plt.subplots()
fig.subplots_adjust(left=0.18, bottom=0.17, right=0.98, top=0.92)
ax.plot(data1["sigma"], data1["time"], "-", markersize=1, linestyle="--", c="black", linewidth=0.25, label="2do orden")
ax.plot(data2["sigma"], data2["time"], "-", markersize=1, linestyle="-.", c="black", linewidth=0.25, label="3er orden")
ax.set_xlabel("radio del soporte ($m$)")
ax.set_ylabel("tiempo ($s$)")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_major_formatter(yFormatter)
plt.legend(loc=0)
fig.set_size_inches(width, height)
fig.savefig("sigma_vs_time.pdf")

src_freq = 7.5
src_delay = 1 / src_freq


def createWavelet(time):
    dispv = 0
    half = 0.5
    two = 2.0
    tau = time - src_delay
    twopi = two * np.pi
    a1 = -2.0
    a2 = -half * half * (twopi * src_freq) ** 2
    temp = a2 * tau * tau
    if dispv == 0:
        wout = np.exp(temp)
    else:
        a1 = two * a1
        wout = a1 * (tau) * np.exp(temp)
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.12, bottom=0.17, right=0.98, top=0.92)
    plt.xlabel(r"tiempo ($s$)")
    plt.ylabel("Amplitud")
    plt.plot(time, wout, color="black", linewidth=1)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_major_formatter(yFormatter)
    fig.set_size_inches(width, height)
    fig.savefig("wavelet.pdf")


time = np.linspace(0, 2.0, 1000)
createWavelet(time)
