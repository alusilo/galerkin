#!/usr/bin/env python3
"""
Visualize wave field and seismograms from a 2D run.

Usage (from project root):
  uv run python scripts/visualize.py
  uv run python scripts/visualize.py --project erase
  uv run python scripts/visualize.py --project erase --no-movie   # seismograms only
"""
import argparse
import json
import logging
import os

import numpy as np
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize wave field and seismograms")
    parser.add_argument(
        "--project",
        default="erase",
        help="Project name (output subfolder under resources/output)",
    )
    parser.add_argument(
        "--no-movie",
        action="store_true",
        help="Skip wave field animation, only plot seismograms",
    )
    parser.add_argument(
        "--snapshot",
        type=int,
        default=None,
        metavar="N",
        help="If set, show only snapshot N instead of animation",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(root, "resources", "output", args.project)
    param_path = os.path.join(out_dir, "model.param")

    if not os.path.exists(param_path):
        logger.error("not found %s", param_path)
        logger.error("Run the 2D solver first (e.g. uv run python src/galerkin/main.py)")
        return

    with open(param_path) as f:
        param = json.load(f)

    xmin, xmax, ymin, ymax = param["limits"]
    duration = param["duration"]
    gather = param["gather"]
    source = param["source"]["position"]

    # Use km for display if domain is large
    scale = 1000.0 if max(xmax - xmin, ymax - ymin) > 100 else 1.0

    # ---- Seismograms ----
    traces_vx_path = os.path.join(out_dir, "tracesVx.npy")
    traces_vy_path = os.path.join(out_dir, "tracesVy.npy")
    if os.path.exists(traces_vx_path):
        data = np.load(traces_vx_path)
        time = data[0]
        traces_vx = data[1:]
        ntraces = traces_vx.shape[0]
        has_vy = os.path.exists(traces_vy_path)
        traces_vy = np.load(traces_vy_path)[1:] if has_vy else None
        nrows = ntraces * 2 if has_vy else ntraces
        fig, axes = plt.subplots(nrows, 1, sharex=True, figsize=(8, 1.5 * nrows))
        if nrows == 1:
            axes = [axes]
        for i in range(ntraces):
            ax = axes[2 * i] if has_vy else axes[i]
            ax.plot(time, traces_vx[i], "k-", lw=0.8)
            ax.set_ylabel("Vx")
            ax.set_title(f"Receiver at {gather[i]}" if i < len(gather) else f"Receiver {i}")
            ax.grid(True, alpha=0.3)
            if has_vy:
                axes[2 * i + 1].plot(time, traces_vy[i], "k-", lw=0.8)
                axes[2 * i + 1].set_ylabel("Vy")
                axes[2 * i + 1].grid(True, alpha=0.3)
        axes[-1].set_xlabel("Time (s)")
        fig.suptitle("Seismograms")
        plt.tight_layout()
        plt.show()
    else:
        logger.warning("Seismogram file not found: %s", traces_vx_path)

    # ---- Wave field (movie) ----
    if args.no_movie:
        return

    movie_vx_path = os.path.join(out_dir, "movieVx.npy")
    if not os.path.exists(movie_vx_path):
        logger.warning("Movie file not found: %s", movie_vx_path)
        return

    data = np.load(movie_vx_path)
    nframes = len(data)
    vmin, vmax = float(data.min()), float(data.max())
    dt = duration / max(nframes - 1, 1)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    extent = [xmin / scale, xmax / scale, ymax / scale, ymin / scale]
    im = ax.imshow(
        data[0].T,
        extent=extent,
        aspect="equal",
        cmap="RdBu_r",
        origin="upper",
        vmin=vmin,
        vmax=vmax,
    )
    for gth in gather:
        ax.scatter(gth[0] / scale, gth[1] / scale, c="green", marker="v", s=80, zorder=2)
    ax.scatter(source[0] / scale, source[1] / scale, c="yellow", marker="*", s=200, zorder=2)
    ax.set_xlabel("x (km)" if scale == 1000 else "x")
    ax.set_ylabel("z (km)" if scale == 1000 else "y")
    plt.colorbar(im, ax=ax, label="Vx")
    ax.set_title("Wave field Vx (t=0 s)")

    if args.snapshot is not None:
        idx = max(0, min(args.snapshot, nframes - 1))
        im.set_data(data[idx].T)
        ax.set_title(f"Wave field Vx (t={idx * dt:.3f} s)")
        plt.show()
        return

    def update(frame):
        im.set_data(data[frame].T)
        ax.set_title(f"Wave field Vx (t={frame * dt:.3f} s)")
        return (im,)

    from matplotlib import animation

    ani = animation.FuncAnimation(
        fig, update, frames=nframes, interval=50, repeat=True, blit=False
    )
    plt.show()


if __name__ == "__main__":
    main()
