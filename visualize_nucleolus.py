#!/usr/bin/env python3
"""
visualize_nucleolus.py

Reads a nucleolus trajectory file produced by run_nucleolus and generates
an animated matplotlib figure showing particle positions coloured by polymer
type, with a background gradient overlay illustrating the chemical gradient.

Usage
-----
    python visualize_nucleolus.py [OPTIONS] <traj_file>

Options
    --output FILE       Save animation to FILE (mp4/gif).  Default: display.
    --fps N             Frames per second (default 5).
    --gradient-length L Draw gradient from x=0 to x=L (default: auto).
    --skip N            Only render every N-th frame (default 1).
    --title TEXT        Figure title.

Trajectory format (one line per particle per frame)
----------------------------------------------------
    FRAME <n>
    <id> <poly_type> <x> <y> <copy>
    ...

(This is the format written by writeFrame() in run_nucleolus.cpp.)
"""

import argparse
import sys
from collections import defaultdict

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import numpy as np


# --------------------------------------------------------------------------- #
# Colour palette: 4 polymer types × 4 copies → 16 distinct colours.
# We use hue = polymer type, saturation/value vary by copy index.
# --------------------------------------------------------------------------- #
POLYMER_COLORS = [
    ["#e6194b", "#f58231", "#ffe119", "#bfef45"],   # type 0 – reds/oranges
    ["#3cb44b", "#42d4f4", "#4363d8", "#911eb4"],   # type 1 – greens/blues
    ["#f032e6", "#a9a9a9", "#fabed4", "#aaffc3"],   # type 2 – pinks/greys
    ["#ffd8b1", "#dcbeff", "#469990", "#000075"],   # type 3 – pastels/darks
]

BACKBONE_COLOR = "#222222"
PARTICLE_RADIUS = 0.4
BOND_LW = 1.5


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

def parse_traj(path):
    """Return list of frames.  Each frame is a list of (id, poly, x, y, copy)."""
    frames = []
    current = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("FRAME"):
                if current:
                    frames.append(current)
                current = []
            else:
                parts = line.split()
                if len(parts) >= 5:
                    pid, poly, x, y, copy_ = (
                        int(parts[0]), int(parts[1]),
                        float(parts[2]), float(parts[3]),
                        int(parts[4]),
                    )
                    current.append((pid, poly, x, y, copy_))
    if current:
        frames.append(current)
    return frames


def build_bonds(particles, n0=16, n_seg=4):
    """
    Infer backbone bonds from particle positions: connect consecutive ids
    within the same polymer (same copy, same polymer type) if distance ≤ √2.
    Returns list of ((x1,y1),(x2,y2)) pairs.
    """
    by_id = {p[0]: p for p in particles}
    bonds = []
    for p in particles:
        pid, poly, x, y, copy_ = p
        next_id = pid + 1
        # Check that next_id is in same polymer (same poly type, same copy)
        if next_id in by_id:
            np2 = by_id[next_id]
            if np2[2] == copy_:           # same copy
                dx = abs(np2[2+1] - x)    # x distance
                dy = abs(np2[2+2] - y)    # y distance
                # Allow wrapping in y isn't drawn here (keep simple)
                dist = (np2[3] - y)**2 + (np2[2] - x)**2  # note: np2[2]=x, np2[3]=y
                if dist < 2.5:
                    bonds.append(((x, y), (np2[2], np2[3])))
    return bonds


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description="Animate nucleolus trajectory.")
    ap.add_argument("traj", help="Trajectory file (_traj.txt)")
    ap.add_argument("--output", default=None, help="Save to file (mp4/gif)")
    ap.add_argument("--fps", type=int, default=5)
    ap.add_argument("--gradient-length", type=float, default=None,
                    help="Column length L for gradient overlay")
    ap.add_argument("--skip", type=int, default=1,
                    help="Render every N-th frame")
    ap.add_argument("--title", default="Nucleolus simulation")
    ap.add_argument("--width", type=float, default=None,
                    help="Column width W (y-axis limit)")
    args = ap.parse_args()

    print(f"Reading {args.traj} ...", flush=True)
    all_frames = parse_traj(args.traj)
    if not all_frames:
        print("No frames found.", file=sys.stderr)
        sys.exit(1)
    frames = all_frames[::args.skip]
    print(f"  {len(all_frames)} frames total, rendering {len(frames)}.", flush=True)

    # Auto-detect axis limits from data
    all_x = [p[2] for fr in frames for p in fr]
    all_y = [p[3] for fr in frames for p in fr]
    x_min = min(all_x) - 1.0
    x_max = max(all_x) + 1.0
    y_min = min(all_y) - 1.0
    y_max = max(all_y) + 1.0

    if args.width is not None:
        y_min = -0.5
        y_max = args.width + 0.5

    gradient_L = args.gradient_length if args.gradient_length else (x_max * 0.6)

    # ---------------------------------------------------------------------- #
    # Figure setup
    # ---------------------------------------------------------------------- #
    fig, ax = plt.subplots(figsize=(max(8, (x_max - x_min) * 0.4),
                                    max(5, (y_max - y_min) * 0.6)))
    ax.set_aspect("equal")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("x (column axis)")
    ax.set_ylabel("y")
    ax.set_title(args.title)

    # Gradient background: column from x=0 to x=gradient_L
    gradient_img_width = 256
    xs = np.linspace(0, gradient_L, gradient_img_width)
    gamma_vals = np.clip(xs / gradient_L, 0, 1)
    gradient_data = np.tile(gamma_vals[np.newaxis, :], (10, 1))
    ax.imshow(
        gradient_data,
        extent=[0, gradient_L, y_min, y_max],
        origin="lower",
        aspect="auto",
        cmap="Blues",
        alpha=0.18,
        zorder=0,
    )
    ax.axvline(0, color="black", lw=1.5, zorder=1, label="hard wall")
    ax.axvline(gradient_L, color="steelblue", lw=1.5, linestyle="--",
               zorder=1, label=f"x = L = {gradient_L:.0f}")

    # Legend patches for polymer types
    legend_patches = []
    for t in range(4):
        patch = mpatches.Patch(color=POLYMER_COLORS[t][0], label=f"Polymer {t}")
        legend_patches.append(patch)
    legend_patches.append(
        mpatches.Patch(color="none", label="")
    )
    ax.legend(handles=legend_patches, loc="upper right", fontsize=7)

    # Particle scatter + bond lines — will be updated each frame
    scat = ax.scatter([], [], s=200, zorder=3, edgecolors="k", linewidths=0.5)
    bond_lines = []

    frame_text = ax.text(0.02, 0.97, "", transform=ax.transAxes,
                         fontsize=9, va="top", family="monospace")

    # ---------------------------------------------------------------------- #
    # Animation update function
    # ---------------------------------------------------------------------- #
    def update(frame_idx):
        nonlocal bond_lines
        particles = frames[frame_idx]

        xs_ = [p[2] for p in particles]
        ys_ = [p[3] for p in particles]
        colors = [POLYMER_COLORS[p[1] % 4][p[4] % 4] for p in particles]

        scat.set_offsets(np.column_stack([xs_, ys_]))
        scat.set_color(colors)

        # Remove old bond lines
        for ln in bond_lines:
            ln.remove()
        bond_lines = []

        # Draw backbone bonds: consecutive ids within same copy, dist ≤ √2
        by_id = {p[0]: p for p in particles}
        for p in particles:
            pid, poly, px, py, copy_ = p
            nid = pid + 1
            if nid in by_id:
                np2 = by_id[nid]
                if np2[4] == copy_:
                    dist2 = (np2[2] - px)**2 + (np2[3] - py)**2
                    if dist2 < 2.5:
                        ln, = ax.plot([px, np2[2]], [py, np2[3]],
                                      color=BACKBONE_COLOR, lw=BOND_LW, zorder=2)
                        bond_lines.append(ln)

        frame_text.set_text(
            f"frame {frame_idx * args.skip}  |  N={len(particles)}"
        )
        return [scat, frame_text] + bond_lines

    ani = animation.FuncAnimation(
        fig, update,
        frames=len(frames),
        interval=1000 // args.fps,
        blit=False,
    )

    if args.output:
        print(f"Saving animation to {args.output} ...", flush=True)
        writer = (
            animation.FFMpegWriter(fps=args.fps)
            if args.output.endswith(".mp4")
            else animation.PillowWriter(fps=args.fps)
        )
        ani.save(args.output, writer=writer)
        print("Done.")
    else:
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
