#!/usr/bin/env python3
"""
visualize_nucleolus.py

Reads a nucleolus trajectory file produced by run_nucleolus and generates a
3-panel animated figure:
  Left:   particle positions coloured by polymer type with gradient overlay
  Top-right:  energy vs simulation step
  Bottom-right: cumulative exited complexes vs simulation step

Usage
-----
    python3 visualize_nucleolus.py <traj_file> [OPTIONS]

Options
    --output FILE            Save animation to FILE (.mp4 or .gif). Default: display.
    --fps N                  Frames per second (default 5).
    --gradient-length L      Draw gradient overlay from x=0 to x=L (default: auto).
    --width W                Column width in y (default: auto).
    --skip N                 Render only every N-th frame (default 1).
    --title TEXT             Figure title.

Trajectory format (written by run_nucleolus)
--------------------------------------------
    <N_particles>
    step=S energy=E exited=X L=LL W=WW nCopies=C
    <id> <poly_type> <x> <y> <copy>
    ...   (repeated N_particles times per frame)
"""

import argparse
import re
import sys

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import numpy as np


# --------------------------------------------------------------------------- #
# Colour palette: 4 polymer types, 4 copies each
# --------------------------------------------------------------------------- #
_POLY_COLORS = [
    ["#e6194b", "#f58231", "#ffe119", "#bfef45"],   # type 0
    ["#3cb44b", "#42d4f4", "#4363d8", "#911eb4"],   # type 1
    ["#f032e6", "#a9a9a9", "#fabed4", "#aaffc3"],   # type 2
    ["#ffd8b1", "#dcbeff", "#469990", "#000075"],   # type 3
]
_BACKBONE_COLOR = "#333333"
_BOND_LW = 1.5


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

def _kv(header, key, default=None):
    """Extract key=value from a header string."""
    m = re.search(rf'{key}=([^\s]+)', header)
    if m:
        return m.group(1)
    return default


def parse_traj(path):
    """
    Return a list of frame dicts, each with:
        step, energy, exited, L, W, nCopies,
        particles: list of (id, poly_type, x, y, copy)
    """
    frames = []
    with open(path) as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            n = int(line)
        except ValueError:
            i += 1
            continue

        # Header line
        i += 1
        if i >= len(lines):
            break
        hdr = lines[i].strip()
        i += 1

        step    = int(_kv(hdr, 'step', 0))
        energy  = float(_kv(hdr, 'energy', 0.0))
        exited  = int(_kv(hdr, 'exited', 0))
        L       = float(_kv(hdr, 'L', 60))
        W       = float(_kv(hdr, 'W', 10))
        nCopies = int(_kv(hdr, 'nCopies', 4))

        particles = []
        for _ in range(n):
            if i >= len(lines):
                break
            parts = lines[i].strip().split()
            i += 1
            if len(parts) >= 5:
                particles.append((
                    int(parts[0]),   # id
                    int(parts[1]),   # poly_type
                    float(parts[2]), # x
                    float(parts[3]), # y
                    int(parts[4]),   # copy
                ))

        frames.append(dict(
            step=step, energy=energy, exited=exited,
            L=L, W=W, nCopies=nCopies, particles=particles,
        ))

    return frames


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description="Animate nucleolus trajectory.")
    ap.add_argument("traj", help="Trajectory file (output of run_nucleolus)")
    ap.add_argument("--output", default=None,
                    help="Save animation to file (.mp4 or .gif)")
    ap.add_argument("--fps", type=int, default=5)
    ap.add_argument("--gradient-length", type=float, default=None,
                    help="Column length L for gradient overlay")
    ap.add_argument("--width", type=float, default=None,
                    help="Column width W (y-axis range)")
    ap.add_argument("--skip", type=int, default=1,
                    help="Render every N-th frame")
    ap.add_argument("--title", default="Nucleolus assembly simulation")
    args = ap.parse_args()

    print(f"Reading {args.traj} ...", flush=True)
    all_frames = parse_traj(args.traj)
    if not all_frames:
        print("No frames found in trajectory.", file=sys.stderr)
        sys.exit(1)
    frames = all_frames[::args.skip]
    print(f"  {len(all_frames)} frames total, rendering {len(frames)}.", flush=True)

    # Pull time-series for background plots (all frames, not just rendered ones)
    ts_steps   = [f['step']   for f in all_frames]
    ts_energy  = [f['energy'] for f in all_frames]
    ts_exited  = [f['exited'] for f in all_frames]

    # Axis limits from data
    all_x = [p[2] for fr in frames for p in fr['particles']]
    all_y = [p[3] for fr in frames for p in fr['particles']]
    x_min = -1.0
    x_max = max(all_x) + 2.0

    grad_L = args.gradient_length if args.gradient_length else frames[0]['L']
    col_W  = args.width if args.width else frames[0]['W']
    y_min  = -0.5
    y_max  = col_W + 0.5

    # ---------------------------------------------------------------------- #
    # Figure: left = animation, right = stacked plots
    # ---------------------------------------------------------------------- #
    fig = plt.figure(figsize=(14, 6))
    fig.suptitle(args.title, fontsize=11)

    gs = fig.add_gridspec(2, 2, width_ratios=[2.2, 1], hspace=0.45, wspace=0.35)
    ax_anim  = fig.add_subplot(gs[:, 0])   # full left column
    ax_ener  = fig.add_subplot(gs[0, 1])   # top-right
    ax_exits = fig.add_subplot(gs[1, 1])   # bottom-right

    # --- Animation panel ---
    ax_anim.set_aspect("equal")
    ax_anim.set_xlim(x_min, x_max)
    ax_anim.set_ylim(y_min, y_max)
    ax_anim.set_xlabel("x  (column axis)", fontsize=9)
    ax_anim.set_ylabel("y", fontsize=9)

    # Gradient background
    grad_img = np.tile(
        np.clip(np.linspace(0, 1, 256), 0, 1)[np.newaxis, :],
        (10, 1)
    )
    ax_anim.imshow(grad_img,
                   extent=[0, grad_L, y_min, y_max],
                   origin="lower", aspect="auto",
                   cmap="Blues", alpha=0.20, zorder=0)
    ax_anim.axvline(0,      color="black",     lw=1.5, zorder=1)
    ax_anim.axvline(grad_L, color="steelblue", lw=1.5, linestyle="--", zorder=1)
    ax_anim.text(grad_L + 0.3, y_min + 0.3, f"L={grad_L:.0f}",
                 color="steelblue", fontsize=7)

    # Polymer-type legend
    legend_patches = [
        mpatches.Patch(color=_POLY_COLORS[t][0], label=f"Polymer type {t}")
        for t in range(4)
    ]
    ax_anim.legend(handles=legend_patches, loc="upper right", fontsize=7,
                   framealpha=0.7)

    scat = ax_anim.scatter([], [], s=180, zorder=3,
                            edgecolors="k", linewidths=0.4)
    bond_lines = []
    frame_text = ax_anim.text(0.02, 0.97, "", transform=ax_anim.transAxes,
                               fontsize=8, va="top", family="monospace")

    # --- Energy panel ---
    ax_ener.plot(ts_steps, ts_energy, color="#555555", lw=0.8, alpha=0.4,
                 label="all steps")
    ener_marker, = ax_ener.plot([], [], 'o', color="#e6194b", ms=5, zorder=5)
    ener_vline   = ax_ener.axvline(0, color="#e6194b", lw=0.8, ls="--", alpha=0.7)
    ax_ener.set_xlabel("Step", fontsize=8)
    ax_ener.set_ylabel("Energy", fontsize=8)
    ax_ener.set_title("System energy", fontsize=8)
    ax_ener.tick_params(labelsize=7)
    ax_ener.set_xlim(min(ts_steps), max(ts_steps))

    # --- Exits panel ---
    ax_exits.plot(ts_steps, ts_exited, color="#555555", lw=0.8, alpha=0.4)
    exits_marker, = ax_exits.plot([], [], 'o', color="#3cb44b", ms=5, zorder=5)
    exits_vline   = ax_exits.axvline(0, color="#3cb44b", lw=0.8, ls="--", alpha=0.7)
    ax_exits.set_xlabel("Step", fontsize=8)
    ax_exits.set_ylabel("Count", fontsize=8)
    ax_exits.set_title("Cumulative exited complexes", fontsize=8)
    ax_exits.tick_params(labelsize=7)
    ax_exits.set_xlim(min(ts_steps), max(ts_steps))

    # ---------------------------------------------------------------------- #
    # Update function
    # ---------------------------------------------------------------------- #
    def update(frame_idx):
        nonlocal bond_lines
        fr = frames[frame_idx]
        pts = fr['particles']

        xs     = [p[2] for p in pts]
        ys     = [p[3] for p in pts]
        colors = [_POLY_COLORS[p[1] % 4][p[4] % 4] for p in pts]

        scat.set_offsets(np.column_stack([xs, ys]))
        scat.set_color(colors)

        # Backbone bonds: consecutive ids within same copy, dist² ≤ 2.5
        for ln in bond_lines:
            ln.remove()
        bond_lines = []
        by_id = {p[0]: p for p in pts}
        for p in pts:
            pid, _, px, py, copy_ = p
            nxt = by_id.get(pid + 1)
            if nxt and nxt[4] == copy_:
                dist2 = (nxt[2] - px)**2 + (nxt[3] - py)**2
                if dist2 < 2.5:
                    ln, = ax_anim.plot([px, nxt[2]], [py, nxt[3]],
                                       color=_BACKBONE_COLOR, lw=_BOND_LW, zorder=2)
                    bond_lines.append(ln)

        frame_text.set_text(
            f"step {fr['step']}  E={fr['energy']:.1f}  exited={fr['exited']}"
        )

        # Update scalar plot markers
        s = fr['step']
        ener_marker.set_data([s], [fr['energy']])
        ener_vline.set_xdata([s, s])
        exits_marker.set_data([s], [fr['exited']])
        exits_vline.set_xdata([s, s])

        return [scat, frame_text, ener_marker, ener_vline,
                exits_marker, exits_vline] + bond_lines

    ani = animation.FuncAnimation(
        fig, update,
        frames=len(frames),
        interval=1000 // args.fps,
        blit=False,
    )

    if args.output:
        print(f"Saving to {args.output} ...", flush=True)
        if args.output.endswith(".mp4"):
            writer = animation.FFMpegWriter(fps=args.fps)
        else:
            writer = animation.PillowWriter(fps=args.fps)
        ani.save(args.output, writer=writer)
        print("Saved.")
    else:
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
