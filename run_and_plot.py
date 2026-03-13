#!/usr/bin/env python3
"""
Run and plot a HierarchicalAssembly simulation.

Usage:
    python run_and_plot.py                          # run with defaults
    python run_and_plot.py --e1 6.0 --n0 64        # custom parameters
    python run_and_plot.py --p 128 --L 50           # specify particles & box size directly
    python run_and_plot.py --no-run                 # just plot existing output

Particle / box-size shortcuts:
    --p   total number of particles (must be a multiple of --n0); sets ncopies = p / n0
    --L   box side length (periodic boundary conditions); sets dens = p / L²
    If --p or --L are omitted the legacy --ncopies / --dens flags are used instead.

Dependencies: numpy, matplotlib
    pip install numpy matplotlib
"""

import argparse
import subprocess
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.patches as patches
import matplotlib.animation as animation
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Run and plot HierarchicalAssembly")

    # Command-line argument passed directly to run_hier
    p.add_argument("--e1", type=float, default=8.0,
                   help="Strongest binding energy (default: 8.0)")

    # Input-file parameters
    p.add_argument("--n0", type=int, default=16, choices=[4, 16, 64, 256, 1024],
                   help="Particles in base cell / assembly target size (default: 16)")
    p.add_argument("--ncopies", type=int, default=6,
                   help="Number of copies of each particle type (default: 6)")
    p.add_argument("--nsteps", type=int, default=100,
                   help="Simulation steps; output is written each step (default: 100)")
    p.add_argument("--nsweep", type=int, default=400,
                   help="Monte Carlo sweeps per step (default: 400)")
    p.add_argument("--dens", type=float, default=0.05,
                   help="Volume fraction / density (default: 0.05)")
    p.add_argument("--filehead", type=str, default="hier",
                   help="Prefix for output files (default: hier)")

    # Convenience shortcuts: specify particles + box size directly
    p.add_argument("--p", type=int, default=None,
                   help="Total particles in simulation (must be a multiple of --n0); "
                        "overrides --ncopies")
    p.add_argument("--L", type=float, default=None,
                   help="Box side length (periodic boundary conditions); overrides --dens")

    # Custom bond matrix
    p.add_argument("--bond-file", type=str, default=None,
                   help="Path to a custom bond file; uses run_custom instead of run_hier")
    p.add_argument("--gen-bonds", action="store_true",
                   help="Generate a random bond file (normal distribution) and use run_custom")
    p.add_argument("--bond-seed", type=int, default=None,
                   help="Random seed for --gen-bonds (default: random)")
    p.add_argument("--bond-std", type=float, default=0.3,
                   help="Std dev of bond strengths as a fraction of e1 for --gen-bonds (default: 0.3)")

    # Script behaviour
    p.add_argument("--no-run", action="store_true",
                   help="Skip the simulation and just plot existing output files")
    p.add_argument("--save-only", action="store_true",
                   help="Save plot to PNG without opening an interactive window")

    return p.parse_args()


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def write_input_file(path, filehead, n0, ncopies, nsteps, nsweep, dens):
    with open(path, "w") as f:
        f.write(f"{filehead}      # filehead\n")
        f.write(f"{n0}        # n\n")
        f.write(f"{ncopies}         # number of copies\n")
        f.write(f"{nsteps}        # number of steps\n")
        f.write(f"{nsweep}       # number of sweeps per step\n")
        f.write(f"{dens}      # density\n")


def run_simulation(exe, input_file, e1):
    cmd = [exe, input_file, str(e1)]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f"Simulation failed (return code {result.returncode})")
    print("Simulation complete.\n")


def generate_bond_file(n0, e1, path, std_frac=0.3, seed=None):
    """Write a bond file with normally distributed random bond strengths.

    Each adjacent pair in the n0-particle target grid gets an independent
    sample from N(e1, e1*std_frac), clipped to (0.01, inf).
    Returns a dict {(i,j): energy} mirroring what compute_bond_table would return.
    """
    rng = np.random.default_rng(seed)
    l0 = round(math.sqrt(n0))
    bonds = {}

    def gidx(col, row): return l0*row + col

    for col in range(l0 - 1):        # east bonds
        for row in range(l0):
            val = max(0.01, rng.normal(e1, e1 * std_frac))
            p1, p2 = gidx(col, row), gidx(col+1, row)
            bonds[(p1, p2)] = val
            bonds[(p2, p1)] = val

    for col in range(l0):             # north bonds
        for row in range(l0 - 1):
            val = max(0.01, rng.normal(e1, e1 * std_frac))
            p1, p2 = gidx(col, row), gidx(col, row+1)
            bonds[(p1, p2)] = val
            bonds[(p2, p1)] = val

    with open(path, "w") as f:
        f.write(f"# Custom bond file  n0={n0}  e1={e1}  std_frac={std_frac}"
                f"  seed={seed}\n")
        f.write("# particle_i  particle_j  energy\n")
        written = set()
        for (p1, p2), val in bonds.items():
            if (p2, p1) not in written:
                f.write(f"{p1} {p2} {val:.6f}\n")
                written.add((p1, p2))

    print(f"Bond file written to {path}  ({len(written)} bonds)")
    return bonds


def load_bond_file(path, n0):
    """Read a bond file and return {(i,j): energy} (symmetric)."""
    l0 = round(math.sqrt(n0))
    bonds = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            pi, pj, val = int(parts[0]), int(parts[1]), float(parts[2])
            bonds[(pi, pj)] = val
            bonds[(pj, pi)] = val
    return bonds


def parse_stats(statsfile):
    """Return (steps, energy, fragment_hist) as numpy arrays.

    Stats file format (after 1-line header):
        step  total_energy  hist[1]  hist[2]  ...  hist[n0]
    """
    rows = []
    with open(statsfile) as f:
        f.readline()  # skip header
        for line in f:
            line = line.strip()
            if line:
                rows.append([float(v) for v in line.split()])
    data = np.array(rows)
    return data[:, 0], data[:, 1], data[:, 2:]


def parse_traj(trajfile):
    """Return (n_particles, box_length, n0, frames).

    Trajectory file format:
        <description line>
        n_particles  box_length  n0          <- first frame header (3 values)
                                             <- blank line
        0  x  y  0.0000                      <- n_particles lines
        ...
        n_particles                          <- subsequent frame headers (1 value)
                                             <- blank line
        0  x  y  0.0000
        ...
    """
    frames = []
    with open(trajfile) as f:
        f.readline()  # description

        # First frame header has 3 values: n, L, n0
        first_header = f.readline().split()
        n_particles = int(first_header[0])
        box_length = float(first_header[1])
        n0 = int(first_header[2])
        f.readline()  # blank line

        frame = _read_frame(f, n_particles)
        if frame is not None:
            frames.append(frame)

        while True:
            header = f.readline()
            if not header:
                break
            header = header.strip()
            if not header:
                header = f.readline().strip()
            if not header:
                break
            f.readline()  # blank line
            frame = _read_frame(f, n_particles)
            if frame is None:
                break
            frames.append(frame)

    return n_particles, box_length, n0, frames


def _read_frame(f, n_particles):
    """Read n_particles lines; return (n_particles, 2) array of [x, y] or None."""
    coords = []
    for _ in range(n_particles):
        line = f.readline()
        if not line:
            return None
        parts = line.split()
        if len(parts) < 3:
            return None
        coords.append([float(parts[1]), float(parts[2])])
    return np.array(coords)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

# Approximate MATLAB cmocean 'thermal' palette (particle identity colours)
_THERMAL_COLOURS = [
    "#04082e", "#2a1654", "#5e1269", "#8f0f6e",
    "#bb2564", "#dc4a50", "#f17733", "#f5a535",
    "#f7d03c", "#fcf4ae",
]
THERMAL = LinearSegmentedColormap.from_list("thermal", _THERMAL_COLOURS, N=256)

def _slot_colours(n0):
    return [THERMAL(i / max(n0 - 1, 1)) for i in range(n0)]


def compute_bond_table(n0, e1, custom_bonds=None):
    """Return {(id1, id2): energy} for every native bond in the target structure.
    If custom_bonds dict is provided it is used directly instead of the
    hierarchical formula (mirrors run_hier.cpp)."""
    if custom_bonds is not None:
        return custom_bonds
    l0 = round(math.sqrt(n0))
    e2, e3, e4, e5 = e1/2, e1/4, e1/8, e1/16
    bonds = {}

    def gidx(i, j): return l0*j + i

    for i in range(l0 - 1):          # East bonds
        for j in range(l0):
            if   (i+1)%16 == 0 and n0 >= 1024: val = e5
            elif (i+1)%8  == 0 and n0 >= 256:  val = e4
            elif (i+1)%4  == 0 and n0 >= 64:   val = e3
            elif (i+1)%2  == 0 and n0 >= 16:   val = e2
            else:                               val = e1
            p1, p2 = gidx(i, j), gidx(i+1, j)
            bonds[(p1, p2)] = val
            bonds[(p2, p1)] = val

    for i in range(l0):               # North bonds
        for j in range(l0 - 1):
            if   (j+1)%16 == 0 and n0 >= 1024: val = e5
            elif (j+1)%8  == 0 and n0 >= 256:  val = e4
            elif (j+1)%4  == 0 and n0 >= 64:   val = e3
            elif (j+1)%2  == 0 and n0 >= 16:   val = e2
            else:                               val = e1
            p1, p2 = gidx(i, j), gidx(i, j+1)
            bonds[(p1, p2)] = val
            bonds[(p2, p1)] = val

    return bonds


def _bond_segments_pbc(cx1, cy1, cx2, cy2, L):
    """Return 1 or 2 line segments for a bond, respecting periodic boundaries."""
    dx = cx2 - cx1
    dy = cy2 - cy1
    if dx >  L/2: dx -= L
    elif dx < -L/2: dx += L
    if dy >  L/2: dy -= L
    elif dy < -L/2: dy += L
    wraps = abs((cx2 - cx1) - dx) > 0.1 or abs((cy2 - cy1) - dy) > 0.1
    if not wraps:
        return [[(cx1, cy1), (cx2, cy2)]]
    # Two half-segments exiting toward opposite walls (axes will clip them)
    return [
        [(cx1, cy1), (cx1 + dx, cy1 + dy)],
        [(cx2, cy2), (cx2 - dx, cy2 - dy)],
    ]


def make_plots(steps, energy, fragment_hist, n_particles, box_length, n0,
               frames, filehead, e1, custom_bonds=None):

    bond_table = compute_bond_table(n0, e1, custom_bonds)

    # Coupling matrix: entry [i,j] = native bond energy between identities i and j
    coupling_matrix = np.zeros((n0, n0))
    for (i, j), val in bond_table.items():
        coupling_matrix[i, j] = val

    colours = _slot_colours(n0)
    n_frames = len(frames)
    L = float(box_length)

    # Shared colormap for bond lines and coupling matrix
    bond_vals = sorted(set(bond_table.values()))
    bond_cmap = plt.cm.plasma
    bond_norm = mcolors.Normalize(vmin=bond_vals[0] if bond_vals else 0, vmax=e1)

    # ---- Layout: lattice (left) | energy (top-right) | matrix (bottom-right) --
    fig = plt.figure(figsize=(15, 7))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1], hspace=0.45, wspace=0.35)
    ax_lat = fig.add_subplot(gs[:, 0])
    ax_en  = fig.add_subplot(gs[0, 1])
    ax_cm  = fig.add_subplot(gs[1, 1])

    fig.suptitle(
        f"HierarchicalAssembly  |  n0={n0}  nParticles={n_particles}"
        f"  e1={e1}  dens={n_particles / L**2:.3f}",
        fontsize=11,
    )

    # ---- Lattice panel ----------------------------------------------------
    ax_lat.set_facecolor("#111111")
    ax_lat.set_xlim(0, L)
    ax_lat.set_ylim(0, L)
    ax_lat.set_aspect("equal")
    ax_lat.set_xlabel("x")
    ax_lat.set_ylabel("y")
    lat_title = ax_lat.set_title("Step 0", fontsize=10)

    # Bond lines drawn below circles
    bond_lc = LineCollection([], linewidths=2.5, alpha=0.85, zorder=1)
    ax_lat.add_collection(bond_lc)

    # Circles for particles
    circles = []
    for i, (x, y) in enumerate(frames[0]):
        circ = patches.Circle(
            (x + 0.5, y + 0.5), radius=0.35,
            facecolor=colours[i % n0], edgecolor="white", linewidth=0.5, zorder=2,
        )
        ax_lat.add_patch(circ)
        circles.append(circ)

    # ---- Energy panel -----------------------------------------------------
    ax_en.plot(steps, energy, lw=1, color="steelblue", alpha=0.7, label="Energy")
    running_avg = np.cumsum(energy) / (np.arange(len(energy)) + 1)
    ax_en.plot(steps, running_avg, lw=1.5, color="tomato", label="Running avg")
    ax_en.set_xlabel("Step")
    ax_en.set_ylabel("Total energy")
    ax_en.set_title("Energy vs time")
    ax_en.legend(fontsize=8)
    ax_en.grid(True, alpha=0.3)
    vline = ax_en.axvline(steps[0], color="gold", lw=1.0, linestyle="--", alpha=0.9)

    # ---- Coupling matrix panel --------------------------------------------
    masked_cm = np.ma.masked_where(coupling_matrix == 0, coupling_matrix)
    ax_cm.set_facecolor("#222222")
    cm_img = ax_cm.imshow(
        masked_cm, cmap=bond_cmap, norm=bond_norm,
        origin="upper", aspect="auto", interpolation="nearest",
    )
    fig.colorbar(cm_img, ax=ax_cm, label="Bond energy")
    ax_cm.set_title("Bond strength matrix", fontsize=9)
    ax_cm.set_xlabel("Particle identity")
    ax_cm.set_ylabel("Particle identity")
    tick_step = max(1, n0 // 8)
    ticks = list(range(0, n0, tick_step))
    ax_cm.set_xticks(ticks)
    ax_cm.set_yticks(ticks)

    plt.tight_layout()

    # ---- Animation --------------------------------------------------------
    MAX_FRAMES = 200
    indices = (np.linspace(0, n_frames - 1, MAX_FRAMES, dtype=int)
               if n_frames > MAX_FRAMES else np.arange(n_frames))
    iL = int(round(L))

    def update(k):
        frame_idx = indices[k]
        coords = frames[frame_idx]

        # Move circles
        for i, circ in enumerate(circles):
            circ.set_center((coords[i, 0] + 0.5, coords[i, 1] + 0.5))

        # Rebuild bond segments for this frame
        pos_map = {(int(coords[i, 0]), int(coords[i, 1])): i for i in range(n_particles)}
        segments, seg_colors = [], []
        for i in range(n_particles):
            xi, yi = int(coords[i, 0]), int(coords[i, 1])
            id1 = i % n0
            # Cardinal + diagonal offsets; each covers unique pairs (no double-counting)
            for dxi, dyi in ((1, 0), (0, 1), (1, 1), (1, -1)):
                j = pos_map.get(((xi + dxi) % iL, (yi + dyi) % iL))
                if j is None:
                    continue
                val = bond_table.get((id1, j % n0), 0)
                if val <= 0:
                    continue
                for seg in _bond_segments_pbc(xi+0.5, yi+0.5,
                                              coords[j,0]+0.5, coords[j,1]+0.5, L):
                    segments.append(seg)
                    seg_colors.append(bond_cmap(bond_norm(val)))

        bond_lc.set_segments(segments)
        bond_lc.set_colors(seg_colors)
        lat_title.set_text(f"Step {int(steps[frame_idx])}")
        vline.set_xdata([steps[frame_idx], steps[frame_idx]])
        return circles + [bond_lc, lat_title, vline]

    anim = animation.FuncAnimation(
        fig, update, frames=len(indices), interval=50, blit=False,
    )

    # Spacebar toggles pause/play
    is_paused = [False]
    def on_key(event):
        if event.key == ' ':
            is_paused[0] = not is_paused[0]
            if is_paused[0]:
                anim.pause()
            else:
                anim.resume()
    fig.canvas.mpl_connect('key_press_event', on_key)

    plt.show()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    exe = os.path.join(script_dir, "run_hier")
    input_file = os.path.join(script_dir, f"input_{args.filehead}.txt")
    statsfile = os.path.join(script_dir, f"{args.filehead}_stats.txt")
    trajfile = os.path.join(script_dir, f"{args.filehead}_traj.txt")

    # Resolve ncopies and dens from --p / --L if supplied
    ncopies = args.ncopies
    dens = args.dens

    if args.p is not None:
        if args.p % args.n0 != 0:
            sys.exit(f"Error: --p ({args.p}) must be a multiple of --n0 ({args.n0})")
        ncopies = args.p // args.n0

    if args.L is not None:
        n_particles = args.n0 * ncopies
        dens = n_particles / args.L ** 2
        print(f"Box: L={args.L}, nParticles={n_particles}, dens={dens:.6f}")

    # Resolve bond file / custom mode
    custom_bonds = None
    bond_file = args.bond_file

    if args.gen_bonds:
        bond_file = os.path.join(script_dir, f"{args.filehead}_bonds.txt")
        custom_bonds = generate_bond_file(
            args.n0, args.e1, bond_file,
            std_frac=args.bond_std, seed=args.bond_seed,
        )
    elif bond_file is not None:
        custom_bonds = load_bond_file(bond_file, args.n0)

    use_custom = bond_file is not None
    exe_to_use = os.path.join(script_dir, "run_custom") if use_custom else exe

    if not args.no_run:
        write_input_file(input_file, args.filehead, args.n0, ncopies,
                         args.nsteps, args.nsweep, dens)
        if use_custom:
            run_simulation(exe_to_use, input_file, bond_file)
        else:
            run_simulation(exe, input_file, args.e1)

    print("Parsing output files...")
    steps, energy, fragment_hist = parse_stats(statsfile)
    n_particles, box_length, n0, frames = parse_traj(trajfile)
    print(f"  {len(frames)} frames  |  {n_particles} particles  |  box={box_length}  |  n0={n0}")

    make_plots(steps, energy, fragment_hist, n_particles, box_length, n0,
               frames, args.filehead, args.e1, custom_bonds=custom_bonds)


if __name__ == "__main__":
    main()
