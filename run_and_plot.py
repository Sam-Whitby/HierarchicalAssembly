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
    p.add_argument("--polymer", type=int, default=None,
                   help="Polymer size: auto-generate backbone bonds and run run_polymer. "
                        "n0 is set to this value, one copy.")
    p.add_argument("--weak-e", type=float, default=1.0,
                   help="Mean weak coupling energy magnitude (default: 1.0)")
    p.add_argument("--weak-std", type=float, default=0.5,
                   help="Std dev of weak couplings (default: 0.5). Can be > mean (allows negative = repulsive).")
    p.add_argument("--weak-seed", type=int, default=None,
                   help="Random seed for weak coupling generation")
    p.add_argument("--sim-seed", type=int, default=None,
                   help="Random seed for the C++ simulation RNG")
    p.add_argument("--hilbert-A", type=float, default=None,
                   help="Coupling strength A for Hilbert-curve cardinal-adjacent pairs "
                        "(d=1 in Hilbert grid). Activates Hilbert-local coupling mode.")
    p.add_argument("--hilbert-B", type=float, default=None,
                   help="Coupling strength B for Hilbert-curve diagonal-adjacent pairs "
                        "(d=√2 in Hilbert grid). Used with --hilbert-A.")
    p.add_argument("--hier-red", type=float, default=None,
                   help="Coupling strength for RED bonds (finest level: pairs within same "
                        "group of 4 consecutive Hilbert-chain positions). Activates "
                        "hierarchical Hilbert coupling mode.")
    p.add_argument("--hier-green", type=float, default=None,
                   help="Coupling strength for GREEN bonds (middle level: same top-level "
                        "quarter of chain but different group of 4).")
    p.add_argument("--hier-blue", type=float, default=None,
                   help="Coupling strength for BLUE bonds (coarsest level: different "
                        "top-level quarters of the chain).")
    p.add_argument("--boltzmann", action="store_true",
                   help="Enable Boltzmann validation and canonical-state diagram (polymer mode only)")
    p.add_argument("--denature", type=int, default=None,
                   help="Number of equilibration steps with all weak bonds zeroed (backbone "
                        "confinement only). Runs before the main simulation; the final "
                        "configuration is used as the starting point for the main run.")

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


def snake_path(n0):
    """Return the snake-path bond list for a polymer of n0 particles in a sqrt(n0)×sqrt(n0) grid.
    Returns list of (i, j) pairs that are adjacent in the snake walk."""
    l0 = round(math.sqrt(n0))
    assert l0 * l0 == n0, f"n0={n0} must be a perfect square for grid polymer"
    # Build snake path: row by row, alternating direction
    path = []
    for row in range(l0):
        cols = range(l0) if row % 2 == 0 else range(l0-1, -1, -1)
        for col in cols:
            path.append(l0 * row + col)
    return [(path[k], path[k+1]) for k in range(len(path)-1)]


def generate_polymer_bondfile(n0, e_backbone, e_weak, std_weak, path, seed=None):
    """Generate an extended bond file using confinement-based backbone representation.

    Backbone bonds are NOT encoded as directional attractions.  Instead, for each
    backbone-bonded pair (i,j) a large repulsive entry is added to wD2 and wDsq5:
        coupling[i,j] -= e_backbone   (coupling << 0  →  physical energy >> 0  →  repulsive)
    Combined with the hard-sphere exclusion (d=0 → INF) this confines every bonded
    pair to d ∈ {1, √2} without any extra C++ machinery.

    A wD0 matrix is returned for display purposes only (d=0 is already INF in C++).

    Returns dict: {'D0': mat, 'D1': mat, 'Dsq2': mat, 'D2': mat, 'Dsq5': mat}
    All matrices are symmetric n0×n0 numpy arrays (physical energy = -coupling).
    """
    rng = np.random.default_rng(seed)

    backbone_bonds = snake_path(n0)

    def sym_matrix(n, mean, std):
        mat = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                val = rng.normal(mean, std)
                mat[i, j] = val
                mat[j, i] = val
        return mat

    wD1   = sym_matrix(n0, e_weak, std_weak)
    wDsq2 = sym_matrix(n0, e_weak, std_weak)
    wD2   = sym_matrix(n0, e_weak, std_weak)
    wDsq5 = sym_matrix(n0, e_weak, std_weak)
    wD0   = np.zeros((n0, n0))  # display-only; C++ returns INF at d=0 regardless

    # Bake confinement into wD2 and wDsq5 for every backbone-bonded pair
    for (pi, pj) in backbone_bonds:
        wD2[pi, pj]   -= e_backbone;  wD2[pj, pi]   -= e_backbone
        wDsq5[pi, pj] -= e_backbone;  wDsq5[pj, pi] -= e_backbone
        wD0[pi, pj]   -= e_backbone;  wD0[pj, pi]   -= e_backbone  # visual only

    def write_matrix(f, mat, tag):
        f.write(f"{tag}\n")
        for i in range(n0):
            for j in range(i, n0):
                f.write(f"{i} {j} {mat[i,j]:.6f}\n")
        f.write(f"{tag}_END\n\n")

    with open(path, 'w') as f:
        f.write(f"# Extended polymer bond file  n0={n0}  e_backbone={e_backbone}"
                f"  e_weak={e_weak}  std_weak={std_weak}  seed={seed}\n")
        f.write("# Backbone confinement: repulsive wD2/wDsq5 entries for bonded pairs\n")
        f.write("# (no directional attraction; hard sphere + repulsion confines chain)\n\n")

        f.write("BACKBONE\n")
        f.write("BACKBONE_END\n\n")  # empty — confinement handled via wD2/wDsq5

        write_matrix(f, wD1,   "WEAK_D1")
        write_matrix(f, wDsq2, "WEAK_DSQRT2")
        write_matrix(f, wD2,   "WEAK_D2")
        write_matrix(f, wDsq5, "WEAK_DSQRT5")

    print(f"Extended bond file written to {path}  (confinement via wD2/wDsq5)")
    return {'D0': wD0, 'D1': wD1, 'Dsq2': wDsq2, 'D2': wD2, 'Dsq5': wDsq5}


def generate_hilbert_bondfile(n0, e_backbone, A, B, path):
    """Generate a bond file with Hilbert-curve-local couplings only.

    Non-zero weak coupling is assigned only between particles that are spatially
    adjacent in the Hilbert curve grid used for initialisation:
      - Hilbert cardinal neighbours (grid distance = 1):   wD1[i,j]   = A
      - Hilbert diagonal neighbours (grid distance = √2):  wDsq2[i,j] = B

    Backbone-bonded pairs (always at Hilbert distance 1, chain-consecutive) are
    excluded from A/B — they are already confined by repulsion in wD2/wDsq5.
    All other matrix entries are zero.

    Requires l0 = sqrt(n0) to be a power of 2 (Hilbert curve is defined for those).

    Returns dict: {'D0': mat, 'D1': mat, 'Dsq2': mat, 'D2': mat, 'Dsq5': mat}
    """
    l0 = round(math.sqrt(n0))
    if not _is_power_of_2(l0):
        raise ValueError(f"Hilbert-local mode requires l0=sqrt(n0) to be a power of 2; "
                         f"got n0={n0}, l0={l0}")

    backbone_bonds = snake_path(n0)
    backbone_set = set()
    for (pi, pj) in backbone_bonds:
        backbone_set.add((pi, pj))
        backbone_set.add((pj, pi))

    # Map particle_id -> (hx, hy) in the Hilbert grid
    chain_order = _chain_order(n0)
    sfc_pos = _space_filling_positions(l0)  # Hilbert positions in chain traversal order
    hpos = {}
    for k, pid in enumerate(chain_order):
        hpos[pid] = sfc_pos[k]

    wD1   = np.zeros((n0, n0))
    wDsq2 = np.zeros((n0, n0))
    wD2   = np.zeros((n0, n0))
    wDsq5 = np.zeros((n0, n0))
    wD0   = np.zeros((n0, n0))  # display-only

    # Assign A/B to Hilbert-adjacent non-backbone pairs
    n_d1, n_dsq2 = 0, 0
    for i in range(n0):
        for j in range(i + 1, n0):
            if (i, j) in backbone_set:
                continue
            dx = hpos[j][0] - hpos[i][0]
            dy = hpos[j][1] - hpos[i][1]
            d2 = dx * dx + dy * dy
            if d2 == 1:
                wD1[i, j] = wD1[j, i] = A;  n_d1 += 1
            elif d2 == 2:
                wDsq2[i, j] = wDsq2[j, i] = B;  n_dsq2 += 1

    # Backbone confinement: repulsion at d=2 and d=√5
    for (pi, pj) in backbone_bonds:
        wD2[pi, pj]   -= e_backbone;  wD2[pj, pi]   -= e_backbone
        wDsq5[pi, pj] -= e_backbone;  wDsq5[pj, pi] -= e_backbone
        wD0[pi, pj]   -= e_backbone;  wD0[pj, pi]   -= e_backbone  # visual only

    def write_matrix(f, mat, tag):
        f.write(f"{tag}\n")
        for i in range(n0):
            for j in range(i, n0):
                f.write(f"{i} {j} {mat[i,j]:.6f}\n")
        f.write(f"{tag}_END\n\n")

    with open(path, 'w') as f:
        f.write(f"# Hilbert-local bond file  n0={n0}  e_backbone={e_backbone}"
                f"  A={A}  B={B}\n")
        f.write(f"# Non-zero pairs: {n_d1} cardinal (A), {n_dsq2} diagonal (B)\n\n")
        f.write("BACKBONE\n")
        f.write("BACKBONE_END\n\n")
        write_matrix(f, wD1,   "WEAK_D1")
        write_matrix(f, wDsq2, "WEAK_DSQRT2")
        write_matrix(f, wD2,   "WEAK_D2")
        write_matrix(f, wDsq5, "WEAK_DSQRT5")

    print(f"Hilbert-local bond file written to {path}  "
          f"({n_d1} cardinal pairs with A={A}, {n_dsq2} diagonal pairs with B={B})")
    return {'D0': wD0, 'D1': wD1, 'Dsq2': wDsq2, 'D2': wD2, 'Dsq5': wDsq5}


def generate_hilbert_hier_bondfile(n0, e_backbone, e_red, e_green, e_blue, path):
    """Generate a bond file with 3-level hierarchical Hilbert-curve couplings.

    Every pair (i, j) receives a coupling value determined by where in the Hilbert
    chain traversal order they sit:
      - RED  (level 0): same group of 4 chain positions   → e_red
      - GREEN(level 1): same top-level quarter, different  → e_green
      - BLUE (level 2): different top-level quarters       → e_blue

    Coupling is placed in both wD1 (d=1) and wDsq2 (d=√2) matrices so that it
    fires at any short-range physical contact.  Backbone confinement (wD2, wDsq5)
    is added as in the Hilbert-local mode.

    Returns (matrices_dict, bond_colors_dict).  bond_colors_dict maps (id_i, id_j)
    → colour string for Hilbert-grid-adjacent non-backbone pairs (used for animation).
    Backbone pairs are NOT included in bond_colors_dict; the caller colours them black.
    """
    l0 = round(math.sqrt(n0))
    if not _is_power_of_2(l0):
        raise ValueError(f"Hilbert-hier mode requires l0=sqrt(n0) to be a power of 2; "
                         f"got n0={n0}, l0={l0}")

    backbone_bonds = snake_path(n0)
    backbone_set = set()
    for (pi, pj) in backbone_bonds:
        backbone_set.add((pi, pj))
        backbone_set.add((pj, pi))

    chain_order = _chain_order(n0)
    chain_pos = {pid: k for k, pid in enumerate(chain_order)}

    # Hilbert grid positions for each particle (for animation adjacency)
    sfc_pos = _space_filling_positions(l0)
    hpos = {}
    for k, pid in enumerate(chain_order):
        hpos[pid] = sfc_pos[k]

    level_energies = [e_red, e_green, e_blue]

    wD1   = np.zeros((n0, n0))
    wDsq2 = np.zeros((n0, n0))
    wD2   = np.zeros((n0, n0))
    wDsq5 = np.zeros((n0, n0))
    wD0   = np.zeros((n0, n0))

    # Assign hierarchical coupling only to Hilbert-grid-adjacent non-backbone pairs,
    # mirroring the logic of generate_hilbert_bondfile but with level-dependent values:
    #   cardinal Hilbert neighbours (d=1 in Hilbert grid)  -> wD1   only
    #   diagonal Hilbert neighbours (d=sqrt2 in Hilbert grid) -> wDsq2 only
    #
    # Placing each coupling at its native Hilbert distance means stretching a cardinal
    # bond from d=1 to d=sqrt2 costs +val (DeltaU > 0) -> VMMC recruits the neighbour.
    # Likewise, stretching a diagonal bond from d=sqrt2 to d=2 costs +val -> recruited.
    # This keeps each group rigid while allowing inter-group conformational changes.
    # Backbone pairs are excluded (confined by wD2/wDsq5; freely articulated at d=1/sqrt2).
    n_d1, n_dsq2 = 0, 0
    for i in range(n0):
        for j in range(i + 1, n0):
            if (i, j) in backbone_set:
                continue
            dx = hpos[j][0] - hpos[i][0]
            dy = hpos[j][1] - hpos[i][1]
            d2 = dx * dx + dy * dy
            if d2 > 2:
                continue  # not Hilbert-adjacent; no coupling
            cp_i = chain_pos[i]
            cp_j = chain_pos[j]
            val = level_energies[_hilbert_bond_level(cp_i, cp_j, n0)]
            if d2 == 1:
                wD1[i, j] = wD1[j, i] = val;    n_d1 += 1
            else:  # d2 == 2
                wDsq2[i, j] = wDsq2[j, i] = val; n_dsq2 += 1

    # Backbone confinement: repulsion at d=2 and d=√5
    for (pi, pj) in backbone_bonds:
        wD2[pi, pj]   -= e_backbone;  wD2[pj, pi]   -= e_backbone
        wDsq5[pi, pj] -= e_backbone;  wDsq5[pj, pi] -= e_backbone
        wD0[pi, pj]   -= e_backbone;  wD0[pj, pi]   -= e_backbone

    def write_matrix(f, mat, tag):
        f.write(f"{tag}\n")
        for i in range(n0):
            for j in range(i, n0):
                f.write(f"{i} {j} {mat[i,j]:.6f}\n")
        f.write(f"{tag}_END\n\n")

    with open(path, 'w') as f:
        f.write(f"# Hierarchical Hilbert bond file  n0={n0}  e_backbone={e_backbone}"
                f"  e_red={e_red}  e_green={e_green}  e_blue={e_blue}\n\n")
        f.write("BACKBONE\n")
        f.write("BACKBONE_END\n\n")
        write_matrix(f, wD1,   "WEAK_D1")
        write_matrix(f, wDsq2, "WEAK_DSQRT2")
        write_matrix(f, wD2,   "WEAK_D2")
        write_matrix(f, wDsq5, "WEAK_DSQRT5")

    print(f"Hierarchical Hilbert bond file written to {path}  "
          f"(e_red={e_red}, e_green={e_green}, e_blue={e_blue}, "
          f"{n_d1} cardinal pairs, {n_dsq2} diagonal pairs)")

    # Build bond_colors for Hilbert-grid-adjacent non-backbone pairs (animation only)
    bond_colors = {}
    for i in range(n0):
        for j in range(i + 1, n0):
            if (i, j) in backbone_set:
                continue
            dx = hpos[j][0] - hpos[i][0]
            dy = hpos[j][1] - hpos[i][1]
            if dx * dx + dy * dy <= 2:  # Hilbert-adjacent (d=1 or d=√2 in Hilbert grid)
                level = _hilbert_bond_level(chain_pos[i], chain_pos[j], n0)
                col = _HIER_COLOURS[level]
                bond_colors[(i, j)] = col
                bond_colors[(j, i)] = col

    matrices = {'D0': wD0, 'D1': wD1, 'Dsq2': wDsq2, 'D2': wD2, 'Dsq5': wDsq5}
    return matrices, bond_colors


def generate_denatured_bondfile(n0, e_backbone, path):
    """Generate a bond file with only backbone confinement — no weak coupling.

    All weak matrices (D1, Dsq2) are zero.  The D2 and Dsq5 matrices retain
    only the −e_backbone confinement entries for backbone-bonded pairs, so the
    chain topology is maintained while the polymer freely explores conformations.
    Used for the denaturation (equilibration) phase before the main simulation.
    """
    backbone_bonds = snake_path(n0)
    wD2   = np.zeros((n0, n0))
    wDsq5 = np.zeros((n0, n0))
    wD0   = np.zeros((n0, n0))  # display-only

    for (pi, pj) in backbone_bonds:
        wD2[pi, pj]   -= e_backbone;  wD2[pj, pi]   -= e_backbone
        wDsq5[pi, pj] -= e_backbone;  wDsq5[pj, pi] -= e_backbone
        wD0[pi, pj]   -= e_backbone;  wD0[pj, pi]   -= e_backbone

    def write_matrix(f, mat, tag):
        f.write(f"{tag}\n")
        for i in range(n0):
            for j in range(i, n0):
                f.write(f"{i} {j} {mat[i,j]:.6f}\n")
        f.write(f"{tag}_END\n\n")

    with open(path, 'w') as f:
        f.write(f"# Denatured bond file  n0={n0}  e_backbone={e_backbone}\n")
        f.write("# Weak couplings zeroed; only backbone confinement active\n\n")
        f.write("BACKBONE\n")
        f.write("BACKBONE_END\n\n")
        write_matrix(f, np.zeros((n0, n0)), "WEAK_D1")
        write_matrix(f, np.zeros((n0, n0)), "WEAK_DSQRT2")
        write_matrix(f, wD2,               "WEAK_D2")
        write_matrix(f, wDsq5,             "WEAK_DSQRT5")

    print(f"Denatured bond file written to {path}  (backbone confinement only, weak bonds zeroed)")
    return {'D0': wD0, 'D1': np.zeros((n0, n0)), 'Dsq2': np.zeros((n0, n0)),
            'D2': wD2, 'Dsq5': wDsq5}


def write_conf_from_last_frame(trajfile, confpath):
    """Read the last frame of a trajectory file and write it as an initial config."""
    _, _, _, frames = parse_traj(trajfile)
    if not frames:
        sys.exit(f"Error: no frames found in {trajfile}")
    last = frames[-1]
    with open(confpath, 'w') as f:
        for x, y in last:
            f.write(f"{int(round(x))} {int(round(y))}\n")
    print(f"Wrote denatured final config ({len(last)} particles) to {confpath}")


def _chain_order(n0):
    """Return particle indices in chain sequence order (endpoints first, path traversal)."""
    bonds = snake_path(n0)
    adj = {}
    for (a, b) in bonds:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    # Find an endpoint (degree 1)
    start = next(k for k, v in adj.items() if len(v) == 1)
    order, prev, cur = [start], None, start
    while len(order) < n0:
        nexts = [x for x in adj.get(cur, []) if x != prev]
        if not nexts:
            break
        nxt = nexts[0]
        order.append(nxt)
        prev, cur = cur, nxt
    return order


def _is_power_of_2(n):
    """Return True if n is a positive power of 2."""
    return n > 0 and (n & (n - 1)) == 0


def _hilbert_bond_level(cp_i, cp_j, n0):
    """Return hierarchical level for a pair of chain positions.

    0 = RED  : same group of 4 consecutive Hilbert-chain positions
    1 = GREEN : same top-level quarter (n0//4 positions) but different group of 4
    2 = BLUE  : different top-level quarters
    """
    top_size = max(n0 // 4, 4)
    if cp_i // 4 == cp_j // 4:
        return 0
    elif cp_i // top_size == cp_j // top_size:
        return 1
    else:
        return 2


def _hilbert_d2xy(n, d):
    """Convert Hilbert curve index d to (x, y) in n×n grid (n must be a power of 2).
    Consecutive indices always map to positions at Manhattan distance 1 (no overlaps).
    """
    x = y = 0
    s = 1
    t = d
    while s < n:
        rx = 1 & (t // 2)
        ry = 1 & (t ^ rx)
        if ry == 0:
            if rx == 1:
                x = s - 1 - x
                y = s - 1 - y
            x, y = y, x
        x += s * rx
        y += s * ry
        t //= 4
        s *= 2
    return x, y


def _space_filling_positions(l0):
    """Return (x, y) positions for l0×l0 grid in space-filling curve order.

    For power-of-2 l0: uses a Hilbert curve (Moore-like space-filling curve) so that
    consecutive positions are always at Manhattan distance 1 — all backbone bonds are
    satisfied immediately and particles fill the grid without overlap.
    For other l0: uses a boustrophedon (snake) path with the same adjacency property.
    """
    if _is_power_of_2(l0):
        return [_hilbert_d2xy(l0, d) for d in range(l0 * l0)]
    else:
        positions = []
        for row in range(l0):
            cols = range(l0) if row % 2 == 0 else range(l0 - 1, -1, -1)
            for col in cols:
                positions.append((col, row))
        return positions


def write_conf_file(path, n0, box_length):
    """Write initial polymer config using a space-filling curve.

    For power-of-2 l0 (= sqrt(n0)): uses a Hilbert space-filling curve so that
    particles consecutive in chain order are at Manhattan distance 1 (backbone bonds
    immediately satisfied) and the whole chain fits in an l0×l0 footprint with no
    overlaps.  For other l0 the same is achieved with a snake path.
    The footprint is centred in the simulation box.
    """
    L = int(round(box_length))
    l0 = round(math.sqrt(n0))
    order = _chain_order(n0)

    sfc_pos = _space_filling_positions(l0)  # n0 grid positions in curve order

    # Centre the l0×l0 footprint in the box
    x_offset = (L - l0) // 2
    y_offset = (L - l0) // 2

    # Map particle ID → grid position: order[k] is placed at sfc_pos[k]
    pos = {}
    for chain_idx, pid in enumerate(order):
        sx, sy = sfc_pos[chain_idx]
        pos[pid] = (sx + x_offset, sy + y_offset)

    with open(path, 'w') as f:
        for pid in range(n0):
            f.write(f"{pos[pid][0]} {pos[pid][1]}\n")
    curve_name = "Hilbert" if _is_power_of_2(l0) else "snake"
    print(f"Initial config written to {path}  ({n0} particles, {curve_name} space-filling curve)")


def run_polymer_sim(exe, input_file, bond_file, conf_file=None, seed=None):
    cmd = [exe, input_file, bond_file]
    if conf_file:
        cmd.append(conf_file)
    if seed is not None:
        cmd.append(str(seed))
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f"Simulation failed (return code {result.returncode})")
    print("Simulation complete.\n")


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

# Hierarchical bond level colours: RED (finest), GREEN (mid), BLUE (coarsest)
_HIER_COLOURS = ['#d62728', '#2ca02c', '#1f77b4']
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
               frames, filehead, e1, custom_bonds=None, coupling_matrices=None,
               bond_colors=None):
    """Animated simulation viewer.

    If coupling_matrices is provided (polymer mode) the right-hand panel shows the
    five coupling matrices (D0, D1, Dsq2, D2, Dsq5) as physical-energy heatmaps
    (blue = attractive, red = repulsive).  Otherwise the classic bond-strength matrix
    is shown for hierarchical-assembly mode.
    """
    bond_table = compute_bond_table(n0, e1, custom_bonds)

    colours = _slot_colours(n0)
    n_frames = len(frames)
    L = float(box_length)

    # Colormap for bond lines in the lattice panel
    bond_vals = sorted(set(bond_table.values()))
    bond_cmap = plt.cm.plasma
    bond_norm = mcolors.Normalize(vmin=bond_vals[0] if bond_vals else 0, vmax=e1)

    # ---- Figure layout -------------------------------------------------------
    if coupling_matrices is not None:
        # Polymer mode: lattice | energy (top-right) + 5 coupling matrices (bottom-right)
        _CM_KEYS   = ['D0',      'D1',   'Dsq2',    'D2',   'Dsq5']
        _CM_TITLES = ['d=0\n(hard sphere)', 'd=1', 'd=√2', 'd=2', 'd=√5']
        fig = plt.figure(figsize=(24, 7))
        gs = fig.add_gridspec(2, 7,
                              width_ratios=[2.2, 1, 1, 1, 1, 1, 0.15],
                              hspace=0.55, wspace=0.45)
        ax_lat = fig.add_subplot(gs[:, 0])
        ax_en  = fig.add_subplot(gs[0, 1:6])
        mat_axes = [fig.add_subplot(gs[1, c]) for c in range(1, 6)]

        for ax, key, title in zip(mat_axes, _CM_KEYS, _CM_TITLES):
            mat = coupling_matrices.get(key)
            if mat is None:
                ax.axis('off')
                continue
            phys = -mat  # physical energy: negative = attractive
            vmax = np.abs(phys).max() or 1.0
            im = ax.imshow(phys, cmap='RdBu', vmin=-vmax, vmax=vmax,
                           origin='upper', aspect='auto', interpolation='nearest')
            ax.set_title(title, fontsize=7, pad=3)
            ax.set_xlabel("id", fontsize=6)
            ax.set_ylabel("id", fontsize=6)
            tick_step = max(1, n0 // 4)
            ticks = list(range(0, n0, tick_step))
            ax.set_xticks(ticks); ax.set_yticks(ticks)
            ax.tick_params(labelsize=5)
            fig.colorbar(im, ax=ax, shrink=0.85, pad=0.04,
                         label='E phys\n(kBT)').ax.tick_params(labelsize=5)

    else:
        # Classic mode: lattice | energy | bond-strength matrix
        coupling_matrix = np.zeros((n0, n0))
        for (i, j), val in bond_table.items():
            coupling_matrix[i, j] = val

        fig = plt.figure(figsize=(15, 7))
        gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1], hspace=0.45, wspace=0.35)
        ax_lat = fig.add_subplot(gs[:, 0])
        ax_en  = fig.add_subplot(gs[0, 1])
        ax_cm  = fig.add_subplot(gs[1, 1])

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

    fig.suptitle(
        f"HierarchicalAssembly  |  n0={n0}  nParticles={n_particles}"
        f"  e1={e1}  dens={n_particles / L**2:.3f}",
        fontsize=11,
    )

    # ---- Lattice panel -------------------------------------------------------
    ax_lat.set_facecolor("white" if coupling_matrices is not None else "#111111")
    ax_lat.set_xlim(0, L)
    ax_lat.set_ylim(0, L)
    ax_lat.set_aspect("equal")
    ax_lat.set_xlabel("x")
    ax_lat.set_ylabel("y")
    lat_title = ax_lat.set_title("Step 0", fontsize=10)

    bond_lc = LineCollection([], linewidths=2.5, alpha=0.85, zorder=1)
    ax_lat.add_collection(bond_lc)

    circles = []
    for i, (x, y) in enumerate(frames[0]):
        circ = patches.Circle(
            (x + 0.5, y + 0.5), radius=0.35,
            facecolor=colours[i % n0],
            edgecolor="black" if coupling_matrices is not None else "white",
            linewidth=0.5, zorder=2,
        )
        ax_lat.add_patch(circ)
        circles.append(circ)

    # ---- Energy panel --------------------------------------------------------
    ax_en.plot(steps, energy, lw=1, color="steelblue", alpha=0.7, label="Energy")
    running_avg = np.cumsum(energy) / (np.arange(len(energy)) + 1)
    ax_en.plot(steps, running_avg, lw=1.5, color="tomato", label="Running avg")
    ax_en.set_xlabel("Step")
    ax_en.set_ylabel("Total energy")
    ax_en.set_title("Energy vs time")
    ax_en.legend(fontsize=8)
    ax_en.grid(True, alpha=0.3)
    vline = ax_en.axvline(steps[0], color="gold", lw=1.0, linestyle="--", alpha=0.9)

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
                    if bond_colors is not None:
                        seg_colors.append(bond_colors.get((id1, j % n0), 'black'))
                    elif coupling_matrices is not None:
                        seg_colors.append('black')
                    else:
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
# Boltzmann validation
# ---------------------------------------------------------------------------

# --- Symmetry helpers for rotation+reflection invariance ---

def _rot90(v):
    """Rotate 2D vector 90° counter-clockwise."""
    return (-v[1], v[0])

def _reflx(v):
    """Reflect 2D vector across the x-axis."""
    return (v[0], -v[1])

def _symmetry_group(positions):
    """Return all ≤8 distinct dihedral transforms of a tuple of (x,y) pairs."""
    versions = set()
    cur = positions
    for _ in range(4):
        versions.add(cur)
        versions.add(tuple(_reflx(p) for p in cur))
        cur = tuple(_rot90(p) for p in cur)
    return versions

def _canonical_key(rel_positions):
    """Canonical (lexicographically minimum) form under translation+rotation+reflection."""
    return min(_symmetry_group(rel_positions))

def _conformation_degeneracy(canonical_rel):
    """Orbit size: number of distinct labeled conformations related by rotation/reflection."""
    return len(_symmetry_group(canonical_rel))


def conformation_key(coords_copy, L):
    """Conformation key invariant under translation, rotation, and reflection.
    coords_copy: numpy array shape (n0, 2).
    """
    r0 = coords_copy[0]
    rel = []
    for i in range(1, len(coords_copy)):
        d = coords_copy[i] - r0
        d = d - L * np.round(d / L)  # minimum image
        rel.append((round(d[0]), round(d[1])))
    return _canonical_key(tuple(rel))


def compute_conformation_energy(key, n0, coupling_matrices):
    """Sum of weak coupling values for all pairs in the conformation.

    Sign convention: positive value = attractive (C++ returns -value as pair energy).
    Boltzmann weight = exp(+energy) [NOT exp(-energy)].

    key: canonical relative position tuple for particles 1..n0-1 (relative to particle 0).
    """
    positions = [(0, 0)]  # particle 0 at origin
    for dx, dy in key:
        positions.append((dx, dy))

    wD1   = coupling_matrices.get('D1')
    wDsq2 = coupling_matrices.get('Dsq2')
    wD2   = coupling_matrices.get('D2')
    wDsq5 = coupling_matrices.get('Dsq5')

    energy = 0.0
    for i in range(n0):
        for j in range(i+1, n0):
            dx = positions[j][0] - positions[i][0]
            dy = positions[j][1] - positions[i][1]
            d2r = round(dx*dx + dy*dy)
            if d2r == 1 and wD1 is not None:
                energy += wD1[i, j]
            elif d2r == 2 and wDsq2 is not None:
                energy += wDsq2[i, j]
            elif d2r == 4 and wD2 is not None:
                energy += wD2[i, j]
            elif d2r == 5 and wDsq5 is not None:
                energy += wDsq5[i, j]
    return energy


def enumerate_canonical_states(n0):
    """Enumerate all canonically distinct conformations for the n0-bead polymer chain.

    Returns list of (canonical_key, example_rel_positions, degeneracy) tuples,
    sorted by canonical_key.  example_rel_positions gives one spatial realisation.

    Conformations are self-avoiding walks where each step is one of the 8 nearest
    lattice directions (cardinal + diagonal), following the snake-path chain sequence.
    """
    chain_seq = _chain_order(n0)  # e.g. [0,1,3,2] for n0=4
    steps = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)]
    canonical_set = {}  # key -> example rel_positions

    def dfs(geom_positions, occupied):
        if len(geom_positions) == n0:
            # geom_positions[k] = spatial position of chain_seq[k]-th particle
            pos_by_pid = {chain_seq[k]: geom_positions[k] for k in range(n0)}
            p0 = pos_by_pid[0]
            rel = tuple((pos_by_pid[pid][0]-p0[0], pos_by_pid[pid][1]-p0[1])
                        for pid in range(1, n0))
            key = _canonical_key(rel)
            if key not in canonical_set:
                canonical_set[key] = rel
            return
        last = geom_positions[-1]
        for dx, dy in steps:
            nxt = (last[0]+dx, last[1]+dy)
            if nxt not in occupied:
                dfs(geom_positions + [nxt], occupied | {nxt})

    dfs([(0,0)], {(0,0)})

    result = [(k, v, _conformation_degeneracy(k)) for k, v in canonical_set.items()]
    result.sort(key=lambda x: x[0])
    return result


def boltzmann_validation(frames, n0, L, coupling_matrices, steps_array):
    """Compute Pearson correlation between observed and Boltzmann frequencies vs step count.

    Sign convention: C++ stores +coupling as attractive (returns -coupling as pair energy),
    so Boltzmann weight ∝ degeneracy × exp(+E_python) where E_python = sum of coupling values.

    Returns:
        step_indices, correlations, obs_freq, boltz_freq, E_vals
    """
    import collections

    nframes = len(frames)
    counts = collections.Counter()
    step_indices = []
    correlations = []

    check_at = set(np.unique(np.round(
        np.logspace(0, np.log10(max(nframes-1, 1)), 50)).astype(int)))
    check_at.add(nframes - 1)

    for fi, frame in enumerate(frames):
        coords_copy = frame[:n0]
        key = conformation_key(coords_copy, L)
        counts[key] += 1

        if fi in check_at and len(counts) >= 2:
            energies = {k: compute_conformation_energy(k, n0, coupling_matrices) for k in counts}
            E_vals  = np.array([energies[k] for k in counts])
            degens  = np.array([_conformation_degeneracy(k) for k in counts], dtype=float)
            obs_counts = np.array([counts[k] for k in counts], dtype=float)

            # Boltzmann weight ∝ g × exp(+E_python)  [positive coupling = attractive in C++]
            boltz = degens * np.exp(E_vals)
            boltz /= boltz.sum()
            obs_freq = obs_counts / obs_counts.sum()

            if len(obs_freq) >= 2 and np.std(obs_freq) > 0 and np.std(boltz) > 0:
                r = float(np.corrcoef(obs_freq, boltz)[0, 1])
                step_indices.append(fi)
                correlations.append(r)

    # Final state for scatter plot
    energies   = {k: compute_conformation_energy(k, n0, coupling_matrices) for k in counts}
    E_vals     = np.array([energies[k] for k in counts])
    degens     = np.array([_conformation_degeneracy(k) for k in counts], dtype=float)
    obs_counts = np.array([counts[k] for k in counts], dtype=float)
    boltz      = degens * np.exp(E_vals)
    boltz     /= boltz.sum()
    obs_freq   = obs_counts / obs_counts.sum()

    return (np.array(step_indices), np.array(correlations), obs_freq, boltz, E_vals)


def plot_boltzmann_validation(step_indices, correlations, obs_freq, boltz_freq, state_energies, n0):
    """Plot Boltzmann validation: Pearson r vs steps + observed vs predicted scatter."""
    if len(step_indices) < 2:
        print("Not enough data for Boltzmann validation plot.")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Boltzmann Validation  |  n0={n0}", fontsize=12)

    # Left: Pearson correlation vs number of simulation steps
    ax1.plot(step_indices, correlations, 'o-', lw=1.5, ms=4, color='steelblue')
    ax1.axhline(1.0, color='tomato', ls='--', lw=1, label='Perfect correlation')
    ax1.set_xlabel("Simulation frame index")
    ax1.set_ylabel("Pearson r (observed vs Boltzmann)")
    ax1.set_title("Convergence to Boltzmann distribution")
    ax1.set_ylim(-0.1, 1.1)
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    if len(step_indices) > 1:
        ax1.set_xscale('log')

    # Right: Observed vs predicted frequencies (final)
    # Colour by physical energy (negative = attractive = blue)
    phys_energies = -state_energies
    e_norm = plt.Normalize(phys_energies.min(), phys_energies.max())
    ax2.scatter(boltz_freq, obs_freq, c=phys_energies, cmap='RdBu', norm=e_norm,
                s=40, alpha=0.8, zorder=3)
    lims = [0, max(boltz_freq.max(), obs_freq.max()) * 1.1]
    ax2.plot(lims, lims, 'r--', lw=1, label='y=x (perfect)')
    ax2.set_xlabel("Boltzmann predicted frequency")
    ax2.set_ylabel("Observed frequency")
    ax2.set_title(f"Observed vs Predicted  (final r={correlations[-1]:.3f})")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    sm = plt.cm.ScalarMappable(cmap='RdBu', norm=e_norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax2, label='Physical energy (kBT)\n(blue = attractive)')
    plt.tight_layout()
    plt.show()


def plot_weak_coupling_matrices(matrices, n0):
    """Show the 4 weak coupling matrices in a 2x2 grid."""
    labels = [('D1', 'Distance 1'), ('Dsq2', 'Distance sqrt(2)'), ('D2', 'Distance 2'), ('Dsq5', 'Distance sqrt(5)')]
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    fig.suptitle(f"Weak Coupling Matrices  |  n0={n0}", fontsize=12)
    for ax, (key, title) in zip(axes.flat, labels):
        mat = matrices.get(key)
        if mat is None:
            ax.set_visible(False)
            continue
        # Display as physical energy (negative = attractive): phys = −coupling
        phys = -mat
        vmax = np.abs(phys).max() or 1.0
        im = ax.imshow(phys, cmap='RdBu', vmin=-vmax, vmax=vmax,
                       origin='upper', aspect='auto', interpolation='nearest')
        fig.colorbar(im, ax=ax, label='Physical energy (kBT)\n(blue = attractive)')
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Particle identity")
        ax.set_ylabel("Particle identity")
        tick_step = max(1, n0 // 8)
        ticks = list(range(0, n0, tick_step))
        ax.set_xticks(ticks); ax.set_yticks(ticks)
    plt.tight_layout()
    plt.show()


def plot_all_states(n0, coupling_matrices=None):
    """Plot every canonically distinct conformation of the n0-bead polymer chain.

    Each state is shown as a small diagram: coloured circles (one per particle)
    connected by lines along the chain backbone.  The title shows the state number,
    degeneracy g (number of rotation/reflection copies), and, if coupling matrices
    are provided, the weak-coupling energy.

    Linear (palindromic) conformations have g=4; general shapes have g=8.
    The Boltzmann weight of each state ∝ g × exp(+E).
    """
    print("Enumerating canonical states...")
    states = enumerate_canonical_states(n0)
    n_states = len(states)
    print(f"  Found {n_states} canonical conformations.")

    if coupling_matrices is not None:
        # Sort by physical energy ascending (most negative = most favourable first)
        states = sorted(
            states,
            key=lambda s: -compute_conformation_energy(s[0], n0, coupling_matrices),
        )

    chain_seq = _chain_order(n0)          # e.g. [0,1,3,2] for n0=4
    colours = _slot_colours(n0)

    ncols = 7
    nrows = math.ceil(n_states / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(2.2 * ncols, 2.4 * nrows))
    fig.suptitle(
        f"All canonical conformations  |  n0={n0}  "
        f"({n_states} states, g=4 linear / g=8 general)",
        fontsize=11,
    )

    for idx, (key, example_rel, deg) in enumerate(states):
        ax = axes.flat[idx]

        # Build {particle_id: (x, y)} from example_rel
        positions = {0: (0, 0)}
        for pid, (dx, dy) in zip(range(1, n0), example_rel):
            positions[pid] = (dx, dy)

        # Bond lines along chain sequence
        for k in range(len(chain_seq) - 1):
            p1, p2 = chain_seq[k], chain_seq[k + 1]
            x1, y1 = positions[p1]
            x2, y2 = positions[p2]
            ax.plot([x1, x2], [y1, y2], 'k-', lw=2, zorder=1)

        # Particle circles
        for pid in range(n0):
            x, y = positions[pid]
            ax.scatter(x, y, s=180,
                       c=[colours[pid]], edgecolors='k',
                       linewidth=0.6, zorder=2)

        # Title — show physical energy (negative = favourable)
        if coupling_matrices is not None:
            E_phys = -compute_conformation_energy(key, n0, coupling_matrices)
            title = f"{idx+1}. g={deg}  E={E_phys:.2f}"
        else:
            title = f"{idx+1}. g={deg}"
        ax.set_title(title, fontsize=6.5, pad=2)

        # Axes limits
        xs = [positions[p][0] for p in range(n0)]
        ys = [positions[p][1] for p in range(n0)]
        m = 0.9
        ax.set_xlim(min(xs) - m, max(xs) + m)
        ax.set_ylim(min(ys) - m, max(ys) + m)
        ax.set_aspect('equal')
        ax.axis('off')

    # Hide empty cells
    for idx in range(n_states, nrows * ncols):
        axes.flat[idx].axis('off')

    plt.tight_layout()
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

    # --- Polymer mode ---
    polymer_matrices = None
    hier_bond_colors = None
    exe_polymer = os.path.join(script_dir, "run_polymer")

    if args.polymer is not None:
        # Override n0 and ncopies
        args.n0 = args.polymer
        ncopies = 1  # one polymer copy for Boltzmann validation
        n_particles = args.n0

        # Recompute dens if L given
        if args.L is not None:
            dens = n_particles / args.L ** 2

        # Generate bond file
        bond_file = os.path.join(script_dir, f"{args.filehead}_polymer_bonds.txt")
        if args.hier_red is not None or args.hier_green is not None or args.hier_blue is not None:
            e_red   = args.hier_red   if args.hier_red   is not None else 0.0
            e_green = args.hier_green if args.hier_green is not None else 0.0
            e_blue  = args.hier_blue  if args.hier_blue  is not None else 0.0
            polymer_matrices, hier_bond_colors = generate_hilbert_hier_bondfile(
                args.n0, args.e1, e_red, e_green, e_blue, bond_file)
        elif args.hilbert_A is not None or args.hilbert_B is not None:
            A = args.hilbert_A if args.hilbert_A is not None else 0.0
            B = args.hilbert_B if args.hilbert_B is not None else 0.0
            polymer_matrices = generate_hilbert_bondfile(args.n0, args.e1, A, B, bond_file)
        else:
            polymer_matrices = generate_polymer_bondfile(
                args.n0, args.e1, args.weak_e, args.weak_std,
                bond_file, seed=args.weak_seed
            )

        # Write initial config
        box_len_for_conf = args.L if args.L else round(math.sqrt(n_particles / dens))
        conf_file = os.path.join(script_dir, f"{args.filehead}_init.conf")
        write_conf_file(conf_file, args.n0, box_len_for_conf)

        if not args.no_run:
            write_input_file(input_file, args.filehead, args.n0, ncopies,
                             args.nsteps, args.nsweep, dens)

            if args.denature is not None:
                # Phase 1: denaturation — run with weak bonds zeroed
                denature_head = args.filehead + "_denature"
                denature_bond_file = os.path.join(script_dir, f"{denature_head}_bonds.txt")
                denature_input     = os.path.join(script_dir, f"input_{denature_head}.txt")
                denature_conf      = os.path.join(script_dir, f"{denature_head}_final.conf")
                denature_traj      = os.path.join(script_dir, f"{denature_head}_traj.txt")

                generate_denatured_bondfile(args.n0, args.e1, denature_bond_file)
                write_input_file(denature_input, denature_head, args.n0, ncopies,
                                 args.denature, args.nsweep, dens)

                print(f"\n--- Denaturation phase ({args.denature} steps, weak bonds off) ---")
                run_polymer_sim(exe_polymer, denature_input, denature_bond_file,
                                conf_file, args.sim_seed)

                # Extract final config from denaturation run for use as main-phase seed
                write_conf_from_last_frame(denature_traj, denature_conf)

                print(f"\n--- Main simulation phase ({args.nsteps} steps) ---")
                run_polymer_sim(exe_polymer, input_file, bond_file, denature_conf, args.sim_seed)
            else:
                run_polymer_sim(exe_polymer, input_file, bond_file, conf_file, args.sim_seed)

        print("Parsing output files...")
        steps, energy, fragment_hist = parse_stats(statsfile)
        n_particles, box_length, n0, frames = parse_traj(trajfile)
        print(f"  {len(frames)} frames  |  {n_particles} particles  |  box={box_length}  |  n0={n0}")

        # Build backbone bond table for animation (snake-path bonds at energy e1)
        backbone_dict = {}
        for (pi, pj) in snake_path(n0):
            backbone_dict[(pi, pj)] = args.e1
            backbone_dict[(pj, pi)] = args.e1

        anim_bond_colors = None
        if hier_bond_colors is not None:
            # Add Hilbert-adjacent non-backbone bonds to the draw table
            for (pi, pj) in hier_bond_colors:
                if (pi, pj) not in backbone_dict:
                    backbone_dict[(pi, pj)] = 1.0
            # Per-bond colour: backbone → black, others → hierarchical colour
            anim_bond_colors = {
                (pi, pj): hier_bond_colors.get((pi, pj), 'black')
                for (pi, pj) in backbone_dict
            }

        # Animation window: lattice + energy + all 5 coupling matrices
        make_plots(steps, energy, fragment_hist, n_particles, box_length, n0,
                   frames, args.filehead, args.e1,
                   custom_bonds=backbone_dict, coupling_matrices=polymer_matrices,
                   bond_colors=anim_bond_colors)

        if args.boltzmann:
            # Show all canonical conformations ranked by physical energy
            plot_all_states(n0, coupling_matrices=polymer_matrices)

            # Run Boltzmann validation
            print("Running Boltzmann validation...")
            step_indices, correlations, obs_freq, boltz_freq, state_energies = boltzmann_validation(
                frames, n0, float(box_length), polymer_matrices, steps
            )
            plot_boltzmann_validation(step_indices, correlations, obs_freq, boltz_freq, state_energies, n0)

    else:
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
