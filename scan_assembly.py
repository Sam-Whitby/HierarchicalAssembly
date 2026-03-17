#!/usr/bin/env python3
"""Scan RED bond energy and plot thermodynamic/kinetic yield curves.

Demonstrates the kinetic-vs-thermodynamic trapping discrepancy for the first
hierarchical level (RED bonds, coupling particles within the same group of 4
consecutive Hilbert-chain positions).

Each simulation starts from a denatured (randomly explored) polymer conformation.
For each RED energy value the simulation runs to the longest requested time, then
yield is sampled at each of the requested time-points.  Multiple independent
repeats are averaged, with standard-error error bars.

YIELD DEFINITION
    For n0=64 the native Hilbert structure has 16 "closing bonds" — one per
    group of 4 chain positions.  Each closing bond is the non-backbone pair
    within the 2×2 sub-square that closes the U-shaped Hilbert sub-path.
    Yield = fraction of those 16 bonds that are physically formed (d ≤ √2 in
    the simulation box, accounting for periodic boundaries).

Usage
-----
    python scan_assembly.py scan_config.ini          # run and plot
    python scan_assembly.py scan_config.ini --plot-only  # replot from CSV

Config file format (INI)
------------------------
    [simulation]
    polymer       = 64      # n0: particles (sqrt must be a power-of-2)
    L             = 20      # box side length
    nsweep        = 1       # MC sweeps per output step
    hier_green    = 0.0     # GREEN bond energy (0 to isolate RED level)
    hier_blue     = 0.0     # BLUE bond energy  (0 to isolate RED level)
    spring_k      = 1.0     # backbone spring stiffness (None for hard-wall)
    e1            = 8.0     # backbone energy for hard-wall denaturation phase
    sim_seed      = 42      # base RNG seed (repeat r uses seed + r)
    unspecific    = true    # non-specific bonding mode
    denature_steps = 500    # equilibration steps with weak bonds off

    [scan]
    hier_red_values = 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0
    nsteps_curves   = 200, 1000, 5000   # time points; max is the sim length

    [output]
    output_dir = scan_results
    n_repeats  = 3
"""

import sys
import os
import math
import argparse
import configparser
import subprocess
import contextlib
import io

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---------------------------------------------------------------------------
# Locate the run_and_plot module (same directory as this script)
# ---------------------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

try:
    from run_and_plot import (
        generate_hilbert_hier_bondfile,
        write_conf_file,
        write_conf_from_last_frame,
        write_input_file,
        parse_traj,
        snake_path,
        _chain_order,
        _space_filling_positions,
    )
except ImportError as exc:
    sys.exit(f"ERROR: Could not import from run_and_plot.py: {exc}\n"
             f"Make sure this script is in the same directory as run_and_plot.py.")

_EXE = os.path.join(_SCRIPT_DIR, "run_polymer")


# ---------------------------------------------------------------------------
# Quiet helpers — suppress all output from C++ and Python library calls
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def _silence():
    """Redirect stdout to /dev/null for both Python prints and C++ subprocesses."""
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def _run_quiet(exe, input_file, bond_file, conf_file=None, seed=None, n_neighbours=4):
    """Run the C++ polymer executable with all output suppressed."""
    cmd = [exe, input_file, bond_file]
    if conf_file:
        cmd.append(conf_file)
    if seed is not None:
        cmd.append(str(seed))
    cmd += ['--neighbours', str(n_neighbours)]
    result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        raise RuntimeError(f"Simulation failed (exit code {result.returncode}): {' '.join(cmd)}")


# ---------------------------------------------------------------------------
# Progress bar
# ---------------------------------------------------------------------------

def _progress(done, total, phase, e_red, rep, n_repeats):
    """Overwrite the current terminal line with a compact progress bar."""
    bar_width = 22
    filled = int(bar_width * done / total)
    bar = '█' * filled + '░' * (bar_width - filled)
    line = (f"\r[{bar}]  {done}/{total}  "
            f"e_red={e_red:<5.2g}  rep={rep+1}/{n_repeats}  {phase:<12}")
    print(line, end='', flush=True)


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(path):
    if not os.path.exists(path):
        sys.exit(f"ERROR: config file not found: {path}")
    cfg = configparser.ConfigParser()
    cfg.read(path)

    def gf(sec, key, fallback):
        return float(cfg.get(sec, key, fallback=str(fallback)))

    def gi(sec, key, fallback):
        return int(cfg.get(sec, key, fallback=str(fallback)))

    def gb(sec, key, fallback):
        return cfg.get(sec, key, fallback=str(fallback)).strip().lower() in ('true', '1', 'yes')

    spring_k_str = cfg.get('simulation', 'spring_k', fallback='1.0').strip()
    spring_k = None if spring_k_str.lower() == 'none' else float(spring_k_str)

    sim_seed_str = cfg.get('simulation', 'sim_seed', fallback='None').strip()
    sim_seed = None if sim_seed_str.lower() == 'none' else int(sim_seed_str)

    red_str   = cfg.get('scan', 'hier_red_values')
    steps_str = cfg.get('scan', 'nsteps_curves')

    d = {
        'polymer':        gi('simulation', 'polymer', 64),
        'L':              gf('simulation', 'L', 20.0),
        'nsweep':         gi('simulation', 'nsweep', 1),
        'hier_green':     gf('simulation', 'hier_green', 0.0),
        'hier_blue':      gf('simulation', 'hier_blue', 0.0),
        'spring_k':       spring_k,
        'e1':             gf('simulation', 'e1', 8.0),
        'sim_seed':       sim_seed,
        'unspecific':     gb('simulation', 'unspecific', True),
        'denature_steps': gi('simulation', 'denature_steps', 500),
        'hier_red_values': [float(x.strip()) for x in red_str.split(',')],
        'nsteps_curves':   sorted([int(x.strip()) for x in steps_str.split(',')]),
        'output_dir':      cfg.get('output', 'output_dir', fallback='scan_results'),
        'n_repeats':       gi('output', 'n_repeats', 1),
        'n_neighbours':    gi('simulation', 'n_neighbours', 4),
    }

    if d['n_neighbours'] not in (4, 8):
        sys.exit(f"ERROR: n_neighbours must be 4 or 8, got {d['n_neighbours']}")

    return d


# ---------------------------------------------------------------------------
# Yield metric: native RED closing bonds
# ---------------------------------------------------------------------------

def find_native_red_closing_pairs(n0):
    """Return all native RED non-backbone contacts in the fully-folded Hilbert structure.

    In the native structure each group of 4 consecutive chain positions occupies a
    2×2 sub-square.  Three backbone bonds trace the U-shaped sub-path; the remaining
    non-backbone contacts within the group that are physically close in the native
    Hilbert grid are the 'native RED contacts':

      • 1 closing edge at grid d=1   (connects start and end of the U)
      • 2 diagonal contacts at d=√2  (across the diagonals of the 2×2 square)

    This gives 3 non-backbone native RED contacts per group of 4, i.e.
    3 × (n0 // 4) total pairs for a polymer of n0 particles.

    Yield = fraction of these pairs physically at d ≤ √2 in the simulation.
    In the fully folded native state: yield = 1.0.
    In a random denatured coil: yield ≈ 0.
    """
    l0 = round(math.sqrt(n0))
    chain_order = _chain_order(n0)                          # particle IDs in chain order
    chain_pos   = {pid: k for k, pid in enumerate(chain_order)}  # pid → position in chain
    sfc_pos     = _space_filling_positions(l0)              # grid (x,y) at each chain index
    native_pos  = {pid: sfc_pos[k] for k, pid in enumerate(chain_order)}

    bb_bonds = snake_path(n0)
    backbone_set = set()
    for (a, b) in bb_bonds:
        backbone_set.add((a, b))
        backbone_set.add((b, a))

    closing_pairs = []
    for i in range(n0):
        for j in range(i + 1, n0):
            # Must be in the same group of 4 (RED level)
            if chain_pos[i] // 4 != chain_pos[j] // 4:
                continue
            # Must NOT be a backbone bond
            if (i, j) in backbone_set:
                continue
            # Must be natively adjacent (d=1 or d=√2 in the Hilbert grid)
            dx = native_pos[j][0] - native_pos[i][0]
            dy = native_pos[j][1] - native_pos[i][1]
            if dx * dx + dy * dy <= 2:
                closing_pairs.append((i, j))

    return closing_pairs


def _pairs_formed(frame, pairs, L):
    """Return which pairs in `pairs` are physically at d ≤ √2 (PBC-aware)."""
    formed = []
    for (i, j) in pairs:
        dx = frame[j, 0] - frame[i, 0]
        dy = frame[j, 1] - frame[i, 1]
        dx -= L * round(dx / L)
        dy -= L * round(dy / L)
        if dx * dx + dy * dy <= 2.01:
            formed.append((i, j))
    return formed


def assembled_particles(frame, closing_pairs, L, n0):
    """Return the set of particle IDs belonging to a fully assembled 2×2 sub-group.

    A group of 4 is 'fully assembled' when all 3 of its native RED contacts are
    simultaneously formed.  Every particle in such a group is considered to be in
    the correct Hilbert position and is coloured green in the lattice plot.
    """
    chain_order = _chain_order(n0)
    chain_pos   = {pid: k for k, pid in enumerate(chain_order)}

    from collections import defaultdict
    group_pairs = defaultdict(list)
    for (i, j) in closing_pairs:
        group_pairs[chain_pos[i] // 4].append((i, j))

    green = set()
    for g, pairs in group_pairs.items():
        if len(_pairs_formed(frame, pairs, L)) == len(pairs):   # all formed
            for pid in chain_order[4 * g: 4 * g + 4]:
                green.add(pid)
    return green


def compute_yield(frame, closing_pairs, L):
    """Fraction of closing-bond pairs that are physically at d ≤ √2 (PBC-aware)."""
    if not closing_pairs:
        return 0.0
    return len(_pairs_formed(frame, closing_pairs, L)) / len(closing_pairs)


# ---------------------------------------------------------------------------
# Single-condition runner
# ---------------------------------------------------------------------------

def run_condition(cfg, e_red, work_dir, repeat=0, status_fn=None):
    """Run denaturation + main simulation for one (e_red, repeat) condition.

    Denaturation uses the SAME backbone physics as the main run (same spring_k)
    but with all weak bond energies set to zero, so the polymer freely explores
    conformations without any assembly bias.

    Returns a dict {nsteps: yield} sampled from the trajectory.
    status_fn(phase_str) is called before each phase to update the progress bar.
    """
    def _status(s):
        if status_fn:
            status_fn(s)

    os.makedirs(work_dir, exist_ok=True)

    n0           = cfg['polymer']
    L            = cfg['L']
    nsweep       = cfg['nsweep']
    max_t        = cfg['nsteps_curves'][-1]
    e_green      = cfg['hier_green']
    e_blue       = cfg['hier_blue']
    spring_k     = cfg['spring_k']
    e1           = cfg['e1']
    unspec       = cfg['unspecific']
    den_steps    = cfg['denature_steps']
    dens         = n0 / (L ** 2)
    n_neighbours = cfg['n_neighbours']

    seed = cfg['sim_seed']
    if seed is not None:
        seed = seed + repeat

    fh = os.path.join(work_dir, 'sim')

    bond_file     = fh + '_bonds.txt'
    den_bond_file = fh + '_den_bonds.txt'
    inp_file      = fh + '_input.txt'
    den_inp_file  = fh + '_den_input.txt'
    init_conf     = fh + '_init.conf'
    den_conf      = fh + '_den_final.conf'
    den_traj      = fh + '_den_traj.txt'
    traj_file     = fh + '_traj.txt'

    # 1. Bond files — both use the same backbone physics (same spring_k).
    #    Denaturation bond file has all weak energies = 0: the polymer diffuses
    #    freely connected but with no assembly driving force.
    with _silence():
        generate_hilbert_hier_bondfile(
            n0, e1, e_red, e_green, e_blue, bond_file,
            spring_k=spring_k, unspecific=unspec,
        )
        generate_hilbert_hier_bondfile(
            n0, e1, 0.0, 0.0, 0.0, den_bond_file,
            spring_k=spring_k, unspecific=unspec,
        )
        write_conf_file(init_conf, n0, L)

    # 2. Denaturation: zero weak bonds, same spring backbone
    _status("denature")
    write_input_file(den_inp_file, fh + '_den', n0, 1, den_steps, nsweep, dens)
    _run_quiet(_EXE, den_inp_file, den_bond_file, init_conf, seed, n_neighbours)
    with _silence():
        write_conf_from_last_frame(den_traj, den_conf)

    # 3. Main simulation from the denatured configuration
    _status("simulate")
    write_input_file(inp_file, fh, n0, 1, max_t, nsweep, dens)
    _run_quiet(_EXE, inp_file, bond_file, den_conf, seed, n_neighbours)

    # 4. Sample yield at each requested time-point.
    #    frames[0] = denatured start, frames[t] = after t steps.
    _, _, _, frames = parse_traj(traj_file)
    closing_pairs   = find_native_red_closing_pairs(n0)

    yields = {}
    for t in cfg['nsteps_curves']:
        idx       = min(t, len(frames) - 1)
        yields[t] = compute_yield(frames[idx], closing_pairs, L)

    return yields


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(e_reds, nsteps_list, results, cfg, out_dir):
    """Plot mean yield ± SEM vs e_red, one curve per time-point."""
    fig, ax = plt.subplots(figsize=(8, 5))

    colours = cm.viridis(np.linspace(0.1, 0.9, len(nsteps_list)))

    for colour, t in zip(colours, nsteps_list):
        means, sems = [], []
        for e in e_reds:
            vals = [v for v in results[t][e] if not math.isnan(v)]
            if vals:
                means.append(np.mean(vals))
                sems.append(np.std(vals) / math.sqrt(len(vals)) if len(vals) > 1 else 0.0)
            else:
                means.append(float('nan'))
                sems.append(0.0)

        ax.errorbar(
            e_reds, means, yerr=sems,
            marker='o', markersize=5, linewidth=1.8, capsize=3,
            label=f"t = {t} steps",
            color=colour,
        )

    ax.set_xlabel("RED bond energy  $e_\\mathrm{red}$  (kBT)", fontsize=12)
    ax.set_ylabel("RED bond yield\n(fraction of 2×2 closing bonds formed)", fontsize=11)

    subtitle = (
        f"n={cfg['polymer']},  L={cfg['L']},  green={cfg['hier_green']},  "
        f"blue={cfg['hier_blue']},  spring_k={cfg['spring_k']},  "
        f"unspecific={cfg['unspecific']},  {cfg['n_repeats']} repeat(s)"
    )
    ax.set_title(
        "Kinetic vs thermodynamic yield  —  RED level\n" + subtitle,
        fontsize=10,
    )

    ax.set_xlim(left=0)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(title="Simulation time", fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    plot_path = os.path.join(out_dir, 'red_yield_scan.png')
    fig.savefig(plot_path, dpi=150)
    print(f"Plot → {plot_path}")
    plt.show()


def plot_final_states(cfg, e_reds, out_dir, closing_pairs, results, nsteps_list):
    """Grid of lattice snapshots — one panel per e_red value.

    Each snapshot shows the final frame (longest requested time, rep=0) of the
    polymer on the lattice.  Particles that belong to a fully assembled 2×2
    sub-group (all 3 native RED contacts formed) are coloured green; all others
    are white.  Backbone bonds are drawn as thin grey lines.
    """
    import matplotlib.patches as mpatches

    n0   = cfg['polymer']
    L    = cfg['L']
    t_final = nsteps_list[-1]

    bb_bonds    = snake_path(n0)
    chain_order = _chain_order(n0)

    # Grid layout
    N     = len(e_reds)
    ncols = math.ceil(math.sqrt(N))
    nrows = math.ceil(N / ncols)

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.2 * ncols, 3.4 * nrows),
                             squeeze=False)
    axes_flat = axes.flatten()

    for panel_idx, e_red in enumerate(e_reds):
        ax = axes_flat[panel_idx]
        e_str = f"{e_red:.4g}".replace('.', 'p').replace('-', 'm')

        # Load final frame from the first available repeat
        frame = None
        for rep in range(cfg['n_repeats']):
            traj = os.path.join(out_dir, f"e{e_str}_r{rep}", 'sim_traj.txt')
            if os.path.exists(traj):
                try:
                    _, _, _, frames = parse_traj(traj)
                    frame = frames[min(t_final, len(frames) - 1)]
                    break
                except Exception:
                    continue

        if frame is None:
            ax.set_title(f"$e_{{\\rm red}}$={e_red:.3g}\n(no data)", fontsize=8)
            ax.axis('off')
            continue

        # Which particles are fully assembled?
        green_set = assembled_particles(frame, closing_pairs, L, n0)

        ax.set_facecolor('#f5f5f5')
        ax.set_xlim(0, L)
        ax.set_ylim(0, L)
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])

        # Backbone bonds (thin grey, PBC-aware minimum-image)
        for (pi, pj) in bb_bonds:
            xi, yi = frame[pi, 0] + 0.5, frame[pi, 1] + 0.5
            xj, yj = frame[pj, 0] + 0.5, frame[pj, 1] + 0.5
            dx, dy = xj - xi, yj - yi
            dx -= L * round(dx / L)
            dy -= L * round(dy / L)
            ax.plot([xi, xi + dx], [yi, yi + dy],
                    color='#888888', linewidth=0.7, alpha=0.7, zorder=1)

        # Formed native RED contacts (green dashed lines, PBC-aware)
        formed = _pairs_formed(frame, closing_pairs, L)
        for (pi, pj) in formed:
            xi, yi = frame[pi, 0] + 0.5, frame[pi, 1] + 0.5
            xj, yj = frame[pj, 0] + 0.5, frame[pj, 1] + 0.5
            dx, dy = xj - xi, yj - yi
            dx -= L * round(dx / L)
            dy -= L * round(dy / L)
            ax.plot([xi, xi + dx], [yi, yi + dy],
                    color='#27ae60', linewidth=1.0, linestyle='--',
                    alpha=0.85, zorder=1)

        # Particles
        for pid in range(n0):
            x, y   = frame[pid, 0] + 0.5, frame[pid, 1] + 0.5
            colour = '#27ae60' if pid in green_set else 'white'
            circ = mpatches.Circle((x, y), radius=0.35,
                                   facecolor=colour,
                                   edgecolor='#333333', linewidth=0.5,
                                   zorder=2)
            ax.add_patch(circ)

        # Mean yield across repeats for this energy
        vals = [v for v in results.get(t_final, {}).get(e_red, [])
                if not math.isnan(v)]
        mean_y = np.mean(vals) if vals else float('nan')
        ax.set_title(f"$e_{{\\rm red}}$ = {e_red:.3g}   yield = {mean_y:.2f}",
                     fontsize=8, pad=3)

    # Hide unused panels
    for ax in axes_flat[N:]:
        ax.axis('off')

    fig.suptitle(
        f"Final states (t = {t_final} steps)  —  green fill = fully assembled 2×2 sub-group  |  green dashed = formed RED contact",
        fontsize=10, y=1.01,
    )
    plt.tight_layout()

    path = os.path.join(out_dir, 'final_states.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f"Plot → {path}")
    plt.show()


# ---------------------------------------------------------------------------
# CSV save / load
# ---------------------------------------------------------------------------

def save_csv(path, e_reds, nsteps_list, results):
    header = 'e_red,repeat,' + ','.join(f'yield_t{t}' for t in nsteps_list)
    n_repeats = max(len(results[nsteps_list[0]][e]) for e in e_reds)
    with open(path, 'w') as f:
        f.write(header + '\n')
        for e in e_reds:
            for rep in range(n_repeats):
                vals = ','.join(
                    f"{results[t][e][rep]:.6f}" if rep < len(results[t][e])
                    else 'nan'
                    for t in nsteps_list
                )
                f.write(f"{e},{rep},{vals}\n")


def load_csv(path, e_reds, nsteps_list):
    """Reconstruct results dict from a saved CSV."""
    results = {t: {e: [] for e in e_reds} for t in nsteps_list}
    with open(path) as f:
        header = f.readline().strip().split(',')
        t_cols = {int(h.replace('yield_t', '')): i
                  for i, h in enumerate(header) if h.startswith('yield_t')}
        for line in f:
            parts = line.strip().split(',')
            if not parts or not parts[0]:
                continue
            e = float(parts[0])
            if e not in results[nsteps_list[0]]:
                continue
            for t, col in t_cols.items():
                if t in results:
                    val = float(parts[col])
                    results[t][e].append(val)
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Scan RED bond energy and plot kinetic/thermodynamic yield curves."
    )
    parser.add_argument('config', help="Path to scan config file (INI format)")
    parser.add_argument('--plot-only', action='store_true',
                        help="Skip simulations; replot from existing yields.csv")
    parser.add_argument('--plot-lattice', action='store_true', default=False,
                        help="Also show the final-state lattice grid (one panel per e_red). "
                             "Default: off.")
    args = parser.parse_args()

    cfg = load_config(args.config)

    # Resolve output directory (relative paths are relative to the script)
    out_dir = cfg['output_dir']
    if not os.path.isabs(out_dir):
        out_dir = os.path.join(_SCRIPT_DIR, out_dir)
    os.makedirs(out_dir, exist_ok=True)

    csv_path = os.path.join(out_dir, 'yields.csv')

    e_reds      = cfg['hier_red_values']
    nsteps_list = cfg['nsteps_curves']
    n_repeats   = cfg['n_repeats']

    # Brief header
    closing_pairs = find_native_red_closing_pairs(cfg['polymer'])
    n_groups = cfg['polymer'] // 4
    print(f"Scan: n0={cfg['polymer']}  "
          f"{len(e_reds)} energies × {n_repeats} repeats = {len(e_reds)*n_repeats} runs  |  "
          f"yield metric: {len(closing_pairs)} native RED contacts ({n_groups}×3)")

    # Check executable
    if not args.plot_only and not os.path.isfile(_EXE):
        sys.exit(
            f"ERROR: run_polymer executable not found at {_EXE}\n"
            f"Build it with:  cd {_SCRIPT_DIR} && make run_polymer"
        )

    # ---- Run or load results ----
    results = {t: {e: [] for e in e_reds} for t in nsteps_list}

    if args.plot_only:
        if not os.path.exists(csv_path):
            sys.exit(f"ERROR: --plot-only requested but no CSV found at {csv_path}")
        print(f"Loading results from {csv_path} ...")
        results = load_csv(csv_path, e_reds, nsteps_list)

    else:
        total = len(e_reds) * n_repeats
        done  = 0

        with open(csv_path, 'w') as f:
            f.write('e_red,repeat,' + ','.join(f'yield_t{t}' for t in nsteps_list) + '\n')

        for e_red in e_reds:
            e_str = f"{e_red:.4g}".replace('.', 'p').replace('-', 'm')
            for rep in range(n_repeats):
                work_dir = os.path.join(out_dir, f"e{e_str}_r{rep}")

                # Status callback updates the progress bar in place
                def _status(phase, _e=e_red, _r=rep, _d=done):
                    _progress(_d + 1, total, phase, _e, _r, n_repeats)

                _progress(done, total, "starting", e_red, rep, n_repeats)

                try:
                    yields = run_condition(cfg, e_red, work_dir, repeat=rep,
                                          status_fn=_status)
                except Exception as ex:
                    yields = {t: float('nan') for t in nsteps_list}
                    _progress(done + 1, total, f"FAILED: {ex}"[:20], e_red, rep, n_repeats)

                done += 1
                _progress(done, total, "done", e_red, rep, n_repeats)

                for t in nsteps_list:
                    results[t][e_red].append(yields.get(t, float('nan')))

                with open(csv_path, 'a') as f:
                    row_vals = ','.join(f"{yields.get(t, float('nan')):.6f}" for t in nsteps_list)
                    f.write(f"{e_red},{rep},{row_vals}\n")

        print(f"\nDone.  {csv_path}")

    # ---- Plots ----
    plot_results(e_reds, nsteps_list, results, cfg, out_dir)
    if args.plot_lattice:
        plot_final_states(cfg, e_reds, out_dir, closing_pairs, results, nsteps_list)


if __name__ == "__main__":
    main()
