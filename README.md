# HierarchicalAssembly

Simulates systems of sticky squares on a 2D square lattice using the Virtual Move Monte Carlo (VMMC) algorithm. Based on:

> Holmes-Cerfon, M. and Wyart, M., 2025. Hierarchical self-assembly for high-yield addressable complexity at fixed conditions. [arXiv:2501.02611](https://arxiv.org/abs/2501.02611)

The `StickySquares` class is adapted from [vmmc.xyz](http://vmmc.xyz/), extended for lattice particles with direction-specific patchy interactions.

---

## Quick start

```bash
mkdir -p obj          # once only

make                  # builds run_hier (hierarchical assembly)
make run_custom       # builds run_custom (arbitrary bond strengths)
make run_polymer      # builds run_polymer (polymer / confinement mode)

# Random 4-bead polymer — animation + coupling matrices
python run_and_plot.py --polymer 4 --L 12 --nsteps 100000 --nsweep 1 \
  --e1 1000 --weak-e 1.0 --weak-std 0.5 --weak-seed 42 --sim-seed 7

# Same run + Boltzmann validation and canonical-state diagram
python run_and_plot.py --polymer 4 --L 12 --nsteps 1000000 --nsweep 1 \
  --e1 1000 --weak-e 1.0 --weak-std 0.5 --weak-seed 42 --sim-seed 7 --boltzmann

# 64-bead polymer with Hilbert-local couplings (A=cardinal, B=diagonal)
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --e1 1000 --hilbert-A 1.0 --hilbert-B 0.5 --sim-seed 42
```

Press **spacebar** to pause/resume the animation.

---

## Drivers

| Binary | Source | Description |
|--------|--------|-------------|
| `run_hier` | `run_hier.cpp` | Hierarchical self-assembly with bond strengths e1, e1/2, e1/4, … at sub-block boundaries |
| `run_custom` | `run_custom.cpp` | Arbitrary omni-directional + diagonal bonds from a bond file |
| `run_polymer` | `run_polymer.cpp` | Polymer mode: extended bond file with BACKBONE + WEAK_D* sections, seeded RNG, pre-assembled initial config |

---

## Polymer mode (`--polymer N`)

`--polymer N` fully automates polymer simulation of an N-bead chain on a √N × √N grid.

### What it does

1. **Snake-path chain** — defines which consecutive particle-identity pairs are backbone-bonded.
2. **Confinement backbone** — backbone bonds encoded as repulsive coupling entries (see below). No directional attraction needed.
3. **Coupling matrices** — one symmetric N×N matrix for each of five distance classes: 0, 1, √2, 2, √5.
4. **Hilbert-curve initialisation** — places all N particles in the assembled configuration (no overlaps, all backbone constraints satisfied) using a Hilbert space-filling curve for power-of-2 √N, or a boustrophedon path otherwise.

### Confinement mechanism

Backbone chain bonds are **not** implemented as directional attractions. Instead, for each bonded pair (i, j), strongly repulsive entries are added to `wD2` and `wDsq5`:

```
coupling[i,j] -= e_backbone     (coupling << 0 → physical energy >> 0 → repulsive)
```

Combined with hard-sphere exclusion (d < 1 → ∞ energy), this traps every bonded pair at d ∈ {1, √2}:

| Distance | Energy |
|----------|--------|
| d = 0    | +∞ (hard sphere) |
| d = 1    | weak coupling only |
| d = √2   | weak coupling only |
| d = 2    | +e_backbone (strongly repulsive for bonded pair) |
| d = √5   | +e_backbone (strongly repulsive for bonded pair) |

All bonding physics is encoded in the five coupling matrices. The BACKBONE section of the bond file is left empty.

### Coupling modes

**Random weak couplings** (default): each matrix entry is drawn from N(weak_e, weak_std). Entries can be negative (repulsive). Use `--weak-seed` for reproducibility.

**Hilbert-local couplings** (`--hilbert-A A --hilbert-B B`): non-zero couplings assigned only between particles whose Hilbert-grid positions are adjacent:
- Hilbert cardinal neighbours (grid d=1) → `wD1[i,j] = A`
- Hilbert diagonal neighbours (grid d=√2) → `wDsq2[i,j] = B`
- Backbone-bonded pairs are excluded from A/B (already confined; their short-range coupling is irrelevant)
- All other matrix entries are zero

Requires √N to be a power of 2 (e.g. N = 4, 16, 64, 256). This mode is designed for studying how spatially local coupling within the Hilbert-curve path shapes the Boltzmann distribution over conformations.

### Animation window (polymer mode)

The window contains three regions:

| Panel | Contents |
|-------|----------|
| Left (full height) | Animated lattice: coloured particle discs (black outline), black backbone bond lines, white background |
| Top-right | Total energy and running average vs simulation step |
| Bottom-right (×5) | Coupling matrices D0, D1, D√2, D2, D√5 as physical-energy heatmaps (blue = attractive, red = repulsive) |

The coupling matrices are static. Their colour scale is symmetric about zero so the sign of each entry is immediately visible. D0 is display-only (hard sphere handled by C++).

### Boltzmann validation (`--boltzmann`)

Two additional figures open after the simulation:

**Canonical-state diagram** — all distinct conformations of the N-bead chain, enumerated by depth-first self-avoiding walk on the 8-connected 2D lattice and canonicalised under the full dihedral group (4 rotations × 2 reflections). States sorted from lowest to highest physical energy. Each panel shows: state index, degeneracy g, physical energy E.

**Boltzmann correlation figure**:
- Left: Pearson r between observed and Boltzmann-predicted frequencies vs simulation frame (log scale). Convergence to r ≈ 1 confirms detailed balance.
- Right: Scatter of observed vs predicted frequency at the final frame, coloured by physical energy.

### Energy sign convention

Throughout the code, **negative physical energy = attractive / favourable**.

The C++ engine stores `pair_energy = −coupling_value`, so a positive coupling value is attractive. The Boltzmann weight of a state is:

```
P(state) ∝ g × exp(+E_python) = g × exp(−E_physical)
```

where `E_python = sum of coupling values` and `E_physical = −E_python`. All plots use `E_physical`.

---

## Python wrapper flags

### All modes

| Flag | Description |
|------|-------------|
| `--n0 N` | Target structure size (default: 16) |
| `--e1 E` | Bond / backbone confinement energy (default: 8.0) |
| `--L N` | Box side length |
| `--nsteps N` | Number of output steps |
| `--nsweep N` | MC sweeps per step |
| `--dens D` | Particle density (alternative to `--L`) |
| `--ncopies N` | Number of target copies |
| `--filehead STR` | Prefix for output files (default: `hier`) |
| `--no-run` | Skip simulation, re-plot existing output files |

### Polymer mode

| Flag | Description |
|------|-------------|
| `--polymer N` | Polymer size N (perfect square). Sets n0=N, 1 copy. |
| `--e1 E` | Backbone confinement energy. Use large values (e.g. 1000) to keep chain assembled. |
| `--weak-e E` | Mean weak coupling value (default: 1.0) |
| `--weak-std S` | Std dev of weak couplings (default: 0.5) |
| `--weak-seed N` | RNG seed for coupling matrix generation |
| `--hilbert-A A` | Coupling strength A for Hilbert cardinal neighbours. Activates Hilbert-local mode. |
| `--hilbert-B B` | Coupling strength B for Hilbert diagonal neighbours |
| `--sim-seed N` | RNG seed for the C++ VMMC simulation |
| `--boltzmann` | Enable canonical-state diagram and Boltzmann validation |

---

## Bond file formats

### Simple format (`run_custom` / `--bond-file`)

```
# comment
particle_i  particle_j  energy
```

Bonds are omni-directional and active at both d=1 and d=√2.

### Extended format (`run_polymer`)

```
BACKBONE
# (empty when using confinement encoding)
BACKBONE_END

WEAK_D1
i j energy
...
WEAK_D1_END

WEAK_DSQRT2
...
WEAK_DSQRT2_END

WEAK_D2
...
WEAK_D2_END

WEAK_DSQRT5
...
WEAK_DSQRT5_END
```

Entries in WEAK_D* sections are upper-triangle (i ≤ j); the matrix is symmetric. Positive values are attractive; negative values are repulsive.

---

## Model overview

Each particle is a unit square on a 2D lattice. Pair interactions are computed up to distance √5:

| Distance | Interaction |
|----------|-------------|
| d < 1    | Hard-core exclusion (+∞) |
| d = 1    | Directional backbone (if set) + wD1 weak coupling |
| d = √2   | Direction-agnostic backbone (if set) + wDsq2 weak coupling |
| d = 2    | wD2 weak coupling only |
| d = √5   | wDsq5 weak coupling only |

Dynamics: **Virtual Move Monte Carlo** (cluster moves), efficiently sampling collective rearrangements including polymer conformational changes.

---

## Repository contents

| File/Directory | Description |
|----------------|-------------|
| `makefile` | Build system |
| `run_hier.cpp` | Hierarchical assembly driver |
| `run_custom.cpp` | Custom bond driver |
| `run_polymer.cpp` | Polymer driver |
| `run_and_plot.py` | Python wrapper: simulation runner, animated visualisation, Boltzmann validation |
| `src/` | All library source and header files |
| `input_hier.txt` | Default input file (used by all drivers via Python wrapper) |
| `old/` | Deprecated files (legacy attractive-backbone bond files, MATLAB scripts, old drivers) |

### Input file format (`input_hier.txt`)

```
filehead     # output file prefix
n            # target structure size
ncopies      # number of copies
nsteps       # number of output steps
nsweep       # MC sweeps per step
density      # volume fraction
```

### Output files (gitignored, generated at runtime)

- `<filehead>_traj.txt`: particle trajectories
- `<filehead>_stats.txt`: per-step energy and fragment-size histogram
- `<filehead>_polymer_bonds.txt`: generated extended bond file (polymer mode)
- `<filehead>_init.conf`: generated initial configuration (polymer mode)

---

## Dependencies

- C++11 compiler (`g++`)
- Python 3 with `numpy`, `matplotlib`

```bash
pip install numpy matplotlib
```
