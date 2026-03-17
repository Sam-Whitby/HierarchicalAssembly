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

# 4-bead polymer, hard-confinement backbone, random weak couplings
python run_and_plot.py --polymer 4 --L 12 --nsteps 100000 --nsweep 1 \
  --e1 1000 --weak-e 1.0 --weak-std 0.5 --weak-seed 42 --sim-seed 7

# Boltzmann validation (pre-enumerated, hard-confinement backbone)
python run_and_plot.py --polymer 4 --L 12 --nsteps 1000000 --nsweep 1 \
  --e1 1000 --weak-e 1.0 --weak-std 0.5 --weak-seed 42 --sim-seed 7 --boltzmann

# 4-bead polymer with spring backbone k=1, zero level coupling, live Boltzmann validation
python run_and_plot.py --polymer 4 --L 12 --nsteps 100000 --nsweep 1 \
  --hier-red 0 --spring-k 1 --sim-seed 42 --boltzmann --enumerate-live

# 64-bead polymer with 3-level hierarchical bonds
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 --sim-seed 42

# Same with spring backbone (k=50) and denaturation pre-equilibration
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 --spring-k 50 \
  --denature 10000 --sim-seed 42

# Non-specific bonding: any same-level pair can bond, not just Hilbert neighbours
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 --spring-k 50 --unspecific --sim-seed 42

# Enable rotational moves (50% rotations, 50% translations)
python run_and_plot.py --polymer 4 --L 12 --nsteps 100000 --nsweep 1 \
  --hier-red 0 --spring-k 1 --sim-seed 42 --prob-translate 0.5

# Boltzmann validation with rotational moves (requires --enumerate-live for spring backbone)
python run_and_plot.py --polymer 4 --L 12 --nsteps 1000000 --nsweep 1 \
  --hier-red 0 --spring-k 1 --sim-seed 42 --prob-translate 0.5 --boltzmann --enumerate-live

# Scan RED bond energy from 0..20 with multiple repeats and time curves (edit scan_config.ini first)
python scan_assembly.py scan_config.ini
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

`--polymer N` automates polymer simulation of an N-bead chain.

### Chain topology

Particle identities are arranged in a boustrophedon (snake) path — or Hilbert space-filling curve for power-of-2 √N — so that chain-adjacent particles have positions close in the Hilbert-curve order. The chain order determines which particles are backbone-bonded and which hierarchical coupling level applies to each pair.

### Backbone confinement

Two options, selected via the absence or presence of `--spring-k`:

**Hard-wall confinement** (default, `--e1 E`): strongly repulsive entries are added to `wD2` and `wDsq5` for each backbone-bonded pair, trapping every bond at d ∈ {1, √2}. The BACKBONE section of the bond file is empty; all physics is in the coupling matrices.

**Spring backbone** (`--spring-k K`): replaces hard-wall with a Hookean spring E(d) = k·(d−1)². Implemented **entirely in C++** (`StickySquare::computePairEnergy`) at any bond distance:

1. `computeInteractions` is overridden to always add backbone partners to the neighbour list, bypassing the cell-list cutoff.
2. `computePairEnergy` returns `k·(d−1)²` for backbone pairs (from actual MIC distance), plus `wD1` level coupling at d=1.
3. `computeEnergy` is overridden to use the above two methods.

Choose k small enough that conformational changes are thermally accessible: with VMMC lattice moves, the backbone partner is recruited with probability 1−exp(−k·(√2−1)²) when a bond would stretch, so for large k the chain moves as a rigid body and conformational changes are rare.

### Coupling modes

**Random weak couplings** (default): each matrix entry drawn from N(weak_e, weak_std). Use `--weak-seed` for reproducibility.

**Hilbert-local couplings** (`--hilbert-A A --hilbert-B B`): non-zero couplings only between Hilbert-grid-adjacent non-backbone pairs (A for cardinal neighbours, B for diagonal). Requires √N to be a power of 2.

**Hierarchical Hilbert couplings** (`--hier-red R --hier-green G --hier-blue B`): three-level hierarchy applied to Hilbert-curve chain. Coupling strength is set by where the two particles sit in the chain traversal order:

| Level | Condition | Coupling | Animation colour |
|-------|-----------|----------|-----------------|
| RED   | same group of 4 consecutive chain positions | `--hier-red` | red |
| GREEN | same top-level quarter, different group of 4 | `--hier-green` | green |
| BLUE  | different top-level quarters | `--hier-blue` | blue |

In specific mode (default): each coupling fires only between Hilbert-grid-adjacent pairs at their native Hilbert distance (d=1 → `wD1`; d=√2 → `wDsq2`). Backbone-bonded pairs receive level coupling at d=1 (spring mode) or no coupling (hard-confinement).

In non-specific mode (`--unspecific`): every pair (i,j) can attract at **any** physical d=1 or d=√2, based on their chain-position relationship. The Hilbert curve is no longer guaranteed to be the ground state — compact same-level clusters can have more pairwise contacts. The competition is controlled by the ratio `e_red / k`: for large `k` backbone stiffness prevents clustering; for large `e_red` compact clusters win.

**Denaturation pre-equilibration** (`--denature STEPS`): runs STEPS steps with all weak couplings zeroed (backbone only) before the main simulation, dispersing the initial assembled configuration.

### Animation window

| Panel | Contents |
|-------|----------|
| Left (full height) | Animated lattice: coloured particle discs, two-layer bond drawing |
| Top-right | Total energy and running average vs simulation step |
| Bottom-right (×5) | Coupling matrices D0, D1, D√2, D2, D√5 as static heatmaps |

**Two-layer bond drawing:**
1. **Coloured layer** (wide, semi-transparent): red/green/blue by hierarchical level; backbone pairs also receive their level colour.
2. **Black backbone layer** (thin, on top): drawn for every backbone pair at any physical distance, including stretched bonds across periodic boundaries.
3. **Green dashed lines**: when yield > 0 but no fully assembled groups exist yet, dashed lines connect particle pairs that have formed native RED contacts (i.e. the bonds that contribute to the yield calculation). Useful for diagnosing partial assembly.

**Coupling matrix panels** show static physical energies: red = repulsive, blue = attractive. For spring backbone, backbone pairs show spring penalties k·(d−1)² in the Dsq2/D2/Dsq5 panels; D0 shows a uniform repulsive marker for all pairs (representing universal hard-sphere exclusion).

### Energy sign convention

**Negative physical energy = attractive / favourable.**

The C++ engine stores `pair_energy = −coupling_value` (positive coupling → attractive). Boltzmann weight:

```
P(state) ∝ g × exp(+E_python) = g × exp(−E_physical)
```

where `E_python = sum of coupling values` and `E_physical = −E_python`.

### Boltzmann validation (`--boltzmann`)

Checks that the simulation samples the correct Boltzmann distribution. Two modes:

**Pre-enumerated** (default, `--boltzmann` without `--enumerate-live`):
- Enumerates all canonically distinct conformations by depth-first self-avoiding walk on the 8-connected lattice.
- Energies computed from coupling matrices.
- Shows: canonical-state diagram (all enumerated states sorted by energy) + Boltzmann correlation figure.
- **Suitable for hard-confinement backbone**, where backbone bonds are restricted to d=1 or d=√2. With spring backbone, longer-bond states are missing from the enumeration and the comparison is incomplete.

**Live-enumerated** (`--boltzmann --enumerate-live`):
- Builds the state set entirely from simulation observations; no prior enumeration.
- Backbone spring energy is computed exactly at the observed distance (not from matrix values).
- Shows: live state diagram (up to 50 states; total count shown if more) + Boltzmann correlation figure.
- **Required for spring backbone** (`--spring-k`), since bonds can stretch to d=2, √5, etc. Also handles chains that temporarily span more than L/2 across periodic boundaries via cumulative minimum-image bond-vector unrolling.

Both figures show:
- Left: Pearson r between observed and Boltzmann-predicted frequencies vs simulation frame (log scale). Convergence to r ≈ 1 confirms detailed balance.
- Right: Scatter of observed vs predicted frequency at the final frame, coloured by physical energy.

**Kinetic convergence:** for spring backbone with large k, VMMC recruits the backbone partner with probability ≈ 1−exp(−k·(√2−1)²) whenever a bond would stretch, making the chain rigid and conformational changes rare. Use k ≲ 3–4 for 100k-step validation runs with n=4.

---

## Parameter scan (`scan_assembly.py`)

Sweeps RED bond energy over a range of values and time points, running multiple independent repeats per condition, and produces a yield-vs-energy plot showing the kinetic/thermodynamic crossover.

```bash
python scan_assembly.py scan_config.ini
```

Configure `scan_config.ini`:

```ini
[simulation]
polymer        = 64          # chain length (sqrt must be power of 2)
L              = 40          # box side length
nsweep         = 1           # sweeps per output step
hier_green     = 0.0         # GREEN and BLUE couplings (set to 0 to isolate RED level)
hier_blue      = 0.0
spring_k       = 0.05        # spring backbone stiffness (or None for hard-wall)
e1             = 8.0         # hard-wall backbone energy (used during denaturation phase)
sim_seed       = 42          # base RNG seed (repeat r uses seed + r)
unspecific     = true        # non-specific bonding
n_neighbours   = 8           # VMMC lattice directions: 4 or 8
denature_steps = 1000        # pre-equilibration steps with weak bonds zeroed

[scan]
hier_red_values = 0,1,2,3,4,5,6,7,8,9,10,20   # RED energies to test
nsteps_curves   = 1,10,100,1000                # time points for yield measurement

[output]
output_dir = scan_results    # directory for per-condition subdirs + summary CSV + plot
n_repeats  = 5               # independent repeats per condition
```

**Output:** `scan_results/yields.csv` (per-condition yields) and a yield-vs-energy PNG with one curve per time point and error bars across repeats.

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
| `--e1 E` | Hard-wall backbone confinement energy (not used with `--spring-k`). |
| `--weak-e E` | Mean weak coupling (default: 1.0) |
| `--weak-std S` | Std dev of weak couplings (default: 0.5) |
| `--weak-seed N` | RNG seed for coupling matrix generation |
| `--hilbert-A A` | Coupling for Hilbert cardinal neighbours. Activates Hilbert-local mode. |
| `--hilbert-B B` | Coupling for Hilbert diagonal neighbours |
| `--hier-red R` | RED coupling (finest level). Activates hierarchical Hilbert mode. |
| `--hier-green G` | GREEN coupling (middle level) |
| `--hier-blue B` | BLUE coupling (coarsest level) |
| `--spring-k K` | Spring backbone stiffness k: E(d)=k·(d−1)², computed in C++ at all distances. Replaces `--e1` hard-wall confinement. |
| `--unspecific` | Non-specific bonding: all same-level pairs bond at any d=1 or d=√2, not just Hilbert neighbours. Only meaningful with `--hier-*`. |
| `--denature N` | Pre-equilibration: N steps with weak bonds zeroed before main run |
| `--sim-seed N` | RNG seed for C++ simulation |
| `--prob-translate P` | Fraction of VMMC moves that are translations (default 1.0). Set < 1 to enable rotational moves (e.g. 0.5 for 50/50). |
| `--n-neighbours N` | Lattice directions for VMMC translation moves: 4 = cardinal only, 8 = cardinal + diagonal (default 4). |
| `--boltzmann` | Enable Boltzmann validation (see below) |
| `--enumerate-live` | Use live-observed state enumeration for Boltzmann validation (required for `--spring-k`). Without this flag, uses pre-enumerated DFS states. |

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
# (empty in confinement-encoding mode)
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

# Spring backbone sections (written when --spring-k is set):

SPRING_K
10.000000
SPRING_K_END

SPRING_BACKBONE
0 1
1 3
3 2
...
SPRING_BACKBONE_END
```

Entries in WEAK_D* sections are upper-triangle (i ≤ j); the matrix is symmetric. Positive values are attractive; negative values are repulsive.

When `SPRING_K` and `SPRING_BACKBONE` sections are present, C++ computes E(d) = k·(d−1)² on the fly for backbone pairs and ignores their wDsq2/wD2/wDsq5 matrix entries. Every backbone bond need only be listed once (i j); the C++ reader symmetrises automatically.

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

Dynamics: **Virtual Move Monte Carlo** (cluster moves). Two move types:

- **Translation**: cluster moves by one lattice step in a cardinal (or optionally diagonal) direction.
- **Rotation**: cluster rotates rigidly about the seed particle by 90°, 180°, or 270°. Only multiples of 90° keep all particles on integer sites. A random neighbor is always force-included so that the rotation of a single isotropic particle (which has no effect) is never proposed. The VMMC Stokes acceptance factor for rotational drag scales as the cube of the cluster's RMS distance from the rotation center.

The fraction of moves that are translations vs. rotations is set by `--prob-translate` (default 1.0 = translations only).

---

## Repository contents

| File/Directory | Description |
|----------------|-------------|
| `makefile` | Build system |
| `run_hier.cpp` | Hierarchical assembly driver |
| `run_custom.cpp` | Custom bond driver |
| `run_polymer.cpp` | Polymer driver |
| `run_and_plot.py` | Python wrapper: runner, animation, Boltzmann validation |
| `scan_assembly.py` | Parameter scan over RED bond energies and time points, producing yield curves |
| `scan_config.ini` | Configuration file for `scan_assembly.py` |
| `src/` | Library source and headers |
| `input_hier.txt` | Default input file |
| `old/` | Deprecated files |

### Input file format

```
filehead     # output file prefix
n            # target structure size
ncopies      # number of copies
nsteps       # number of output steps
nsweep       # MC sweeps per step
density      # volume fraction
```

### Output files (gitignored)

- `<filehead>_traj.txt`: particle trajectories
- `<filehead>_stats.txt`: per-step energy and fragment-size histogram
- `<filehead>_polymer_bonds.txt`: generated bond file (polymer mode)
- `<filehead>_init.conf`: initial configuration (polymer mode)
- `<filehead>_denature_bonds.txt`, `_denature_traj.txt`, `_denature_final.conf`: denaturation phase outputs

---

## Dependencies

- C++11 compiler (`g++`)
- Python 3 with `numpy`, `matplotlib`

```bash
pip install numpy matplotlib
```
