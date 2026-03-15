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

# 64-bead polymer with 3-level hierarchical bonds (red > green > blue)
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --e1 1000 --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 --sim-seed 42

# Same but with denaturation pre-equilibration (10000 steps, weak bonds off)
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --e1 1000 --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 \
  --denature 10000 --sim-seed 42

# Hierarchical bonds with a soft spring backbone (k=50) instead of hard confinement
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 --spring-k 50 --sim-seed 42

# Non-specific mode: any two particles that share a colour level can bond (not just
# their specific Hilbert-curve neighbour).  Spring backbone keeps the chain topology.
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 --spring-k 50 --unspecific --sim-seed 42
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

**Hierarchical Hilbert couplings** (`--hier-red R --hier-green G --hier-blue B`): three-level hierarchy inspired by the Holmes-Cerfon & Wyart paper, applied to the Hilbert-curve chain. Couplings are assigned only between Hilbert-grid-adjacent non-backbone pairs, with strength determined by where the two particles sit in the chain traversal order:

| Level | Chain-position condition | Coupling | Animation colour |
|-------|--------------------------|----------|-----------------|
| RED   | same group of 4 consecutive chain positions | `--hier-red` | red |
| GREEN | same top-level quarter (N/4 positions), different group of 4 | `--hier-green` | green |
| BLUE  | different top-level quarters | `--hier-blue` | blue |

For N=64 (order-3 Hilbert curve) this gives 16 groups of 4 (RED), 4 groups of 16 (GREEN within a quarter), and inter-quarter contacts (BLUE). Setting `--hier-blue 0` allows the four quadrant blocks to separate from each other while each block remains internally rigid.

Each coupling fires only at its native Hilbert distance:
- Cardinal Hilbert contacts (grid d=1) are placed in `wD1` — stretching to physical d=√2 costs +val, causing VMMC to recruit the whole bonded group as a rigid cluster.
- Diagonal Hilbert contacts (grid d=√2) are placed in `wDsq2` — stretching to physical d=2 likewise triggers cluster recruitment.

Backbone-bonded pairs are excluded from all weak coupling so the chain remains freely articulated between d=1 and d=√2 at its hinge points.

Requires √N to be a power of 2.

**Spring backbone** (`--spring-k K`): replaces hard-wall confinement with a pure Hookean spring E(d) = k·(d−1)²:

| Distance | Spring energy |
|----------|--------------|
| d = 1    | 0 (equilibrium; level coupling also fires here) |
| d = √2   | k·(√2−1)² ≈ 0.172k |
| d = 2    | k |
| d = √5   | k·(√5−1)² ≈ 1.528k  ← maximum within cell-list range |
| d > √5   | 0 (C++ cell-list cutoff, but barrier ≈ 1.528k makes escape negligible for k >> 1) |

The spring is implemented **entirely in C++** (`StickySquare::computePairEnergy`), not via the coupling matrices. This is necessary because the coupling matrices are only sampled at the five discrete lattice distances; the spring must also act at arbitrarily large separations (e.g. during a VMMC cluster move that temporarily stretches a backbone bond). The C++ override works as follows:

1. `computeInteractions` is overridden to always add backbone partners to the neighbour list, regardless of their current physical distance (bypassing the cell-list cutoff).
2. `computePairEnergy` is overridden to return `k·(d−1)²` for any backbone pair, computed from the actual MIC distance, plus the `wD1` level coupling at d=1.
3. `computeEnergy` is overridden accordingly to use the above two methods.

Backbone pairs are **included** in the hierarchical level coupling at d=1. `--e1` is not used in spring mode.

**Coupling matrix display with spring backbone:** even though the spring energy is computed on the fly in C++, the coupling matrix panels in the animation display are still populated with the correct **static** spring penalties for display purposes:
- `D1`: level coupling only (spring is zero at equilibrium d=1)
- `Dsq2`: repulsive spring penalty k·(√2−1)² for backbone pairs; level coupling for non-backbone Hilbert-diagonal pairs
- `D2`, `Dsq5`: repulsive spring penalties for backbone pairs
- `D0`: uniform repulsive marker for **all** particle pairs (representing the universal hard-sphere exclusion at d=0; C++ handles this as INF regardless)

In the colour scheme: **red = repulsive**, **blue = attractive**.

**Important:** do not substitute a flat hard-wall at d=2 and d=√5. Equal energy at both distances creates a plateau, after which d>√5 (energy 0) is thermodynamically downhill — escape is actually favoured. The pure spring avoids this by making d=√5 a strict maximum.

**Non-specific bonding** (`--unspecific`): in the default (specific) mode, level coupling fires only when particle i is physically adjacent to its specific Hilbert-curve partner j at the correct distance. With `--unspecific`, **any two particles** whose chain positions place them in the same colour level can attract each other whenever they happen to be at physical d=1 or d=√2:

| Level | Condition | Energy at d=1 and d=√2 |
|-------|-----------|------------------------|
| RED   | same group of 4 chain positions (`cp // 4` identical) | `--hier-red` |
| GREEN | same top-level quarter, different group of 4 | `--hier-green` |
| BLUE  | different top-level quarters | `--hier-blue` |

This is the maximum generalisation: every close pair of particles gets some attractive coupling based on their chain-position relationship, not just those that are neighbours in the Hilbert curve.

*Will the Hilbert curve still be the lowest-energy structure?* **Probably not.** In specific mode the Hilbert curve uniquely maximises coupling because each coupling fires only for its one specific partner. In non-specific mode the energy depends on how many same-level contacts happen to be at d=1 or d=√2. Compact cluster configurations — where all members of a RED group (4 particles) are packed into a 2×2 block, providing 6 pairwise contacts at RED energy rather than the 3 available in a linear chain segment — can have lower energy than the Hilbert-curve layout. The competition is controlled by the ratio `e_red / k`: for large `k` the backbone stiffness prevents clustering; for large `e_red` compact clusters win. The Hilbert curve is the ground state only in the limit k → ∞.

In the animation, `--unspecific` causes bond lines to be drawn between **any** physically close pair (not just Hilbert-adjacent pairs), coloured by level: red/green/blue as appropriate.

**Denaturation pre-equilibration** (`--denature STEPS`): before the main simulation, runs `STEPS` steps with all weak coupling matrices zeroed (backbone confinement only). The final configuration is used as the starting point for the main run. This disperses the initial Hilbert-curve assembly so the system explores conformational space from a randomised starting state rather than the ground-state configuration.

### Animation window (polymer mode)

The window contains three regions:

| Panel | Contents |
|-------|----------|
| Left (full height) | Animated lattice: coloured particle discs (black outline), two-layer bond drawing (see below), white background |
| Top-right | Total energy and running average vs simulation step |
| Bottom-right (×5) | Coupling matrices D0, D1, D√2, D2, D√5 as static physical-energy heatmaps |

**Two-layer bond drawing:** bonds are drawn in two overlapping layers so that backbone topology and hierarchical coupling can both be seen at once.
1. **Coloured layer** (wide, semi-transparent, zorder=1): red/green/blue according to the hierarchical level of the pair. Backbone pairs also receive their level colour here, so the chain's colour hierarchy is always visible.
2. **Black backbone layer** (thin, opaque, zorder=2): drawn on top for every backbone pair regardless of physical distance (i.e. stretched backbone bonds are shown crossing the periodic boundary if necessary).

**Coupling matrix colour scale:** symmetric around zero; **red = repulsive (positive physical energy)**, **blue = attractive (negative physical energy)**. The matrix sign convention is: positive stored value = attractive (C++ returns `−stored` as the pair energy). The five panels show:
- `D0`: hard-sphere repulsion marker (uniform for all pairs; C++ returns +∞ regardless)
- `D1`: level coupling values; for spring mode, backbone pairs show their level coupling (zero if `--hier-*` are all 0)
- `Dsq2`, `D2`, `Dsq5`: level coupling for non-backbone Hilbert pairs; in spring mode, backbone pairs show static spring penalties k·(d−1)² (repulsive, red)

The coupling matrices are **static** — they are computed once from the bond file parameters and do not change during the animation.

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
| `--hier-red R` | RED bond strength (finest level: same group of 4 chain positions). Activates hierarchical Hilbert mode. |
| `--hier-green G` | GREEN bond strength (middle level: same top-level quarter, different group of 4) |
| `--hier-blue B` | BLUE bond strength (coarsest level: different top-level quarters) |
| `--spring-k K` | Backbone spring stiffness: pure Hookean spring E(d)=k·(d−1)² implemented in C++ at all distances. Energy increases strictly with d; barrier ≈1.528k at d=√5. Use k >> 1 (e.g. k=10–50) to prevent thermal escape. Backbone pairs also receive level coupling at d=1. `--e1` is not used. |
| `--unspecific` | Non-specific bonding: every particle pair (i,j) can form a bond at the level energy corresponding to their chain positions (RED/GREEN/BLUE), regardless of whether they are Hilbert-curve neighbours. The coupling fires at any physical d=1 or d=√2. Only meaningful with `--hier-*`. |
| `--denature N` | Run N pre-equilibration steps with weak bonds zeroed before the main simulation |
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
# (empty when using confinement encoding; directional backbone bonds listed here otherwise)
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

# Optional spring backbone sections (written when --spring-k is set):

SPRING_K
10.000000
SPRING_K_END

SPRING_BACKBONE
0 1
1 2
2 3
...
SPRING_BACKBONE_END
```

Entries in WEAK_D* sections are upper-triangle (i ≤ j); the matrix is symmetric. Positive values are attractive; negative values are repulsive.

The `SPRING_K` section contains a single float: the spring stiffness k. The `SPRING_BACKBONE` section lists all backbone bond pairs (one per line, both orderings accepted — the C++ reader symmetrises). When these sections are present, the C++ simulation ignores the wDsq2/wD2/wDsq5 entries for backbone pairs and computes E(d) = k·(d−1)² on the fly instead. **Every backbone pair must be listed** — pairs listed as `pi pj` only (not `pj pi`) are fine because the reader adds both directions.

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
- `<filehead>_denature_bonds.txt`: denatured bond file (weak bonds zeroed, `--denature` mode)
- `<filehead>_denature_traj.txt`: trajectory from denaturation phase
- `<filehead>_denature_final.conf`: final configuration of denaturation phase (used as main-run seed)

---

## Dependencies

- C++11 compiler (`g++`)
- Python 3 with `numpy`, `matplotlib`

```bash
pip install numpy matplotlib
```
