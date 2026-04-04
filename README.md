# Nucleolus Assembly / HierarchicalAssembly

Simulates self-assembly of patchy-particle systems on a 2D square lattice using
the Virtual Move Monte Carlo (VMMC) algorithm.

Two simulation modes are available:

| Mode | Binary | Description |
|---|---|---|
| **Hierarchical assembly** | `run_hier`, `run_polymer`, `run_custom` | General-purpose Gō-model polymers with hierarchical couplings, Boltzmann validation, and parameter scanning. See [Hierarchical assembly](#hierarchical-assembly) below. |
| **Nucleolus assembly** | `run_nucleolus` | Nucleolus ribogenesis model from Chapter 2 of Sam Whitby's PhD thesis. Gradient-modulated assembly in a condensate column with removal/reinjection. See [Nucleolus assembly](#nucleolus-assembly) below. |

Based on:
> Holmes-Cerfon, M. and Wyart, M., 2025. Hierarchical self-assembly for high-yield addressable complexity at fixed conditions. [arXiv:2501.02611](https://arxiv.org/abs/2501.02611)

The `StickySquares` class is adapted from [vmmc.xyz](http://vmmc.xyz/), extended for lattice particles with direction-specific patchy interactions.

---

## Build

```bash
mkdir -p obj          # once only

make run_nucleolus    # nucleolus simulation (Chapter 2)
make                  # builds run_hier (hierarchical assembly)
make run_custom       # builds run_custom
make run_polymer      # builds run_polymer
```

Requires a C++11 compiler (`g++`). Python 3 with `numpy` and `matplotlib` for visualisation.

---

# Nucleolus assembly

## Physics

The system contains `nCopies = 4` copies of a 16-particle target complex `T`.
Each complex is made of 4 polymers × 4 segments, arranged on a 4×4 grid following
the n=2 Moore (space-filling) curve and partitioned into quadrants:

```
Polymer 3: (0,2)(1,2)(1,3)(0,3)  |  Polymer 2: (2,2)(3,2)(3,3)(2,3)
───────────────────────────────────┼──────────────────────────────────
Polymer 0: (0,0)(1,0)(1,1)(0,1)  |  Polymer 1: (2,0)(3,0)(3,1)(2,1)
```

The column has a hard wall at `x = 0`, periodic boundaries in `y` with period `W`,
and extends to large `x`. Particles can drift past `x = L` before being removed.

### Interactions (J = 8, ε = 0.5)

All weak coupling is stored in 16×16 matrices indexed by local particle id (0–15).

| Condition | d = 1 | d = √2 | d = 2 |
|---|---|---|---|
| Same polymer type | −J (repulsive) | −J | −εJ |
| Cross-type, Gō neighbours in T at d = 1 | +J | — | — |
| Cross-type, Gō neighbours in T at d = √2 | — | +εJ | — |
| All other cross-type pairs | 0 | 0 | 0 |

Backbone bonds (consecutive segments within each polymer) are stored as strong
`Triples` (value 1000 ≈ ∞) and are **never** scaled by the gradient.
Hard-core overlap (`d < 1`) is always forbidden.

### Chemical gradient

When `--gradient` is set, all weak couplings between particles `i` and `j` are
scaled by `γ(xᵢ) · γ(xⱼ)` where:

```
γ(x) = min(x / L, 1)
```

Near `x = 0` all interactions are suppressed (denatured / disordered regime).
Full coupling is reached at `x ≥ L` (assembled / folded regime).
Backbone bonds are unaffected.

### Saturated-Link (SL) moves

A fraction `φ_SL` of VMMC steps are designated Saturated-Link moves. During an
SL move the algorithm tracks which of the 16 particle types
(`particle_id % 16`) are already in the moving cluster. When considering
whether to recruit a neighbour, if that neighbour's type is already represented
the link is skipped. This removes the kinetic frustration caused by repulsive
same-type interactions and is required for efficient sampling in the repulsive/
attractive coupling scheme.

### Removal and replacement

After every iteration a BFS is run over the interaction graph to find connected
components. A component is **removed and reinserted** if:

1. It contains exactly 16 particles (one complete complex worth).
2. All its particles have `x > L`.
3. It has no non-backbone interactions with anything outside itself.

Removed particles are reinserted as 4 denatured horizontal chains near `x = 1`
in the first free 4×4 block found by scanning upward from `x = 1`. If no space
is available the component waits in place.

## Running a simulation

```
./run_nucleolus [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Total outer iterations. Each iteration = `N_particles` VMMC move attempts + one removal/replacement check. |
| `--snapshots N` | 1000 | Number of trajectory frames to save. Frames are written at evenly spaced intervals so the file stays manageable regardless of step count. |
| `--length L` | 60 | Condensate column length in lattice units. |
| `--width W` | 10 | Column width (y periodic with period W). |
| `--gradient` | off | Enable the linear chemical gradient. |
| `--stokes` | off | Stokes hydrodynamic drag: diffusion ∝ 1/R per cluster. Without this flag all cluster sizes diffuse equally. |
| `--phi-sl φ` | 0.2 | Fraction of SL moves (0 = none, 1 = all). |
| `--phi-rot φ` | 0.2 | Fraction of rotation moves. |
| `--output PREFIX` | `nucleolus` | Output file prefix → `PREFIX_traj.txt` and `PREFIX_stats.txt`. |
| `--seed S` | 0 | RNG seed; 0 uses a time-based seed. |

### Example

```bash
./run_nucleolus \
    --steps 100000  --snapshots 1000 \
    --length 40     --width 10       \
    --gradient      --stokes         \
    --phi-sl 0.2    --phi-rot 0.2    \
    --output my_run
```

Runs 100 000 iterations, saves 1000 frames to `my_run_traj.txt`, with the
chemical gradient and Stokes drag enabled.

## Output files

### `PREFIX_traj.txt` — trajectory

Extended XYZ format. Each frame:

```
<N_particles>
step=S energy=E exited=X L=LL W=WW nCopies=C
<id> <poly_type> <x> <y> <copy>
...
```

The header carries the simulation step, total system energy, cumulative
exited-complex count, and box parameters, so no separate file is needed for
time-series plots.

### `PREFIX_stats.txt` — scalar time series

Tab-separated: `step  energy  nExited  acceptRatio`. One row per saved frame.

## Visualising

```bash
python3 visualize_nucleolus.py PREFIX_traj.txt \
        --gradient-length 40 --width 10
```

Opens a 3-panel animated figure:

- **Left**: particle positions coloured by polymer type (hue = type, shade =
  copy), backbone bonds drawn as dark lines, blue gradient overlay from x=0 to
  x=L, dashed line at the column boundary.
- **Top right**: system energy vs simulation step. Full time series in grey,
  current frame highlighted in red.
- **Bottom right**: cumulative exited complexes vs step, highlighted in green.

| Flag | Default | Description |
|---|---|---|
| `--output FILE` | display | Save to `.mp4` (requires ffmpeg) or `.gif` (requires Pillow). |
| `--fps N` | 5 | Playback frame rate. |
| `--gradient-length L` | from file | Column length for gradient overlay. |
| `--width W` | from file | Column width. |
| `--skip N` | 1 | Render every N-th frame to speed up display of large files. |
| `--title TEXT` | automatic | Figure title. |

### Save to file

```bash
# MP4 (requires ffmpeg)
python3 visualize_nucleolus.py my_run_traj.txt \
        --gradient-length 40 --width 10 \
        --output my_run.mp4 --fps 10

# GIF
python3 visualize_nucleolus.py my_run_traj.txt \
        --gradient-length 40 --output my_run.gif
```

---

# Hierarchical assembly

## Quick start

```bash
# 4-bead polymer, hard-confinement backbone, random weak couplings
python run_and_plot.py --polymer 4 --L 12 --nsteps 100000 --nsweep 1 \
  --e1 1000 --weak-e 1.0 --weak-std 0.5 --weak-seed 42 --sim-seed 7

# Boltzmann validation
python run_and_plot.py --polymer 4 --L 12 --nsteps 1000000 --nsweep 1 \
  --e1 1000 --weak-e 1.0 --weak-std 0.5 --weak-seed 42 --sim-seed 7 --boltzmann

# 4-bead polymer with spring backbone k=1
python run_and_plot.py --polymer 4 --L 12 --nsteps 100000 --nsweep 1 \
  --hier-red 0 --spring-k 1 --sim-seed 42 --boltzmann --enumerate-live

# 64-bead polymer with 3-level hierarchical bonds
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 --sim-seed 42

# Same with spring backbone (k=50) and denaturation pre-equilibration
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 --spring-k 50 \
  --denature 10000 --sim-seed 42

# Non-specific bonding
python run_and_plot.py --polymer 64 --L 20 --nsteps 200000 --nsweep 1 \
  --hier-red 3.0 --hier-green 2.0 --hier-blue 1.0 --spring-k 50 --unspecific --sim-seed 42

# Enable rotational moves (50/50)
python run_and_plot.py --polymer 4 --L 12 --nsteps 100000 --nsweep 1 \
  --hier-red 0 --spring-k 1 --sim-seed 42 --prob-translate 0.5

# Parameter scan
python scan_assembly.py scan_config.ini
```

Press **spacebar** to pause/resume the animation.

## Drivers

| Binary | Source | Description |
|---|---|---|
| `run_hier` | `run_hier.cpp` | Hierarchical self-assembly with bond strengths e1, e1/2, e1/4, … |
| `run_custom` | `run_custom.cpp` | Arbitrary omni-directional + diagonal bonds from a bond file |
| `run_polymer` | `run_polymer.cpp` | Polymer mode: BACKBONE + WEAK_D* sections, seeded RNG, pre-assembled initial config |

## Polymer mode (`--polymer N`)

`--polymer N` simulates an N-bead chain arranged on a Hilbert / boustrophedon
path on the lattice. Backbone-bonded pairs are those consecutive in the chain.

### Backbone

**Hard-wall confinement** (default): strongly repulsive wD2/wDsq5 entries trap
backbone bonds at d ∈ {1, √2}.

**Spring backbone** (`--spring-k K`): Hookean spring E(d) = k·(d−1)² computed
in C++ at any distance. Choose small k (≲ 3–4) for flexible chains with
efficient sampling.

### Coupling modes

| Mode | Flags | Description |
|---|---|---|
| Random | `--weak-e E --weak-std S` | Each matrix entry ~ N(E, S) |
| Hilbert-local | `--hilbert-A A --hilbert-B B` | Non-zero only between Hilbert-grid-adjacent pairs |
| Hierarchical | `--hier-red R --hier-green G --hier-blue B` | Three-level hierarchy: RED = same group of 4, GREEN = same quarter, BLUE = across quarters |

**Non-specific** (`--unspecific`): all same-level pairs attract at any d=1 or d=√2, not just Hilbert neighbours.

**Denaturation** (`--denature N`): N steps with weak bonds zeroed before main run.

### Boltzmann validation (`--boltzmann`)

Checks the simulation samples the correct Boltzmann distribution.

- **Pre-enumerated** (default): enumerates conformations by DFS. Suitable for hard-confinement backbone.
- **Live-enumerated** (`--enumerate-live`): builds state set from observations. Required for spring backbone.

Both modes show a Pearson-r convergence plot and observed-vs-predicted frequency scatter.

## Parameter scan

```bash
python scan_assembly.py scan_config.ini
```

Sweeps RED bond energy over a range of values and time points with multiple
repeats per condition, producing a yield-vs-energy plot with error bars.

## Python wrapper flags

### All modes

| Flag | Description |
|---|---|
| `--n0 N` | Target structure size (default 16) |
| `--e1 E` | Bond / backbone confinement energy (default 8.0) |
| `--L N` | Box side length |
| `--nsteps N` | Number of output steps |
| `--nsweep N` | MC sweeps per step |
| `--dens D` | Particle density (alternative to `--L`) |
| `--ncopies N` | Number of target copies |
| `--filehead STR` | Output file prefix (default `hier`) |
| `--no-run` | Skip simulation, re-plot existing files |

### Polymer mode

| Flag | Description |
|---|---|
| `--polymer N` | Polymer size N |
| `--weak-e E` | Mean weak coupling |
| `--weak-std S` | Std dev of weak couplings |
| `--weak-seed N` | Seed for coupling matrix |
| `--hilbert-A A` | Hilbert cardinal neighbour coupling |
| `--hilbert-B B` | Hilbert diagonal neighbour coupling |
| `--hier-red R` | RED coupling (finest level) |
| `--hier-green G` | GREEN coupling (middle level) |
| `--hier-blue B` | BLUE coupling (coarsest level) |
| `--spring-k K` | Spring backbone stiffness |
| `--unspecific` | Non-specific bonding |
| `--denature N` | Pre-equilibration steps |
| `--sim-seed N` | RNG seed |
| `--prob-translate P` | Fraction of translation moves (default 1.0) |
| `--n-neighbours N` | Lattice directions: 4 or 8 (default 4) |
| `--boltzmann` | Enable Boltzmann validation |
| `--enumerate-live` | Live state enumeration (required for spring backbone) |

## Model overview

Each particle is a unit square on a 2D lattice. Pair interactions are computed
up to distance √5:

| Distance | Interaction |
|---|---|
| d < 1 | Hard-core exclusion (+∞) |
| d = 1 | Directional backbone (if set) + wD1 weak coupling |
| d = √2 | Direction-agnostic backbone (if set) + wDsq2 weak coupling |
| d = 2 | wD2 weak coupling only |
| d = √5 | wDsq5 weak coupling only |

VMMC move types: **translation** (one lattice step, cardinal or diagonal) and
**rotation** (90°/180°/270° about the seed particle). The fraction of each is
set by `--prob-translate`.

---

## Repository contents

| File / Directory | Description |
|---|---|
| `makefile` | Build system |
| `run_nucleolus.cpp` | Nucleolus simulation driver |
| `visualize_nucleolus.py` | Nucleolus animation and plotting script |
| `src/NucleolusModel.h/cpp` | Gradient energy model |
| `run_hier.cpp` | Hierarchical assembly driver |
| `run_custom.cpp` | Custom bond driver |
| `run_polymer.cpp` | Polymer driver |
| `run_and_plot.py` | Python wrapper for hierarchical mode |
| `scan_assembly.py` | Parameter scan script |
| `scan_config.ini` | Scan configuration template |
| `src/` | Library source and headers |
| `info.txt` | PhD thesis Chapter 2 (LaTeX source) |
| `input_hier.txt` | Default input file for `run_hier` |
| `old/` | Deprecated files |

## Dependencies

- C++11 compiler (`g++`)
- Python 3 with `numpy` and `matplotlib`

```bash
pip install numpy matplotlib
```
