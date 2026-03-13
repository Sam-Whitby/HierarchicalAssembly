# HierarchicalAssembly

This code reproduces simulations from the following paper by Miranda Holmes-Cerfon and Matthieu Wyart:

> Holmes-Cerfon, M. and Wyart, M., 2025. Hierarchical self-assembly for high-yield addressable complexity at fixed conditions. [arXiv:2501.02611](https://arxiv.org/abs/2501.02611).

The code simulates systems of sticky squares on a square lattice using the Virtual Move Monte Carlo (VMMC) algorithm. The `StickySquares` class is adapted from [vmmc.xyz](http://vmmc.xyz/), extended to handle lattice particles with direction-specific patchy interactions.

This fork adds:
- **`run_custom.cpp`**: driver for arbitrary user-defined bond strengths via a bond file, with omni-directional and diagonal bonding for floppy polymer simulations
- **`run_polymer.cpp`**: driver for polymer mode — auto-generated backbone bonds, multi-distance weak couplings, seeded RNG, and pre-assembled initial configurations
- **`run_and_plot.py`**: Python wrapper to run simulations and visualise results as an animated plot, including Boltzmann distribution validation


---

## Quick start

```bash
# 1. Create the obj directory (once only)
mkdir -p obj

# 2. Compile drivers
make               # builds run_hier
make run_custom    # builds run_custom
make run_polymer   # builds run_polymer

# 3. Run a 4-particle polymer with Boltzmann validation
python run_and_plot.py --polymer 4 --L 12 --nsteps 100000 --nsweep 1 \
  --e1 1000 --weak-e 1.0 --weak-std 0.5 --weak-seed 42 --sim-seed 7
```

---

## Polymer mode (new)

The `--polymer N` flag fully automates polymer simulation:

1. **Auto-generates backbone bonds** following a snake path through the √N × √N target grid, ensuring consecutive chain members are always adjacent (Manhattan distance = 1).
2. **Generates random weak coupling matrices** — one symmetric N×N matrix for each of four spatial distance classes (1, √2, 2, √5), drawn from N(weak_e, weak_std). Entries can be negative (repulsive).
3. **Initialises particles in the assembled configuration** — the chain is pre-placed in a straight line so all backbone bonds are immediately satisfied.
4. **Validates Boltzmann sampling** — after the run, plots the Pearson correlation between observed conformation frequencies and the Boltzmann prediction as a function of simulation time (convergence → 1 implies detailed balance is satisfied).

### Polymer mode example

```bash
python run_and_plot.py --polymer 4 --L 12 --nsteps 1000000 --nsweep 1 \
  --e1 1000 --weak-e 1.0 --weak-std 0.5 --weak-seed 42 --sim-seed 7 --filehead hier
```

This opens three figures:
1. **Simulation animation** — lattice, energy over time, backbone bond matrix
2. **Weak coupling matrices** — 2×2 grid showing the four distance-class coupling matrices
3. **Boltzmann validation** — Pearson r vs simulation steps (left) and observed vs predicted frequency scatter (right)

### Polymer mode flags

| Flag | Description |
|------|-------------|
| `--polymer N` | Polymer size N (must be a perfect square: 4, 9, 16, …). Sets n0=N, 1 copy. |
| `--e1 E` | Backbone bond energy (default: 8.0). Should be large (e.g. 1000) to keep chain assembled. |
| `--weak-e E` | Mean weak coupling energy (default: 1.0). Can produce repulsive couplings when combined with large `--weak-std`. |
| `--weak-std S` | Std dev of weak couplings (default: 0.5). Values > mean allow negative (repulsive) couplings. |
| `--weak-seed N` | RNG seed for coupling matrix generation (for reproducibility). |
| `--sim-seed N` | RNG seed for the C++ VMMC simulation (for reproducibility). |
| `--L N` | Box side length. |
| `--nsteps N` | Number of output steps. |
| `--nsweep N` | MC sweeps per step (default: 1 for fine-grained sampling). |
| `--filehead STR` | Prefix for output files (default: `hier`). |
| `--no-run` | Skip simulation, just re-plot existing output files. |

---

## Running the code manually

### Hierarchical assembly (original)
```bash
make
./run_hier input_hier.txt 8.0
```
where `8.0` is the strongest bond energy (e1). Weaker bonds are set as e1/2, e1/4, e1/8, e1/16 at sub-block boundaries.

### Custom bond strengths
```bash
make run_custom
./run_custom input_hier.txt mybonds.txt
```
where `mybonds.txt` specifies bond strengths between particle pairs (see **Bond file format** below).

### Polymer with weak couplings (new)
```bash
make run_polymer
./run_polymer input_hier.txt mybonds_extended.txt [conffile] [seed]
```
where `mybonds_extended.txt` uses the extended format with `BACKBONE` and `WEAK_D*` sections (see below).

### Python wrapper (recommended)
```bash
python run_and_plot.py [options]
```
Runs the simulation and opens an animated visualisation showing the lattice, energy over time, and the bond-strength coupling matrix.

**Key options (all modes):**

| Flag | Description |
|------|-------------|
| `--n0 N` | Target structure size (must be a perfect square: 4, 16, 64, 256, …) |
| `--p N` | Total particles in simulation (multiple of n0) |
| `--L N` | Box side length (overrides `--dens`) |
| `--nsteps N` | Number of output steps |
| `--nsweep N` | MC sweeps per step |
| `--e1 E` | Bond energy for hierarchical/backbone |
| `--dens D` | Particle density (alternative to `--L`) |
| `--ncopies N` | Number of target copies (alternative to `--p`) |
| `--filehead STR` | Prefix for output files (default: `hier`) |
| `--bond-file PATH` | Use a custom bond file (runs `run_custom`) |
| `--gen-bonds` | Auto-generate a random bond file from N(e1, e1·std) |
| `--bond-std F` | Std dev fraction for `--gen-bonds` (default: 0.3) |
| `--bond-seed N` | Random seed for `--gen-bonds` |
| `--no-run` | Skip simulation, just re-plot existing output files |

Press **spacebar** to pause/resume the animation.

---

## Bond file formats

### Simple format (`run_custom` / `--bond-file`)

One bond per non-comment line:

```
# comment lines start with #
particle_i  particle_j  energy
```

- `particle_i` and `particle_j` must be adjacent (Manhattan distance = 1) in the `l0 × l0` target grid.
- Bonds are stored **omni-directionally**: they fire regardless of runtime orientation.
- Bonds are also active at diagonal distance √2.

### Extended format (`run_polymer`)

Used for polymer mode with multi-distance weak couplings:

```
# Extended polymer bond file
BACKBONE
particle_i  particle_j  energy
...
BACKBONE_END

WEAK_D1
particle_i  particle_j  energy   # weak coupling at distance 1
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

- `BACKBONE` section: strong omni-directional bonds (same rules as simple format). These are the chain backbone bonds.
- `WEAK_D*` sections: symmetric N×N coupling matrices, one entry per unique pair `i j energy` (upper triangle). Added to the pair energy at the corresponding distance. Can be positive (attractive) or negative (repulsive).
- The four supported distance classes: **1** (cardinal), **√2** (diagonal), **2** (two lattice steps), **√5** (knight's move).

### Example bond files

| File | Description |
|------|-------------|
| `bonds_polymer_chain_4.txt` | 3-bond linear polymer, n0=4 |
| `bonds_polymer_chain_16.txt` | 15-bond linear polymer, n0=16 |
| `bonds_polymer_chain.txt` | 63-bond linear polymer, n0=64 |
| `bonds_custom_example.txt` | All-pair random bonds, n0=64 |

### Polymer chain examples

```bash
# Flexible 4-particle polymer (legacy bond file, no weak couplings)
python run_and_plot.py --n0 4 --p 4 --L 12 --nsteps 2000 --nsweep 1 \
  --bond-file bonds_polymer_chain_4.txt --filehead hier --e1 1000

# Flexible 16-particle polymer
python run_and_plot.py --n0 16 --p 16 --L 12 --nsteps 2000 --nsweep 1 \
  --bond-file bonds_polymer_chain_16.txt --filehead hier --e1 1000

# Auto-generated polymer with weak couplings and Boltzmann validation
python run_and_plot.py --polymer 4 --L 12 --nsteps 1000000 --nsweep 1 \
  --e1 1000 --weak-e 1.0 --weak-std 0.5 --weak-seed 42 --sim-seed 7
```

---

## Boltzmann validation

When running with `--polymer`, the code checks that the VMMC simulation correctly samples the Boltzmann distribution.

**States** are defined as distinct polymer conformations, invariant under lattice translation (but not rotation). Each state is represented as the set of relative positions of particles 1…N-1 with respect to particle 0.

**Energy** of each state is computed analytically from the weak coupling matrices: sum of all pairwise coupling energies at their actual separations in that conformation. (The backbone energy is constant for an assembled polymer and cancels in the distribution.)

**Pearson correlation** between observed frequencies and Boltzmann weights (exp(−E)/Z) is plotted as a function of simulation frames. It should converge to ~1 for a correctly-implemented sampler with enough statistics.

---

## Model overview

Each particle is a square on a 2D lattice. Interactions are:

- **Hard-core exclusion**: particles cannot overlap (infinite energy penalty).
- **Backbone bonds** (`run_polymer`): strong omni-directional sticky interactions between designated identity pairs, active at distance 1 and √2. Keeps the polymer chain assembled.
- **Weak couplings** (`run_polymer`): distance-class-specific symmetric coupling matrices, active at distances 1, √2, 2, and √5. Differentiates the energy of distinct conformations.
- **Nearest-neighbour sticky interactions** (`run_hier`/`run_custom`): bond strength depends on sub-block boundary level or is specified per identity pair.

Dynamics use the **Virtual Move Monte Carlo** algorithm (cluster moves), which efficiently samples collective rearrangements including conformational changes of flexible polymers.

---

## Repository contents

| File/Directory | Description |
|----------------|-------------|
| `makefile` | Build system. `make` builds `run_hier`; `make run_custom`/`run_polymer` builds those drivers. Requires `g++` and a pre-existing `obj/` directory. |
| `run_hier.cpp` | Original hierarchical assembly driver |
| `run_custom.cpp` | Custom bond driver with omni-directional + diagonal bonding |
| `run_polymer.cpp` | Polymer driver: extended bond format, weak couplings, initial config, seeded RNG |
| `run_and_plot.py` | Python wrapper: runs simulation, animated visualisation, Boltzmann validation |
| `src/` | All library source and header files |
| `input_hier.txt` | Example input file |
| `bonds_polymer_chain_4.txt` | 4-particle polymer bond file |
| `bonds_polymer_chain_16.txt` | 16-particle polymer bond file |
| `bonds_polymer_chain.txt` | 64-particle polymer bond file |
| `bonds_custom_example.txt` | Random bond matrix example (n0=64) |

### Input file format (`input_hier.txt`)

```
filehead    # output file prefix
n           # target structure size (4, 16, 64, 256, or 1024)
ncopies     # number of copies of the target structure
nsteps      # number of output steps
nsweep      # MC sweeps per step
density     # volume fraction
```

### Output files

- `<filehead>_traj.txt`: particle trajectories (XYZ-like format)
- `<filehead>_stats.txt`: per-step energy and fragment-size histogram
- `<filehead>_polymer_bonds.txt`: generated extended bond file (polymer mode)
- `<filehead>_init.conf`: generated initial configuration file (polymer mode)

---

## Dependencies

- C++11 compiler (`g++`)
- Python 3 with: `numpy`, `matplotlib`

```bash
pip install numpy matplotlib
```
