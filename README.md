# HierarchicalAssembly

This code reproduces simulations from the following paper by Miranda Holmes-Cerfon and Matthieu Wyart:

> Holmes-Cerfon, M. and Wyart, M., 2025. Hierarchical self-assembly for high-yield addressable complexity at fixed conditions. [arXiv:2501.02611](https://arxiv.org/abs/2501.02611).

The code simulates systems of sticky squares on a square lattice using the Virtual Move Monte Carlo (VMMC) algorithm. The `StickySquares` class is adapted from [vmmc.xyz](http://vmmc.xyz/), extended to handle lattice particles with direction-specific patchy interactions.

This fork adds:
- **`run_custom.cpp`**: a new driver for arbitrary user-defined bond strengths via a bond file, with omni-directional and diagonal bonding support for floppy polymer simulations
- **`run_and_plot.py`**: a Python script to run simulations and visualise results as an animated plot
- **Example bond files** for polymer chain simulations


---

## Quick start

```bash
# 1. Create the obj directory (once only)
mkdir -p obj

# 2. Compile the hierarchical assembly driver
make

# 3. Run via Python (runs simulation + shows animated plot)
python run_and_plot.py --n0 16 --p 16 --L 12 --nsteps 2000 --nsweep 1 --e1 8.0
```

---

## Running the code manually

### Hierarchical assembly (original)
```bash
make
./run_hier input_hier.txt 8.0
```
where `8.0` is the strongest bond energy (e1). Weaker bonds are set as e1/2, e1/4, e1/8, e1/16 at sub-block boundaries.

### Custom bond strengths (this fork)
```bash
make run_custom
./run_custom input_hier.txt mybonds.txt
```
where `mybonds.txt` specifies bond strengths between particle pairs (see **Bond file format** below).

### Python wrapper (recommended)
```bash
python run_and_plot.py [options]
```
Runs the simulation and opens an animated visualisation showing the lattice, energy over time, and the bond-strength coupling matrix.

**Key options:**

| Flag | Description |
|------|-------------|
| `--n0 N` | Target structure size (must be a perfect square: 4, 16, 64, 256, …) |
| `--p N` | Total particles in simulation (multiple of n0) |
| `--L N` | Box side length (overrides `--dens`) |
| `--nsteps N` | Number of output steps |
| `--nsweep N` | MC sweeps per step |
| `--e1 E` | Bond energy for hierarchical mode |
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

## Bond file format

Used with `run_custom` / `--bond-file`. One bond per non-comment line:

```
# comment lines start with #
particle_i  particle_j  energy
```

- `particle_i` and `particle_j` are particle indices in `[0, n0)`.
- They must be adjacent (Manhattan distance = 1) in the `l0 × l0` target grid, where `l0 = sqrt(n0)` and particle `p` sits at `col = p % l0`, `row = p // l0`.
- Bonds are stored **omni-directionally**: they fire regardless of which direction the two particles are adjacent at runtime, enabling floppy/flexible polymer conformations.
- Bonds are also active at diagonal distance sqrt(2), allowing conformational changes without bond breaking.

### Example bond files included

| File | Description |
|------|-------------|
| `bonds_polymer_chain_4.txt` | 3-bond linear polymer, n0=4 |
| `bonds_polymer_chain_16.txt` | 15-bond linear polymer, n0=16 |
| `bonds_polymer_chain.txt` | 63-bond linear polymer, n0=64 |
| `bonds_custom_example.txt` | All-pair random bonds, n0=64 |

### Polymer chain example

```bash
python run_and_plot.py --n0 16 --p 16 --L 12 --nsteps 2000 --nsweep 1 \
  --bond-file bonds_polymer_chain_16.txt --filehead hier --e1 1000
```

### Random bond matrix example

```bash
python run_and_plot.py --n0 64 --p 64 --L 12 --nsteps 1000 --nsweep 1 \
  --gen-bonds --bond-std 0.3 --bond-seed 42 --e1 8.0 --filehead hier
```

---

## Model overview

Each particle is a square on a 2D lattice. Interactions are:

- **Hard-core exclusion**: particles cannot overlap (infinite energy penalty).
- **Nearest-neighbour sticky interactions**: energy is released when two bonded-identity particles are adjacent (distance = 1) or diagonally adjacent (distance = sqrt(2), for `run_custom`).
- In `run_hier`, bond strength depends on the sub-block boundary level: e1 > e2=e1/2 > e3=e1/4 > …
- In `run_custom`, bond strengths are specified explicitly per identity pair.

Dynamics use the **Virtual Move Monte Carlo** algorithm (cluster moves), which efficiently samples collective rearrangements.

---

## Repository contents

| File/Directory | Description |
|----------------|-------------|
| `makefile` | Build system. `make` builds `run_hier`; `make run_custom` builds the custom driver. Requires `g++` and a pre-existing `obj/` directory. |
| `run_hier.cpp` | Original hierarchical assembly driver |
| `run_custom.cpp` | Custom bond driver with omni-directional + diagonal bonding |
| `run_and_plot.py` | Python wrapper: runs simulation and shows animated visualisation |
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

---

## Dependencies

- C++11 compiler (`g++`)
- Python 3 with: `numpy`, `matplotlib`

```bash
pip install numpy matplotlib
```
