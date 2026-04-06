# Nucleolus Assembly Simulation

VMMC simulation of nucleolus ribogenesis from Chapter 2 of Sam Whitby's PhD thesis. Four copies of a 16-particle target complex T self-assemble in a gradient-modulated condensate column.

Based on the VMMC framework of Holmes-Cerfon & Wyart (2025), [arXiv:2501.02611](https://arxiv.org/abs/2501.02611).

---

## Build

```bash
mkdir -p obj
make run_nucleolus
```

Requires C++11 (`g++`). Visualisation requires Python 3 with `numpy` and `matplotlib`.

---

## Quick start

```bash
# Assembled free diffusion with gradient, from a file of coupling matrices
./run_nucleolus --width 10 --length 60 --gradient \
    --couplings couplings.txt \
    --free-steps 1000 --denature-steps 500 --steps 10000 \
    --snapshots 200 --seed 1 --output my_run

# Visualise
python3 visualize_nucleolus.py my_run_traj.txt --gradient-length 60 --width 10
```

---

## Physics

### Target complex

16 particles arranged as 4 polymers × 4 segments on a 4×4 grid, partitioned by the n=2 Moore curve:

```
y=3:  1 1 1 2       polymer 0 (ids  0- 3):  (0,0)→(0,1)→(1,1)→(1,2)
y=2:  1 0 2 2       polymer 1 (ids  4- 7):  (0,2)→(0,3)→(1,3)→(2,3)
y=1:  0 0 2 3       polymer 2 (ids  8-11):  (2,1)→(2,2)→(3,2)→(3,3)
y=0:  0 3 3 3       polymer 3 (ids 12-15):  (1,0)→(2,0)→(3,0)→(3,1)
```

### Interactions

**Backbone bonds** (consecutive segments within a polymer) are fixed at strength ≈ 1000 (effectively infinite) and are independent of the gradient. Backbone partners that are at d = 1 or d = √2 experience a bonding energy of −1000; at any other distance the energy is +∞, making separation impossible.

**Weak coupling** is read from a user-supplied file (see [Coupling matrices](#coupling-matrices)). All weak energies are scaled by the chemical gradient γ(xᵢ)·γ(xⱼ). Hard-core overlap (d < 1) is always forbidden.

### Chemical gradient

When `--gradient` is set:

```
γ(x) = min(x / L, 1)
```

γ = 0 at x = 0 (disordered), γ = 1 at x ≥ L (fully coupled). Backbone bonds are never scaled.

### Removal and replacement

After every outer iteration, connected components of the interaction graph are found. A component is **removed and reinserted** near x = 0 if all its particles have x > L and it has no non-backbone interactions with other particles. Removed particles are placed as denatured vertical chains.

---

## Coupling matrices

Weak interactions are stored in four 16×16 symmetric matrices — one per distance shell:

| Section | Distance |
|---|---|
| `[d1]`   | d = 1 |
| `[dsq2]` | d = √2 |
| `[d2]`   | d = 2 |
| `[dsq5]` | d = √5 |

Matrix entry [i][j] is the energy between a particle of local type i and local type j (local type = `global_id % 16`). Positive = attractive, negative = repulsive. Diagonal entries are ignored (hard-core). Lines beginning with `#` are comments.

**`couplings.txt`** is provided with the default values (J = 8, ε = 0.5):
- Same polymer type: −8 at d = 1 and d = √2, −4 at d = 2 (repulsive)
- Cross-type Gō d = 1 neighbours in T: +8 (attractive)
- Cross-type Gō d = √2 neighbours in T: +4 (attractive)
- All other cross-type pairs: 0

Pass a custom file with `--couplings FILE`.

---

## VMMC algorithm

Each outer iteration performs `nParticles` VMMC move attempts. Two move types:

**Translation**: the seed particle and all linked cluster members translate by one lattice step in a random cardinal or diagonal direction.

**Rotation**: the seed and a randomly chosen neighbour define a rotation centre. The cluster (seed + pivot + further recruited particles) rotates by 90°, 180°, or 270° on the lattice.

**Saturated-Link (SL) moves** (`--phi-sl φ`): a fraction φ of moves are SL moves. During an SL move, at most one particle of each of the 16 local types may be in the moving cluster. Once a type slot is filled, all further neighbours of that type are skipped (not recruited, not frustrated). This limits the maximum cluster size to 16 (one complete complex) and prevents kinetic trapping due to same-type repulsion. Backbone bonds are handled by the backbone safety check (see below) rather than by exempting them from SL — if a backbone partner's type slot is already filled, the move is rejected rather than breaking the bond.

**Backbone safety check**: after every cluster proposal (before acceptance), every moving particle checks its post-move position against all neighbours (moving and non-moving) using post-move positions. If any backbone-consecutive pair would end up at a non-bonded distance, the move is aborted.

---

## CLI flags

```
./run_nucleolus [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Phase 3 (main) outer iterations |
| `--free-steps N` | 0 | Phase 1: assembled initial state, no replacement |
| `--denature-steps N` | 0 | Phase 2: all weak bonds = 0, backbone intact, replacement active |
| `--snapshots N` | 1000 | Frames to save across all phases |
| `--length L` | 60 | Column length (lattice units) |
| `--width W` | 10 | Column width (y periodic); must be ≥ 4 |
| `--gradient` | off | Enable chemical gradient γ(x) = x/L |
| `--stokes` | off | Stokes drag: diffusion ∝ 1/cluster radius |
| `--phi-sl φ` | 0.2 | Fraction of Saturated-Link moves (0–1) |
| `--phi-rot φ` | 0.2 | Fraction of rotation moves (90°/180°/270°) |
| `--couplings FILE` | built-in | Coupling matrices file |
| `--output PREFIX` | `nucleolus` | Output prefix for `_traj.txt` and `_stats.txt` |
| `--seed S` | 0 | RNG seed; 0 = time-based |

### Three-phase protocol

1. **Phase 1 — assembled free diffusion** (`--free-steps N`): 4 assembled copies of T placed at x = 0, 6, 12, 18. No replacement.
2. **Phase 2 — denaturation** (`--denature-steps N`): β → 0 (weak bonds zeroed, backbone intact). Replacement active.
3. **Phase 3 — main** (`--steps N`): gradient-coupled dynamics, replacement active.

---

## Output files

**`PREFIX_traj.txt`** — extended XYZ trajectory:
```
<N>
step=S energy=E exited=X L=LL W=WW nCopies=C phase=P
<id> <poly_type> <x> <y> <copy>
...
```

**`PREFIX_stats.txt`** — scalar time series: `step  energy  nExited  acceptRatio`

---

## Visualisation

```bash
python3 visualize_nucleolus.py PREFIX_traj.txt --gradient-length L --width W
```

Opens an animated 3-panel figure: particle positions (left), energy vs step (top right), cumulative exited complexes vs step (bottom right).

```bash
# Save as video
python3 visualize_nucleolus.py PREFIX_traj.txt \
        --gradient-length 60 --width 10 --output my_run.mp4 --fps 10
```

| Flag | Default | Description |
|---|---|---|
| `--output FILE` | display | Save to `.mp4` (ffmpeg) or `.gif` (Pillow) |
| `--fps N` | 5 | Frame rate |
| `--skip N` | 1 | Render every N-th frame |

---

## Dependencies

- C++11 compiler (`g++`)
- Python 3, `numpy`, `matplotlib`
