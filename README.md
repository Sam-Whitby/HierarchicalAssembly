# Nucleolus Assembly Simulation

VMMC simulation of nucleolus ribogenesis from Chapter 2 of Sam Whitby's PhD thesis.
Four copies of a 16-particle target complex T self-assemble in a gradient-modulated condensate column.

Based on the VMMC framework of:
> Holmes-Cerfon, M. and Wyart, M., 2025. Hierarchical self-assembly for high-yield addressable complexity at fixed conditions. [arXiv:2501.02611](https://arxiv.org/abs/2501.02611)

---

## Build

```bash
mkdir -p obj
make run_nucleolus
```

Requires C++11 (`g++`). Visualisation requires Python 3 with `numpy` and `matplotlib`.

---

## Physics

### Target complex

The target T is 16 particles arranged as 4 polymers × 4 segments on a 4×4 grid, partitioned by the n=2 Moore curve:

```
y=3:  1 1 1 2       polymer 0 (ids  0- 3):  (0,0)→(0,1)→(1,1)→(1,2)
y=2:  1 0 2 2       polymer 1 (ids  4- 7):  (0,2)→(0,3)→(1,3)→(2,3)
y=1:  0 0 2 3       polymer 2 (ids  8-11):  (2,1)→(2,2)→(3,2)→(3,3)
y=0:  0 3 3 3       polymer 3 (ids 12-15):  (1,0)→(2,0)→(3,0)→(3,1)
```

Numbers indicate polymer type (0–3). Each polymer has 3 backbone bonds between consecutive segments.

### Interactions (J = 8, ε = 0.5)

Pair energies are computed up to distance √5. All weak couplings are stored in 16×16 matrices indexed by local particle id (0–15).

**Backbone bonds** (consecutive segments within a polymer) are encoded as strong `Triples` (value ≈ 1000) and are **distance-enforced as effectively infinite**:
- At d = 1 or d = √2: energy = −1000 (bonded)
- At any other distance: energy = +∞ (forbidden, always rejected)

**Weak coupling** (all scaled by the chemical gradient γ):

| | d = 1 | d = √2 | d = 2 |
|---|---|---|---|
| Same polymer type | −J = −8 | −J = −8 | −εJ = −4 |
| Cross-type, Gō neighbours in T at d = 1 | +J = +8 | — | — |
| Cross-type, Gō neighbours in T at d = √2 | — | +εJ = +4 | — |
| Other cross-type pairs | 0 | 0 | 0 |

Hard-core overlap (d < 1) is always forbidden (+∞).

<details>
<summary>Full coupling matrices (16×16)</summary>

```
weakD1[i][j]  (d=1)  — same-type: -8, Gō d=1 cross-type: +8, else: 0
     0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
 0:  0   -8   -8   -8    0    0    0    0    0    0    0    0   +8    0    0    0
 1: -8    0   -8   -8   +8    0    0    0    0    0    0    0    0    0    0    0
 2: -8   -8    0   -8    0    0    0    0   +8    0    0    0   +8    0    0    0
 3: -8   -8   -8    0   +8    0   +8    0    0   +8    0    0    0    0    0    0
 4:  0   +8    0   +8    0   -8   -8   -8    0    0    0    0    0    0    0    0
 5:  0    0    0    0   -8    0   -8   -8    0    0    0    0    0    0    0    0
 6:  0    0    0   +8   -8   -8    0   -8    0    0    0    0    0    0    0    0
 7:  0    0    0    0   -8   -8   -8    0    0   +8    0   +8    0    0    0    0
 8:  0    0   +8    0    0    0    0    0    0   -8   -8   -8    0   +8    0   +8
 9:  0    0    0   +8    0    0    0   +8   -8    0   -8   -8    0    0    0    0
10:  0    0    0    0    0    0    0    0   -8   -8    0   -8    0    0    0   +8
11:  0    0    0    0    0    0    0   +8   -8   -8   -8    0    0    0    0    0
12: +8    0   +8    0    0    0    0    0    0    0    0    0    0   -8   -8   -8
13:  0    0    0    0    0    0    0    0   +8    0    0    0   -8    0   -8   -8
14:  0    0    0    0    0    0    0    0    0    0    0    0   -8   -8    0   -8
15:  0    0    0    0    0    0    0    0   +8    0   +8    0   -8   -8   -8    0

weakDsq2[i][j]  (d=√2)  — same-type: -8, Gō d=√2 cross-type: +4, else: 0
     0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
 0:  0   -8   -8   -8    0    0    0    0    0    0    0    0    0    0    0    0
 1: -8    0   -8   -8    0    0    0    0    0    0    0    0   +4    0    0    0
 2: -8   -8    0   -8   +4    0    0    0    0   +4    0    0    0   +4    0    0
 3: -8   -8   -8    0    0   +4    0   +4   +4    0    0    0    0    0    0    0
 4:  0    0   +4    0    0   -8   -8   -8    0    0    0    0    0    0    0    0
 5:  0    0    0   +4   -8    0   -8   -8    0    0    0    0    0    0    0    0
 6:  0    0    0    0   -8   -8    0   -8    0   +4    0    0    0    0    0    0
 7:  0    0    0   +4   -8   -8   -8    0    0    0   +4    0    0    0    0    0
 8:  0    0    0   +4    0    0    0    0    0   -8   -8   -8   +4    0   +4    0
 9:  0    0   +4    0    0    0   +4    0   -8    0   -8   -8    0    0    0   +4
10:  0    0    0    0    0    0    0   +4   -8   -8    0   -8    0    0    0    0
11:  0    0    0    0    0    0    0    0   -8   -8   -8    0    0    0    0    0
12:  0   +4    0    0    0    0    0    0   +4    0    0    0    0   -8   -8   -8
13:  0    0   +4    0    0    0    0    0    0    0    0    0   -8    0   -8   -8
14:  0    0    0    0    0    0    0    0   +4    0    0    0   -8   -8    0   -8
15:  0    0    0    0    0    0    0    0    0   +4    0    0   -8   -8   -8    0

weakD2[i][j]  (d=2)  — same-type: -4, else: 0
     0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
 0:  0   -4   -4   -4    0    0    0    0    0    0    0    0    0    0    0    0
 1: -4    0   -4   -4    0    0    0    0    0    0    0    0    0    0    0    0
 2: -4   -4    0   -4    0    0    0    0    0    0    0    0    0    0    0    0
 3: -4   -4   -4    0    0    0    0    0    0    0    0    0    0    0    0    0
 4:  0    0    0    0    0   -4   -4   -4    0    0    0    0    0    0    0    0
 5:  0    0    0    0   -4    0   -4   -4    0    0    0    0    0    0    0    0
 6:  0    0    0    0   -4   -4    0   -4    0    0    0    0    0    0    0    0
 7:  0    0    0    0   -4   -4   -4    0    0    0    0    0    0    0    0    0
 8:  0    0    0    0    0    0    0    0    0   -4   -4   -4    0    0    0    0
 9:  0    0    0    0    0    0    0    0   -4    0   -4   -4    0    0    0    0
10:  0    0    0    0    0    0    0    0   -4   -4    0   -4    0    0    0    0
11:  0    0    0    0    0    0    0    0   -4   -4   -4    0    0    0    0    0
12:  0    0    0    0    0    0    0    0    0    0    0    0    0   -4   -4   -4
13:  0    0    0    0    0    0    0    0    0    0    0    0   -4    0   -4   -4
14:  0    0    0    0    0    0    0    0    0    0    0    0   -4   -4    0   -4
15:  0    0    0    0    0    0    0    0    0    0    0    0   -4   -4   -4    0
```

</details>

### Chemical gradient

When `--gradient` is set, all weak couplings between particles i and j are scaled by γ(xᵢ)·γ(xⱼ) where:

```
γ(x) = min(x / L, 1)
```

γ = 0 at x = 0 (fully disordered), γ = 1 at x ≥ L (fully coupled). Backbone bonds are never scaled.

### Removal and replacement

After every outer iteration (N_particles VMMC steps), connected components of the interaction graph are identified. A component is **removed and reinserted** near x = 0 if all its particles satisfy x > L and it has no non-backbone interactions with other particles. Components of any size may exit; there is no restriction to exactly one full complex.

---

## Running a simulation

```bash
./run_nucleolus [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Phase 3 (main) outer iterations |
| `--free-steps N` | 0 | Phase 1: assembled initial state, no replacement |
| `--denature-steps N` | 0 | Phase 2: all weak bonds = 0, backbone intact, replacement active |
| `--snapshots N` | 1000 | Frames to save across all phases |
| `--length L` | 60 | Column length (lattice units) |
| `--width W` | 10 | Column width (y periodic) |
| `--gradient` | off | Enable chemical gradient γ(x) = x/L |
| `--stokes` | off | Stokes drag: diffusion ∝ 1/cluster radius |
| `--phi-sl φ` | 0.2 | Fraction of Saturated-Link moves |
| `--phi-rot φ` | 0.2 | Fraction of rotation moves (90°/180°/270°) |
| `--output PREFIX` | `nucleolus` | Output prefix for `_traj.txt` and `_stats.txt` |
| `--seed S` | 0 | RNG seed; 0 = time-based |

### Three-phase protocol

1. **Phase 1 — assembled free diffusion** (`--free-steps N`): 4 assembled copies of T placed at even x-spacings. No replacement. Use to observe diffusion of intact complexes.
2. **Phase 2 — denaturation** (`--denature-steps N`): β → 0 everywhere (weak bonds zeroed, backbone intact). Complexes disassemble and individual polymers are recycled near x = 0.
3. **Phase 3 — main** (`--steps N`): gradient-coupled dynamics, replacement active.

### Example (replicates thesis run)

```bash
./run_nucleolus \
    --width 10 --length 60 --gradient \
    --free-steps 1000 --denature-steps 500 --steps 10000 \
    --snapshots 200 --seed 1 --output my_run
```

---

## Visualisation

```bash
python3 visualize_nucleolus.py my_run_traj.txt --gradient-length 60 --width 10
```

Opens an animated 3-panel figure:
- **Left**: particle positions coloured by polymer type, backbone bonds, gradient overlay
- **Top right**: energy vs step
- **Bottom right**: cumulative exited complexes vs step

Save as video:

```bash
python3 visualize_nucleolus.py my_run_traj.txt \
        --gradient-length 60 --width 10 \
        --output my_run.mp4 --fps 10
```

| Flag | Default | Description |
|---|---|---|
| `--output FILE` | display | Save to `.mp4` (ffmpeg) or `.gif` (Pillow) |
| `--fps N` | 5 | Frame rate |
| `--skip N` | 1 | Render every N-th frame |

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

## Implementation notes

VMMC backbone safety is enforced at two levels:

1. **Energy**: `NucleolusModel::computePairEnergy` returns +∞ for any backbone-consecutive pair not at d = 1 or d = √2, regardless of the chemical gradient.

2. **Move rejection** (in `VMMC::step`): after a cluster move is proposed, every moving particle's post-move position is checked against both its moving and non-moving neighbours using post-move positions for all. This catches the case where the rotation seed (which stays put) has a backbone partner that rotates away: the seed is in the cluster but its post-move position equals its pre-move position, so the energy check correctly detects the separation.

Periodic boundary conditions use full modulo (`floor`-based) wrapping to handle clusters that span the box multiple times via accumulated minimum-image offsets.

---

## Dependencies

- C++11 compiler (`g++`)
- Python 3, `numpy`, `matplotlib`
