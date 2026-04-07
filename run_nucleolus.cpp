/*
  run_nucleolus.cpp

  Nucleolus assembly simulation following Chapter 2 of Sam Whitby's PhD thesis
  "Towards a Model of Annealing in Spatial Gradients".

  Models the assembly of 4 copies of a 16-particle (n=2 Moore-curve) target
  complex T inside a column condensate of width W and length L.  A linear
  chemical gradient γ(x) = min(x/L, 1) scales all weak coupling strengths,
  mimicking denaturing conditions near the condensate core.  When a connected
  cluster of particles is fully past x = L and non-interacting with the rest
  of the system, it is removed and replaced as denatured individual polymers
  near x = 0.

  Coupling matrices follow Equation repulsive_total_eq2 in the thesis:
    - same polymer-type (identical local id / 4): repulsive at d=1,√2,2
    - different polymer-type Gō neighbours in T: attractive at d=1,√2
  All weak couplings are scaled by γ(x_i)·γ(x_j).
  Backbone bonds (hard, value 1000) are unscaled.

  Usage:
    ./run_nucleolus [options]

  Options:
    --steps     N       total outer loop iterations (each = nParticles VMMC moves)  [10000]
    --snapshots N       number of trajectory snapshots to save                      [1000]
    --length    L       condensate length in lattice units                           [60]
    --width     W       column width (periodic y direction)                          [10]
    --gradient          enable linear chemical gradient γ(x) = x/L
    --stokes            enable Stokes hydrodynamic drag (D ∝ 1/R)
    --phi-sl    φ       fraction of Saturated-Link moves                             [0.2]
    --phi-rot   φ       fraction of rotation moves                                   [0.2]
    --output    PREFIX  prefix for output files                                      [nucleolus]
    --seed      S       RNG seed (0 = random)                                        [0]
    --couplings FILE    coupling matrices file (default: built-in J=8, eps=0.5)
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "Demo.h"
#include "VMMC.h"
#include "NucleolusModel.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern double INF;
extern double TOL;

// ============================================================
//  Target complex T  (n=2 Moore curve partitioned into 4 polymers)
//
//  Grid layout (y increases upward, x increases rightward):
//
//    y=3:  1 1 1 2
//    y=2:  1 0 2 2
//    y=1:  0 0 2 3
//    y=0:  0 3 3 3
//
//  Polymer 0 (local ids 0-3):   (0,0)→(0,1)→(1,1)→(1,2)   left S-shape
//  Polymer 1 (local ids 4-7):   (0,2)→(0,3)→(1,3)→(2,3)   top-left L
//  Polymer 2 (local ids 8-11):  (2,1)→(2,2)→(3,2)→(3,3)   right S-shape
//  Polymer 3 (local ids 12-15): (1,0)→(2,0)→(3,0)→(3,1)   bottom-right L
//
//  Backbone (within each polymer, consecutive segment pairs):
//    (0,1),(1,2),(2,3) | (4,5),(5,6),(6,7) | (8,9),(9,10),(10,11) | (12,13),(13,14),(14,15)
// ============================================================

static const int N0 = 16;        // particles per complex
static const int N_POLYMER = 4;  // polymers per complex
static const int N_SEG     = 4;  // segments per polymer

// Target positions for local ids 0..15
static const int TARGET_X[N0] = { 0,0,1,1,  0,0,1,2,  2,2,3,3,  1,2,3,3 };
static const int TARGET_Y[N0] = { 0,1,1,2,  2,3,3,3,  1,2,2,3,  0,0,0,1 };

// Backbone consecutive pairs within each polymer (local ids)
static const int BACKBONE_PAIRS[][2] = {
    {0,1},{1,2},{2,3},
    {4,5},{5,6},{6,7},
    {8,9},{9,10},{10,11},
    {12,13},{13,14},{14,15}
};
static const int N_BB_PAIRS = 12;

// Polymer type for a local id: id / N_SEG
inline int polyType(int id) { return id / N_SEG; }

// Squared distance between two target-complex positions
inline double targetDistSqd(int id1, int id2) {
    double dx = TARGET_X[id1] - TARGET_X[id2];
    double dy = TARGET_Y[id1] - TARGET_Y[id2];
    return dx*dx + dy*dy;
}

// ============================================================
//  Load coupling matrices from a text file.
//
//  Format: four named sections [d1] [dsq2] [d2] [dsq5], each
//  followed by N0 lines of N0 whitespace-separated values.
//  Lines beginning with '#' are comments and are ignored.
//  Matrices must be symmetric; a warning is printed if not.
// ============================================================
static void loadCouplingMatrices(
    const string& path,
    vector<vector<double>>& wD1,
    vector<vector<double>>& wDsq2,
    vector<vector<double>>& wD2,
    vector<vector<double>>& wDsq5)
{
    wD1.assign(N0, vector<double>(N0, 0.0));
    wDsq2.assign(N0, vector<double>(N0, 0.0));
    wD2.assign(N0, vector<double>(N0, 0.0));
    wDsq5.assign(N0, vector<double>(N0, 0.0));

    map<string, vector<vector<double>>*> sections = {
        {"d1", &wD1}, {"dsq2", &wDsq2}, {"d2", &wD2}, {"dsq5", &wDsq5}
    };

    ifstream f(path);
    if (!f) {
        cerr << "[ERROR] Cannot open couplings file: " << path << "\n";
        exit(EXIT_FAILURE);
    }

    vector<vector<double>>* current = nullptr;
    int rowsRead = 0;
    string sectionName;
    string line;
    int lineNum = 0;

    while (getline(f, line)) {
        lineNum++;
        // Strip comments
        auto hash = line.find('#');
        if (hash != string::npos) line = line.substr(0, hash);
        // Trim
        size_t s = line.find_first_not_of(" \t\r\n");
        if (s == string::npos) continue;
        line = line.substr(s);

        // Section header?
        if (line[0] == '[') {
            size_t close = line.find(']');
            if (close == string::npos) {
                cerr << "[ERROR] couplings file line " << lineNum << ": malformed section header\n";
                exit(EXIT_FAILURE);
            }
            sectionName = line.substr(1, close - 1);
            auto it = sections.find(sectionName);
            if (it == sections.end()) {
                cerr << "[ERROR] couplings file: unknown section [" << sectionName << "]\n";
                exit(EXIT_FAILURE);
            }
            current = it->second;
            rowsRead = 0;
            continue;
        }

        // Data row
        if (!current) {
            cerr << "[ERROR] couplings file line " << lineNum << ": data outside any section\n";
            exit(EXIT_FAILURE);
        }
        if (rowsRead >= N0) {
            cerr << "[ERROR] couplings file [" << sectionName << "]: more than " << N0 << " data rows\n";
            exit(EXIT_FAILURE);
        }
        istringstream ss(line);
        for (int j = 0; j < N0; j++) {
            if (!(ss >> (*current)[rowsRead][j])) {
                cerr << "[ERROR] couplings file [" << sectionName << "] row " << rowsRead
                     << ": expected " << N0 << " values\n";
                exit(EXIT_FAILURE);
            }
        }
        rowsRead++;
    }

    // Verify all four sections were read
    for (auto& kv : sections) {
        // count non-zero or just check rowsRead somehow — we can't, so trust the parse above
        (void)kv;
    }

    // Symmetry check
    auto checkSym = [&](const vector<vector<double>>& m, const string& name) {
        for (int i = 0; i < N0; i++)
            for (int j = 0; j < N0; j++)
                if (fabs(m[i][j] - m[j][i]) > 1e-9)
                    cerr << "[WARNING] couplings [" << name << "]: not symmetric at ("
                         << i << "," << j << "): " << m[i][j] << " vs " << m[j][i] << "\n";
    };
    checkSym(wD1, "d1"); checkSym(wDsq2, "dsq2");
    checkSym(wD2, "d2"); checkSym(wDsq5, "dsq5");
}

// ============================================================
//  Build default coupling matrices (16×16, indexed by local id)
//  J=8, eps=0.5: same-polymer repulsion, cross-type Gō attraction.
// ============================================================
static void buildCouplingMatrices(
    double J, double eps,
    vector<vector<double>>& wD1,
    vector<vector<double>>& wDsq2,
    vector<vector<double>>& wD2,
    vector<vector<double>>& wDsq5)
{
    wD1.assign(N0, vector<double>(N0, 0.0));
    wDsq2.assign(N0, vector<double>(N0, 0.0));
    wD2.assign(N0, vector<double>(N0, 0.0));
    wDsq5.assign(N0, vector<double>(N0, 0.0));

    for (int i = 0; i < N0; i++) {
        for (int j = 0; j < N0; j++) {
            if (i == j) continue;

            bool sameType = (polyType(i) == polyType(j));
            double dsqd = targetDistSqd(i, j);

            if (sameType) {
                wD1[i][j]   = -J;
                wDsq2[i][j] = -J;
                wD2[i][j]   = -eps * J;
            } else {
                if (dsqd < 1.0 + TOL)      wD1[i][j]   = J;
                else if (dsqd < 2.0 + TOL) wDsq2[i][j] = eps * J;
            }
        }
    }
}

// ============================================================
//  Build backbone Triples for all nCopies copies
//  Each backbone pair (p,q) gets entries in both east[] and north[]
//  (direction-agnostic) with value bbEnergy ≈ 1000.
// ============================================================
static void buildBackboneTriples(int nCopies, double bbEnergy,
                                  vector<Triple>& north, vector<Triple>& east)
{
    north.clear();
    east.clear();
    for (int c = 0; c < nCopies; c++) {
        int base = c * N0;
        for (int k = 0; k < N_BB_PAIRS; k++) {
            int gi = base + BACKBONE_PAIRS[k][0];
            int gj = base + BACKBONE_PAIRS[k][1];
            // Add both orientations and both arrays for direction-agnostic lookup
            east.push_back({gi, gj, bbEnergy});
            east.push_back({gj, gi, bbEnergy});
            north.push_back({gi, gj, bbEnergy});
            north.push_back({gj, gi, bbEnergy});
        }
    }
}

// ============================================================
//  Place particles as denatured (linear) polymers near x = 0.
//  4 copies × 4 polymers each = 16 polymers of 4 particles.
//  Each polymer is placed as a horizontal chain.
//  poly_global = copy * N_POLYMER + poly_within_copy (0..15)
//  x_base = 1 + (poly_global / W) * (N_SEG + 1)   [gap of 1 between groups]
//  y      = poly_global % W
// ============================================================
static void placeParticles(vector<Particle>& particles,
                            CellList& cells, Box& box,
                            int nCopies, int W)
{
    cells.reset();
    int nParticles = nCopies * N0;
    for (int i = 0; i < nParticles; i++) {
        particles[i].index = i;
        particles[i].position.resize(2);
        particles[i].orientation.resize(2);
        particles[i].orientation[0] = 1.0;
        particles[i].orientation[1] = 0.0;
    }

    int nPolymers = nCopies * N_POLYMER;  // 16
    for (int poly = 0; poly < nPolymers; poly++) {
        int copy = poly / N_POLYMER;
        int p    = poly % N_POLYMER;
        // x starting position for this polymer
        int x_base = 1 + (poly / W) * (N_SEG + 1);
        int y      = poly % W;

        for (int s = 0; s < N_SEG; s++) {
            int global_id = copy * N0 + p * N_SEG + s;
            double px = (double)(x_base + s);
            double py = (double)(y % W);
            particles[global_id].position[0] = px;
            particles[global_id].position[1] = py;
            // Apply PBC (y-periodic, x capped by large box)
            box.periodicBoundaries(particles[global_id].position);
            // Register in cell list
            particles[global_id].cell = cells.getCell(particles[global_id]);
            cells.initCell(particles[global_id].cell, particles[global_id]);
        }
    }
}

// ============================================================
//  Place particles as 4 assembled copies of the target complex T.
//  Each copy is placed at a different x offset within the condensate.
//  Copy c is placed at x_offset = 2 + c*6 (copies span x = 2..5, 8..11, 14..17, 20..23).
//  The complex is centered in y: y_offset = (W-4)/2.
// ============================================================
static void placeAssembled(vector<Particle>& particles,
                            CellList& cells, Box& box,
                            int nCopies, int W, int /*L_col*/)
{
    cells.reset();
    int nParticles = nCopies * N0;
    for (int i = 0; i < nParticles; i++) {
        particles[i].index = i;
        particles[i].position.resize(2);
        particles[i].orientation.resize(2);
        particles[i].orientation[0] = 1.0;
        particles[i].orientation[1] = 0.0;
    }

    // Place copies in a row: copy c at x = c*(N_SEG+1), y = 0.
    // Each 4×4 complex occupies x ∈ [c*5, c*5+3] — separated by a gap of 1
    // in x, so no two copies ever share a lattice site or have overlapping
    // interaction ranges (min inter-copy distance = 1 > sqrt(5)/2 would be
    // too close, but gap=1 means nearest cross-copy distance = 1 which is
    // allowed — they are NOT backbone-connected so they interact normally).
    // A gap of 2 guarantees no direct d=1 contact between copies; use gap=2.
    for (int c = 0; c < nCopies; c++) {
        int x_offset = c * (N_SEG + 2);  // gap of 2 between 4-wide complexes → stride=6
        for (int lid = 0; lid < N0; lid++) {
            int gi = c * N0 + lid;
            double px = (double)(x_offset + TARGET_X[lid]);
            double py = (double)(TARGET_Y[lid] % W);
            particles[gi].position[0] = px;
            particles[gi].position[1] = py;
            box.periodicBoundaries(particles[gi].position);
            particles[gi].cell = cells.getCell(particles[gi]);
            cells.initCell(particles[gi].cell, particles[gi]);
        }
    }
}

// ============================================================
//  Occupancy map: set<pair<int,int>> of occupied (x,y) sites
// ============================================================
static set<pair<int,int>> buildOccupancy(const vector<Particle>& particles)
{
    set<pair<int,int>> occ;
    for (const auto& p : particles) {
        int px = (int)round(p.position[0]);
        int py = (int)round(p.position[1]);
        occ.insert({px, py});
    }
    return occ;
}

// ============================================================
//  Try to place a set of particle indices as vertical polymer chains
//  near x = 0, scanning column-major (y varies fastest).
//
//  Placement layout (example: W=10, N_SEG=4):
//    polymer 0 → x=0, y=0,1,2,3
//    polymer 1 → x=0, y=4,5,6,7
//    polymer 2 → x=0, y=8,9 then x=1, y=0,1   (PBC wrap in y)
//    polymer 3 → x=1, y=2,3,4,5
//
//  Consecutive segments within a polymer are connected by backbone bonds
//  at distance 1 (same x) or sqrt(2) (crossing x boundary with PBC wrap).
//  Returns true if successful, false if no space was found.
// ============================================================
static bool replacementPlacement(
    vector<Particle>& particles, CellList& cells, Box& box,
    const vector<int>& globalIds, int W, int L_col,
    vmmc::VMMC& vmmc)
{
    int nParts = (int)globalIds.size();

    // Build occupancy map excluding the particles to be replaced
    set<int> replSet(globalIds.begin(), globalIds.end());
    set<pair<int,int>> occ;
    for (int i = 0; i < (int)particles.size(); i++) {
        if (replSet.count(i)) continue;
        int px = (int)round(particles[i].position[0]);
        int py = (int)round(particles[i].position[1]);
        occ.insert({px, py});
    }

    // Scan linearised column (linear index = x*W + y) from x=0 forward.
    // Find a run of nParts consecutive free sites.
    int maxLinear = (L_col + 1) * W;  // search up to x = L_col

    for (int base = 0; base < maxLinear; base++) {
        bool allFree = true;
        for (int k = 0; k < nParts && allFree; k++) {
            int idx = base + k;
            int x   = idx / W;
            int y   = idx % W;
            if (x > L_col) { allFree = false; break; }
            if (occ.count({x, y})) allFree = false;
        }
        if (!allFree) continue;

        // Place particles in column-major order
        for (int k = 0; k < nParts; k++) {
            int gi  = globalIds[k];
            int idx = base + k;
            double newX = (double)(idx / W);
            double newY = (double)(idx % W);
            particles[gi].position[0] = newX;
            particles[gi].position[1] = newY;
            box.periodicBoundaries(particles[gi].position);
            int newCell = cells.getCell(particles[gi]);
            cells.updateCell(newCell, particles[gi], particles);
            // Sync VMMC's internal preMovePosition to prevent stale-position desync:
            // without this, rotation moves compute cluster displacements from the old
            // (large-x) position, producing y values out of range.
            vmmc.syncPosition(gi, &particles[gi].position[0]);
        }
        return true;
    }
    return false;  // no space available
}

// ============================================================
//  Interaction-graph connected components using computePairEnergy.
//  Builds adjacency for ALL particles (not just within x>L).
//  Returns fragmentID[i] = component id for particle i,
//  and a list of component members.
// ============================================================
static int buildComponents(NucleolusModel& model,
                            vector<Particle>& particles,
                            int nParticles,
                            vector<int>& fragmentID,
                            vector<vector<int>>& components)
{
    fragmentID.assign(nParticles, -1);
    components.clear();
    int nfrag = 0;

    const int maxInt = 30;
    unsigned int nbrs[maxInt];

    for (int i = 0; i < nParticles; i++) {
        if (fragmentID[i] != -1) continue;

        vector<int> comp = {i};
        fragmentID[i] = nfrag;

        for (int ci = 0; ci < (int)comp.size(); ci++) {
            int j = comp[ci];
            int nn = (int)model.computeInteractions(
                j, &particles[j].position[0], &particles[j].orientation[0], nbrs);
            for (int k = 0; k < nn; k++) {
                int nbr = (int)nbrs[k];
                if (fragmentID[nbr] != -1) continue;
                double e = model.computePairEnergy(
                    j,   &particles[j].position[0],   &particles[j].orientation[0],
                    nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
                if (e != 0.0 && e < 1e5) {
                    fragmentID[nbr] = nfrag;
                    comp.push_back(nbr);
                }
            }
        }
        components.push_back(comp);
        nfrag++;
    }
    return nfrag;
}


// ============================================================
//  Check whether a component that has exited is a perfect copy of target T.
//
//  Criteria:
//    1. Exactly N0 = 16 particles, one of each local type 0..15.
//    2. Internal energy (gradient off) equals targetEnergy within ±0.5.
//       The energy matching guarantees all particles are at the correct
//       relative positions with all expected bonds present.
// ============================================================
static bool isPerfectComplex(NucleolusModel& model,
                              vector<Particle>& particles,
                              const vector<int>& comp,
                              double targetEnergy)
{
    if ((int)comp.size() != N0) return false;

    // One of each local type.
    bool typePresent[N0] = {};
    for (int gi : comp) {
        int lid = gi % N0;
        if (typePresent[lid]) return false;
        typePresent[lid] = true;
    }
    for (int t = 0; t < N0; t++)
        if (!typePresent[t]) return false;

    // Sum internal pair energies with gradient disabled.
    bool savedGradient  = model.hasGradient;
    bool savedDenatured = model.denatured;
    model.hasGradient  = false;
    model.denatured    = false;

    double energy = 0.0;
    bool hardCore = false;
    for (int a = 0; a < (int)comp.size() && !hardCore; a++) {
        for (int b = a + 1; b < (int)comp.size() && !hardCore; b++) {
            int i = comp[a], j = comp[b];
            double e = model.computePairEnergy(
                i, &particles[i].position[0], &particles[i].orientation[0],
                j, &particles[j].position[0], &particles[j].orientation[0]);
            if (e > 1e5) { hardCore = true; break; }
            energy += e;
        }
    }

    model.hasGradient  = savedGradient;
    model.denatured    = savedDenatured;

    if (hardCore) return false;
    return (fabs(energy - targetEnergy) < 0.5);
}


// ============================================================
//  Check exit condition and perform remove/replace.
//
//  exitedMass       incremented by comp.size() for every removed component.
//  exitedPerfectMass incremented by N0 for every perfect complex removed.
//  fp_exits         receives one log line per removed component (may be NULL).
//  step             current simulation step (for the log).
//
//  Returns number of components removed this step.
// ============================================================
static int checkAndReplace(NucleolusModel& model,
                            vector<Particle>& particles,
                            CellList& cells, Box& box,
                            int nCopies, int W, int L_col,
                            vmmc::VMMC& vmmc,
                            double targetComplexEnergy,
                            long long& exitedMass,
                            long long& exitedPerfectMass,
                            FILE* fp_exits,
                            long long step)
{
    int nParticles = nCopies * N0;
    vector<int> fragmentID;
    vector<vector<int>> components;
    buildComponents(model, particles, nParticles, fragmentID, components);

    int nReplaced = 0;
    for (auto& comp : components) {
        // Check if ALL particles in this component have x > L_col
        bool allPast = true;
        for (int gi : comp)
            if (particles[gi].position[0] <= (double)L_col) { allPast = false; break; }
        if (!allPast) continue;

        // Verify isolation: no non-backbone edges to particles outside this component
        set<int> compSet(comp.begin(), comp.end());
        bool isolated = true;
        const int maxInt = 30;
        unsigned int nbrs[maxInt];
        for (int gi : comp) {
            int nn = (int)model.computeInteractions(
                gi, &particles[gi].position[0], &particles[gi].orientation[0], nbrs);
            for (int k = 0; k < nn; k++) {
                int nbr = (int)nbrs[k];
                if (!compSet.count(nbr)) {
                    double e = model.computePairEnergy(
                        gi,  &particles[gi].position[0],  &particles[gi].orientation[0],
                        nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
                    if (e != 0.0 && e < 1e5) { isolated = false; break; }
                }
            }
            if (!isolated) break;
        }
        if (!isolated) continue;

        // Sort by global index so polymer order is preserved
        sort(comp.begin(), comp.end());

        // Classify and log BEFORE placement overwrites positions.
        bool perfect = isPerfectComplex(model, particles, comp, targetComplexEnergy);

        if (fp_exits) {
            fprintf(fp_exits, "%lld\t%d\t%d", step, (int)comp.size(), perfect ? 1 : 0);
            for (int gi : comp)
                fprintf(fp_exits, "\t%d:%d:%.4f:%.4f",
                        gi, gi % N0,
                        particles[gi].position[0],
                        particles[gi].position[1]);
            fprintf(fp_exits, "\n");
        }

        if (replacementPlacement(particles, cells, box, comp, W, L_col, vmmc)) {
            nReplaced++;
            exitedMass        += (long long)comp.size();
            exitedPerfectMass += perfect ? (long long)N0 : 0LL;
        }
    }
    return nReplaced;
}

// ============================================================
//  Write one frame to trajectory file (extended XYZ format).
//  Header line 2 carries step, energy, cumulative exit count, and box params
//  so that the visualizer can plot scalar time-series without a separate file.
//  Columns: particle_id  polymer_type  x  y  copy
//
//  exitedMass:    cumulative count of all particles that have exited (each = 1).
//  exitedPerfect: subset of exitedMass that exited inside a perfect complex.
// ============================================================
static void writeFrame(FILE* fp, const vector<Particle>& particles,
                        int nCopies, double L_col, double W,
                        long long step, double energy,
                        long long totalExited,
                        long long exitedMass,
                        long long exitedPerfect,
                        const char* phase = "main")
{
    int nParticles = (int)particles.size();
    fprintf(fp, "%d\n", nParticles);
    fprintf(fp,
        "step=%lld energy=%.6f exited=%lld exitedMass=%lld exitedPerfect=%lld"
        " L=%.1f W=%.1f nCopies=%d phase=%s\n",
        step, energy, totalExited, exitedMass, exitedPerfect,
        L_col, W, nCopies, phase);
    for (int i = 0; i < nParticles; i++) {
        int copy  = i / N0;
        int lid   = i % N0;
        int ptype = polyType(lid);
        fprintf(fp, "%d %d %.4f %.4f %d\n",
                i, ptype,
                particles[i].position[0],
                particles[i].position[1],
                copy);
    }
}

// ============================================================
//  main()
// ============================================================
int main(int argc, char** argv)
{
    // --- Defaults ---
    long long nsteps         = 10000;
    long long nsnaps         = 1000;
    long long freeSteps      = 0;      // phase 1: assembled free diffusion
    long long denatureSteps  = 0;      // phase 2: β→0 denaturation
    int  L_col               = 60;
    int  W                   = 10;
    bool useGradient         = false;
    bool useStokes           = false;
    double phi_sl            = 0.2;
    double phi_rot           = 0.2;
    string outPrefix         = "nucleolus";
    unsigned int seed        = 0;
    string couplingsFile     = "";

    // --- Parse arguments ---
    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"--steps")          && i+1<argc) { nsteps        = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--snapshots")      && i+1<argc) { nsnaps        = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--free-steps")     && i+1<argc) { freeSteps     = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--denature-steps") && i+1<argc) { denatureSteps = atoll(argv[++i]); }
        else if (!strcmp(argv[i],"--length")         && i+1<argc) { L_col         = atoi(argv[++i]); }
        else if (!strcmp(argv[i],"--width")          && i+1<argc) { W             = atoi(argv[++i]); }
        else if (!strcmp(argv[i],"--gradient"))                    { useGradient   = true; }
        else if (!strcmp(argv[i],"--stokes"))                      { useStokes     = true; }
        else if (!strcmp(argv[i],"--phi-sl")         && i+1<argc) { phi_sl        = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--phi-rot")        && i+1<argc) { phi_rot       = atof(argv[++i]); }
        else if (!strcmp(argv[i],"--output")         && i+1<argc) { outPrefix     = argv[++i]; }
        else if (!strcmp(argv[i],"--seed")           && i+1<argc) { seed          = (unsigned int)atoi(argv[++i]); }
        else if (!strcmp(argv[i],"--couplings")      && i+1<argc) { couplingsFile = argv[++i]; }
        else {
            cerr << "Unknown argument: " << argv[i] << "\n"
                 << "Run ./run_nucleolus --help for usage.\n";
        }
    }

    // Target complex is 4 particles tall; column must be at least that wide.
    if (W < N_SEG) {
        cerr << "[ERROR] --width " << W << " is less than " << N_SEG
             << " (the height of one polymer in the target complex). "
             << "A complex can never exist in this column. Increase --width.\n";
        return 1;
    }

    long long totalSteps = freeSteps + denatureSteps + nsteps;
    // Clamp snapshots; distribute evenly over total simulation
    if (nsnaps > totalSteps + 1) nsnaps = totalSteps + 1;
    long long saveEvery = (nsnaps <= 1) ? totalSteps : max(1LL, totalSteps / (nsnaps - 1));

    cout << "=== Nucleolus Assembly Simulation ===" << endl;
    cout << "  free-steps=" << freeSteps
         << "  denature-steps=" << denatureSteps
         << "  main-steps=" << nsteps
         << "  snapshots=" << nsnaps << " (every " << saveEvery << " steps)" << endl;
    cout << "  L=" << L_col << " W=" << W
         << "  gradient=" << useGradient << " stokes=" << useStokes
         << "  phi_sl=" << phi_sl << " phi_rot=" << phi_rot << endl;

    // --- Parameters ---
    const int    nCopies    = 4;
    const int    nParticles = nCopies * N0;
    const double J          = 8.0;
    const double eps        = 0.5;
    const double bbEnergy   = 1000.0;   // backbone bond strength (effective ∞)
    const double X_MAX      = (double)(max(5 * L_col, 300)); // effective box x size

    const unsigned int dimension       = 2;
    const double interactionRange      = 2.5;  // covers d=1,√2,2,√5
    // maxInteractions must be large enough for the densest possible neighbourhood.
    // With a narrow column (small W), the cell-list neighbourhood can cover the
    // entire y-extent, so all nParticles-1 other particles may appear as candidates.
    const unsigned int maxInteractions = (unsigned int)(nParticles);
    const double interactionEnergy     = 0.0;
    const bool isLattice               = true;

    // --- Build coupling matrices ---
    vector<vector<double>> wD1, wDsq2, wD2, wDsq5;
    if (!couplingsFile.empty()) {
        cout << "Loading coupling matrices from: " << couplingsFile << endl;
        loadCouplingMatrices(couplingsFile, wD1, wDsq2, wD2, wDsq5);
    } else {
        buildCouplingMatrices(J, eps, wD1, wDsq2, wD2, wDsq5);
    }

    // --- Build backbone Triples ---
    vector<Triple> north0, east0;
    buildBackboneTriples(nCopies, bbEnergy, north0, east0);

    // Empty north/east0 used only as starting point; buildBackboneTriples
    // already fills for all copies.
    vector<Triple> north = north0;
    vector<Triple> east  = east0;

    // No spring backbone (harmonic spring disabled; use hard backbone above)
    vector<vector<int>> bbPartners(N0);  // empty: no harmonic spring
    double springK = 0.0;

    Interactions interactions(nParticles, N0, north, east,
                               wD1, wDsq2, wD2, wDsq5,
                               springK, bbPartners);

    // --- Simulation box ---
    // x: non-periodic (hard wall at x=0 via boundary callback)
    // y: periodic with period W
    vector<double> boxSize = { X_MAX, (double)W };
    vector<bool>   isPeriodic = { false, true };
    Box box(boxSize, isPeriodic);
    box.isLattice = true;

    // --- Cell list ---
    CellList cells;
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);
    // With a narrow column the cell list may have only 3 y-cells, so the
    // 3x3 neighbourhood covers the whole column and a single cell can hold
    // all nParticles.  Bump maxParticles and resize every cell buffer so
    // the overflow check never fires spuriously.
    if (cells.maxParticles < (unsigned int)nParticles)
        cells.maxParticles = (unsigned int)nParticles;
    for (auto& cell : cells)
        if (cell.particles.size() < cells.maxParticles)
            cell.particles.resize(cells.maxParticles);

    // --- Particles ---
    vector<Particle> particles(nParticles);

    // --- Nucleolus model (with gradient) ---
    NucleolusModel model(box, particles, cells,
                          maxInteractions, interactionEnergy, interactionRange,
                          interactions, (double)L_col, useGradient);

    // --- Compute energy baseline from the fully assembled state ---
    // Baseline is computed with gradient OFF (γ=1 everywhere) so that E=0 at
    // step 0 regardless of where the assembled complexes sit in x.
    placeAssembled(particles, cells, box, nCopies, W, L_col);
    model.hasGradient = false;
    double baselineEnergy = model.getEnergy() * nParticles;
    model.hasGradient = useGradient;  // restore to requested setting
    // Energy of a single perfectly-assembled complex (gradient off).
    // Used by isPerfectComplex() to decide whether an exiting structure
    // has all particles in the correct relative positions.
    double targetComplexEnergy = baselineEnergy / (double)nCopies;

    // Now set the actual initial state for the simulation
    if (!(freeSteps > 0 || denatureSteps > 0)) {
        // No assembled or denaturation phase: start from denatured linear chains
        placeParticles(particles, cells, box, nCopies, W);
    }
    // else: keep the assembled state for phase 1

    // --- VMMC setup ---
    vector<double> coordinates(dimension * nParticles);
    vector<double> orientations_v(dimension * nParticles);
    // Note: std::vector<bool> is specialised and lacks .data(); use a plain array.
    bool isIsotropicArr[nParticles];
    for (int i = 0; i < nParticles; i++) {
        coordinates[2*i]       = particles[i].position[0];
        coordinates[2*i + 1]   = particles[i].position[1];
        orientations_v[2*i]     = 1.0;
        orientations_v[2*i + 1] = 0.0;
        isIsotropicArr[i]       = true;
    }

    double maxTrialTranslation = 1.5;
    double maxTrialRotation    = (phi_rot > 0.0) ? M_PI : 0.0;
    double probTranslate       = 1.0 - phi_rot;
    double referenceRadius     = 0.5;
    bool   isRepulsive         = false;
    int    nLatticeNeighbours  = 8;   // 8 directions (including diagonals)

    using namespace std::placeholders;
    vmmc::CallbackFunctions callbacks;
    callbacks.energyCallback =
        std::bind(&NucleolusModel::computeEnergy, &model, _1, _2, _3);
    callbacks.pairEnergyCallback =
        std::bind(&NucleolusModel::computePairEnergy, &model, _1, _2, _3, _4, _5, _6);
    callbacks.interactionsCallback =
        std::bind(&NucleolusModel::computeInteractions, &model, _1, _2, _3, _4);
    callbacks.postMoveCallback =
        std::bind(&NucleolusModel::applyPostMoveUpdates, &model, _1, _2, _3);

    // Hard wall at x = 0: reject any move that takes a particle to x < 0
    callbacks.boundaryCallback =
        [](unsigned int /*idx*/, const double* pos, const double* /*ori*/) -> bool {
            return pos[0] < -0.5;  // x = -1 on lattice → outside boundary
        };

    vmmc::VMMC vmmc(nParticles, dimension, coordinates.data(), orientations_v.data(),
                     maxTrialTranslation, maxTrialRotation,
                     probTranslate, referenceRadius,
                     maxInteractions, &boxSize[0], isIsotropicArr, isRepulsive,
                     callbacks, isLattice, nLatticeNeighbours,
                     phi_sl, N0 /*slN0: particle type period*/);

    // Stokes: set hydrAlpha = 0 to disable (unit diffusion), 1 to enable
    vmmc.hydrAlpha = useStokes ? 1.0 : 0.0;

    // Set RNG seed (seed=0 means use the time-based default from the constructor)
    if (seed != 0) vmmc.rng.setSeed(seed);
    cout << "RNG seed: " << vmmc.rng.getSeed()
         << (seed == 0 ? "  (random; use --seed N for reproducibility)" : "") << endl;

    // --- Open output files ---
    string trajFile  = outPrefix + "_traj.txt";
    string statFile  = outPrefix + "_stats.txt";
    string exitsFile = outPrefix + "_exits.txt";

    FILE* fp_traj  = fopen(trajFile.c_str(), "w");
    FILE* fp_stat  = fopen(statFile.c_str(), "w");
    FILE* fp_exits = fopen(exitsFile.c_str(), "w");
    if (!fp_traj || !fp_stat || !fp_exits) {
        cerr << "Cannot open output files.\n";
        return 1;
    }
    // Exits log: one tab-separated line per exit event.
    // Columns: step  nParticles  isPerfect  id:lid:x:y (repeated)
    fprintf(fp_exits, "# step\tnParticles\tisPerfect\tparticles(id:lid:x:y)...\n");
    fprintf(fp_stat, "# step  energy  nExited  exitedMass  exitedPerfect  acceptRatio\n");

    // --- Simulation loop (three phases) ---
    cout << "Starting simulation..." << endl;
    clock_t startTime         = clock();
    long long totalExited     = 0;
    long long totalExitedMass = 0;    // cumulative particles that have exited
    long long totalExitedPerfect = 0; // subset from perfect complexes
    long long globalStep      = 0;

    // Write initial frame (step 0)
    double initEnergy = model.getEnergy() * nParticles - baselineEnergy;
    const char* initPhase = (freeSteps > 0) ? "assembled" : (denatureSteps > 0 ? "denature" : "main");
    writeFrame(fp_traj, particles, nCopies, L_col, W, 0, initEnergy, 0, 0, 0, initPhase);
    fprintf(fp_stat, "0  %.4f  0  0  0  0.0000\n", initEnergy);

    // Helper lambda: run one outer iteration and optionally save a frame.
    auto runStep = [&](const char* phase, bool doReplace) {
        globalStep++;
        vmmc += nParticles;

        if (doReplace) {
            int nExited = checkAndReplace(model, particles, cells, box,
                                           nCopies, W, L_col, vmmc,
                                           targetComplexEnergy,
                                           totalExitedMass,
                                           totalExitedPerfect,
                                           fp_exits,
                                           globalStep);
            totalExited += nExited;
        }

        double energy      = model.getEnergy() * nParticles - baselineEnergy;
        double acceptRatio = (double)vmmc.getAccepts() / (double)vmmc.getAttempts();

        bool doSave = (globalStep % saveEvery == 0) || (globalStep == totalSteps);
        if (doSave) {
            writeFrame(fp_traj, particles, nCopies, L_col, W,
                       globalStep, energy, totalExited,
                       totalExitedMass, totalExitedPerfect, phase);
            fprintf(fp_stat, "%lld  %.4f  %lld  %lld  %lld  %.4f\n",
                    globalStep, energy, totalExited,
                    totalExitedMass, totalExitedPerfect, acceptRatio);
        }

        long long logEvery = max(1LL, totalSteps / 20);
        if (globalStep % logEvery == 0) {
            cout << "  [" << phase << "] step " << globalStep << "/" << totalSteps
                 << "  E=" << energy
                 << "  exited=" << totalExited
                 << "  perfect=" << (totalExitedPerfect / N0)
                 << "  accept=" << acceptRatio << "\n";
        }
    };

    // Phase 1: assembled free diffusion — gradient OFF (γ=1 everywhere, full bonding).
    // The purpose is to let the assembled complex diffuse freely before the gradient
    // is introduced; having the gradient on here would suppress bonds near x=0.
    if (freeSteps > 0) {
        cout << "Phase 1: assembled free diffusion (" << freeSteps << " steps, gradient off)..." << endl;
        model.hasGradient = false;
        model.denatured   = false;
        for (long long s = 0; s < freeSteps; s++)
            runStep("assembled", false);
        model.hasGradient = useGradient;  // restore for subsequent phases
    }

    // Phase 2: denaturation (β→0: γ=0 everywhere, all weak coupling zeroed; backbone intact).
    if (denatureSteps > 0) {
        cout << "Phase 2: denaturation (" << denatureSteps << " steps)..." << endl;
        model.denatured = true;
        for (long long s = 0; s < denatureSteps; s++)
            runStep("denature", true);
        model.denatured = false;  // restore coupling for phase 3
    }

    // Phase 3: main simulation — gradient on if --gradient was passed.
    if (nsteps > 0) {
        cout << "Phase 3: main simulation (" << nsteps << " steps)..." << endl;
        for (long long s = 0; s < nsteps; s++)
            runStep("main", true);
    }

    double simTime = (clock() - startTime) / (double)CLOCKS_PER_SEC;
    cout << "Done! Time = " << simTime << " s (" << simTime/60 << " min)" << endl;
    cout << "Total exited complexes: " << totalExited << endl;
    cout << "Total exited mass:      " << totalExitedMass << " particles" << endl;
    cout << "Perfect complex exits:  " << (totalExitedPerfect / N0)
         << " complexes (" << totalExitedPerfect << " particles)" << endl;
    cout << "Acceptance ratio: "
         << (double)vmmc.getAccepts() / (double)vmmc.getAttempts() << endl;

    fclose(fp_traj);
    fclose(fp_stat);
    fclose(fp_exits);

    cout << "Trajectory written to: " << trajFile << endl;
    cout << "Statistics written to:  " << statFile << endl;
    cout << "Exit events written to: " << exitsFile << endl;
    cout << "\nTo visualize:" << endl;
    cout << "  python3 visualize_nucleolus.py " << trajFile
         << " --gradient-length " << L_col << " --width " << W << endl;

    return EXIT_SUCCESS;
}
