#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <sstream>
#include <fstream>

#include "Demo.h"
#include "VMMC.h"
#include "StickySquare.h"

using namespace std;

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

// Parse one section of the extended bond file into a symmetric n0×n0 matrix.
// Returns an n0×n0 matrix filled from lines "i j energy" in the section.
// Call after reading the section header line; reads until END tag.
static vector<vector<double>> readMatrix(ifstream& f, int n0, const string& endTag) {
    vector<vector<double>> mat(n0, vector<double>(n0, 0.0));
    string line;
    while(getline(f, line)) {
        if(line.empty() || line[0] == '#') continue;
        if(line.find(endTag) != string::npos) break;
        istringstream ss(line);
        int pi, pj; double val;
        if(!(ss >> pi >> pj >> val)) continue;
        if(pi < 0 || pi >= n0 || pj < 0 || pj >= n0) continue;
        mat[pi][pj] = val;
        mat[pj][pi] = val;
    }
    return mat;
}

int main(int argc, char** argv)
{
    // Parse optional named flags first, then positional args.
    int nLatticeNeighbours = 4;
    double probTranslate   = 1.0;
    bool usePatch          = false;
    std::vector<string> positional;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        if (a == "--neighbours" && i + 1 < argc) {
            nLatticeNeighbours = atoi(argv[++i]);
            if (nLatticeNeighbours != 4 && nLatticeNeighbours != 8) {
                cerr << "[ERROR] --neighbours must be 4 or 8\n";
                return 1;
            }
        } else if (a == "--prob-translate" && i + 1 < argc) {
            probTranslate = atof(argv[++i]);
            if (probTranslate < 0.0 || probTranslate > 1.0) {
                cerr << "[ERROR] --prob-translate must be in [0, 1]\n";
                return 1;
            }
        } else if (a == "--patches") {
            usePatch = true;
        } else {
            positional.push_back(a);
        }
    }
    if (positional.size() < 2) {
        cout << "Usage: ./run_polymer <inputfile> <bondfile> [conffile] [seed] [--neighbours 4|8] [--prob-translate P] [--patches]\n"
             << "  inputfile        : same format as run_hier\n"
             << "  bondfile         : extended bond file with BACKBONE and WEAK_D* sections\n"
             << "  conffile         : (optional) initial configuration file (x y per line)\n"
             << "  seed             : (optional) integer RNG seed\n"
             << "  --neighbours 4|8 : lattice directions for VMMC (default 4)\n"
             << "  --prob-translate : fraction of moves that are translations (default 1.0)\n"
             << "  --patches        : enable directional patches; bonds only form when patch faces are aligned\n";
        return 1;
    }
    string inputfile = positional[0];
    string bondfile  = positional[1];
    string conffile  = (positional.size() >= 3) ? positional[2] : "";
    unsigned int rng_seed = 0;
    bool use_seed = false;
    if (positional.size() >= 4) { use_seed = true; rng_seed = (unsigned int)atoi(positional[3].c_str()); }

    cout << "inputfile=" << inputfile << " bondfile=" << bondfile;
    if(!conffile.empty()) cout << " conffile=" << conffile;
    if(use_seed) cout << " seed=" << rng_seed;
    cout << endl;

    /* ----------  Parse input file  ---------- */
    string filehead;
    int n0, nCopies, nsteps, nsweep;
    double dens;

    ifstream paramfile;
    string line;
    stringstream stream1;
    paramfile.open(inputfile, ifstream::in);
    if(paramfile.is_open()) {
        getline(paramfile, line, ' ');  filehead = line;
        getline(paramfile, line);
        getline(paramfile, line, '#');  stream1 << line; stream1 >> n0;      stream1.clear(); getline(paramfile, line);
        getline(paramfile, line, '#');  stream1 << line; stream1 >> nCopies; stream1.clear(); getline(paramfile, line);
        getline(paramfile, line, '#');  stream1 << line; stream1 >> nsteps;  stream1.clear(); getline(paramfile, line);
        getline(paramfile, line, '#');  stream1 << line; stream1 >> nsweep;  stream1.clear(); getline(paramfile, line);
        getline(paramfile, line, '#');  stream1 << line; stream1 >> dens;    stream1.clear(); getline(paramfile, line);
        paramfile.close();
        cout << "filehead=" << filehead << " n0=" << n0 << " nCopies=" << nCopies
             << " nsteps=" << nsteps << " nsweep=" << nsweep << " dens=" << dens << endl;
    } else {
        cout << "Failed to open input file; exiting" << endl;
        return 1;
    }

    /* ----------  Set additional parameters  ---------- */
    int l0 = (int)round(sqrt((double)n0));
    int f0 = n0;
    double boxLength = round(sqrt(n0 * nCopies / dens));
    int nParticles = n0 * nCopies;

    string statfile = filehead + "_stats.txt";
    string trajfile = filehead + "_traj.txt";

    string description;
    stringstream os;
    os << "n0=" << n0 << " nCopies=" << nCopies << " nParticles=" << nParticles
       << " dens=" << dens << " nsteps=" << nsteps << " nsweep=" << nsweep
       << " bondfile=" << bondfile;
    description = os.str();
    cout << "-----------\n  " << description << endl;

    vector<double> stats;
    vector<int> fragmenthist;
    int nfrag;

    bool isLattice = true;
    unsigned int dimension = 2;
    double interactionRange = 2.4;   // covers cardinal(1), diagonal(sqrt2), dist2, sqrt5
    unsigned int maxInteractions = 25; // up to 25 neighbors within sqrt(5) + backbone spring partners
    double interactionEnergy = 0;

    MersenneTwister rng;
    if(use_seed) rng.setSeed(rng_seed);

    /* ----------  Read extended bond file  ---------- */
    vector<Triple> north0, east0;
    vector<vector<double>> wD1(n0, vector<double>(n0, 0.0));
    vector<vector<double>> wDsq2(n0, vector<double>(n0, 0.0));
    vector<vector<double>> wD2(n0, vector<double>(n0, 0.0));
    vector<vector<double>> wDsq5(n0, vector<double>(n0, 0.0));
    bool has_weak = false;
    double springK = 0.0;
    vector<vector<int>> bbPartners(n0);

    ifstream bfile;
    bfile.open(bondfile, ifstream::in);
    if(!bfile.is_open()) {
        cout << "Failed to open bond file: " << bondfile << endl;
        return 1;
    }

    int bonds_read = 0;
    bool in_backbone = false;
    while(getline(bfile, line)) {
        if(line.empty() || line[0] == '#') continue;

        // Section headers
        if(line.find("SPRING_K_END") != string::npos) continue;
        if(line.find("SPRING_K") != string::npos && line.find("BACKBONE") == string::npos) {
            string kline;
            while(getline(bfile, kline)) {
                if(kline.empty() || kline[0] == '#') continue;
                if(kline.find("SPRING_K_END") != string::npos) break;
                istringstream ss2(kline); ss2 >> springK;
            }
            continue;
        }
        if(line.find("SPRING_BACKBONE_END") != string::npos) continue;
        if(line.find("SPRING_BACKBONE") != string::npos) {
            string bbline;
            while(getline(bfile, bbline)) {
                if(bbline.empty() || bbline[0] == '#') continue;
                if(bbline.find("SPRING_BACKBONE_END") != string::npos) break;
                istringstream ss3(bbline);
                int pi2, pj2;
                if(!(ss3 >> pi2 >> pj2)) continue;
                if(pi2 < 0 || pi2 >= n0 || pj2 < 0 || pj2 >= n0) continue;
                bbPartners[pi2].push_back(pj2);
                bbPartners[pj2].push_back(pi2);
            }
            continue;
        }
        if(line.find("BACKBONE_END") != string::npos) { in_backbone = false; continue; }
        if(line.find("BACKBONE") != string::npos)     { in_backbone = true;  continue; }
        if(line.find("WEAK_D1_END") != string::npos || line.find("WEAK_DSQRT2_END") != string::npos ||
           line.find("WEAK_D2_END") != string::npos || line.find("WEAK_DSQRT5_END") != string::npos) continue;
        if(line.find("WEAK_D1") != string::npos && line.find("SQRT") == string::npos) {
            wD1 = readMatrix(bfile, n0, "WEAK_D1_END"); has_weak = true; continue; }
        if(line.find("WEAK_DSQRT2") != string::npos) {
            wDsq2 = readMatrix(bfile, n0, "WEAK_DSQRT2_END"); has_weak = true; continue; }
        if(line.find("WEAK_D2") != string::npos && line.find("SQRT") == string::npos) {
            wD2 = readMatrix(bfile, n0, "WEAK_D2_END"); has_weak = true; continue; }
        if(line.find("WEAK_DSQRT5") != string::npos) {
            wDsq5 = readMatrix(bfile, n0, "WEAK_DSQRT5_END"); has_weak = true; continue; }

        if(!in_backbone) continue;

        // Parse backbone bond: particle_i  particle_j  energy
        istringstream ss(line);
        int pi, pj; double val;
        if(!(ss >> pi >> pj >> val)) continue;
        if(pi < 0 || pi >= n0 || pj < 0 || pj >= n0) {
            cout << "Warning: bond (" << pi << "," << pj << ") out of range; skipping.\n";
            continue;
        }
        if(val <= 0.0) continue;

        int col_i = pi % l0, row_i = pi / l0;
        int col_j = pj % l0, row_j = pj / l0;
        int manhattan = abs(col_j - col_i) + abs(row_j - row_i);
        if(manhattan != 1) {
            cout << "Warning: bond (" << pi << "," << pj << ") not adjacent in grid; skipping.\n";
        } else {
            east0.push_back({pi, pj, val});
            east0.push_back({pj, pi, val});
            north0.push_back({pi, pj, val});
            north0.push_back({pj, pi, val});
            bonds_read++;
        }
    }
    bfile.close();
    cout << "Loaded " << bonds_read << " backbone bonds." << endl;
    if(has_weak) cout << "Weak coupling matrices loaded." << endl;

    /* ----------  Duplicate interactions for nCopies  ---------- */
    vector<Triple> north, east;
    for(int k=0; k<(int)north0.size(); k++)
        for(int c1=0; c1<nCopies; c1++)
            for(int c2=0; c2<nCopies; c2++)
                north.push_back({north0[k].i + c1*n0, north0[k].j + c2*n0, north0[k].val});
    for(int k=0; k<(int)east0.size(); k++)
        for(int c1=0; c1<nCopies; c1++)
            for(int c2=0; c2<nCopies; c2++)
                east.push_back({east0[k].i + c1*n0, east0[k].j + c2*n0, east0[k].val});

    // Duplicate weak matrices for nCopies (same identity matrix, just indexed by mod n0)
    // (n0_size handles the modulo, so no duplication needed for the matrices themselves)

    Interactions interactions(nParticles, n0, north, east, wD1, wDsq2, wD2, wDsq5, springK, bbPartners);

    /* ----------  Build directional patches (--patches mode)  ---------- */
    // Each particle identity i has 4 patch slots in its LOCAL frame:
    //   slot 0 = local-east (faces same dir as orientation vector)
    //   slot 1 = local-north (90° CCW from orientation)
    //   slot 2 = local-west  (opposite orientation)
    //   slot 3 = local-south (90° CW from orientation)
    // A slot is activated when wD1[i][j] != 0 and j is at grid-distance 1 from i.
    // All particles start with orientation (1,0), so local frame == world frame at t=0:
    //   direction (+1, 0) → slot 0,  (0,+1) → slot 1,  (-1,0) → slot 2,  (0,-1) → slot 3.
    if (usePatch) {
        interactions.patchesEnabled = true;
        interactions.patchSlots.assign(n0, {false, false, false, false});
        for (int i = 0; i < n0; i++) {
            int col_i = i % l0, row_i = i / l0;
            for (int j = 0; j < n0; j++) {
                if (wD1[i][j] == 0.0) continue;
                int col_j = j % l0, row_j = j / l0;
                int dx = col_j - col_i, dy = row_j - row_i;
                // Only grid-adjacent (Manhattan distance 1) pairs get patches.
                if (abs(dx) + abs(dy) != 1) continue;
                // Backbone pairs bypass the patch check anyway; exclude their
                // directions so backbone-direction slots can't accidentally gate
                // non-backbone bonds after a rotation.
                bool isBackbone = false;
                for (int bp : bbPartners[i]) {
                    if (bp == j) { isBackbone = true; break; }
                }
                if (isBackbone) continue;
                int slot = (dx == 1) ? 0 : (dy == 1) ? 1 : (dx == -1) ? 2 : 3;
                interactions.patchSlots[i][slot] = true;
            }
        }
        int total_patches = 0;
        for (int i = 0; i < n0; i++)
            for (int s = 0; s < 4; s++)
                if (interactions.patchSlots[i][s]) total_patches++;
        cout << "Patch mode enabled: " << total_patches << " patch slots activated across " << n0 << " particle identities." << endl;
    }

    /* ----------  Initialise data structures & classes  ---------- */
    std::vector<Particle> particles(nParticles);
    bool isIsotropic[nParticles];

    std::vector<double> boxSize {boxLength, boxLength};
    Box box(boxSize, isLattice);

    CellList cells;
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    StickySquare StickySquare(box, particles, cells,
        maxInteractions, interactionEnergy, interactionRange, interactions);

    /* ----------  Initialise Particles  ---------- */
    if(!conffile.empty()) {
        // Load initial configuration from file
        InputOutput io_init;
        io_init.loadConfiguration(conffile, box, particles, cells, true);
        cout << "Loaded initial configuration from " << conffile << endl;
    } else {
        Initialise initialise;
        initialise.random(particles, cells, box, rng, false, isLattice);
    }

    double coordinates[dimension*nParticles];
    double orientations[dimension*nParticles];
    for(int i=0; i<nParticles; i++) {
        for(int j=0; j<(int)dimension; j++) {
            coordinates[dimension*i + j] = particles[i].position[j];
            orientations[dimension*i + j] = particles[i].orientation[j];
        }
        // Patch mode: particles are anisotropic (orientations tracked by VMMC).
        // All start facing (1,0) so their local patch frame matches the Hilbert grid layout.
        if (usePatch) {
            orientations[dimension*i + 0] = 1.0;
            orientations[dimension*i + 1] = 0.0;
            particles[i].orientation[0]   = 1.0;
            particles[i].orientation[1]   = 0.0;
        }
        isIsotropic[i] = !usePatch;
    }

    /* ----------  Initialise VMMC  ---------- */
    using namespace std::placeholders;
    vmmc::CallbackFunctions callbacks;
    callbacks.energyCallback =
        std::bind(&StickySquare::computeEnergy, StickySquare, _1, _2, _3);
    callbacks.pairEnergyCallback =
        std::bind(&StickySquare::computePairEnergy, StickySquare, _1, _2, _3, _4, _5, _6);
    callbacks.interactionsCallback =
        std::bind(&StickySquare::computeInteractions, StickySquare, _1, _2, _3, _4);
    callbacks.postMoveCallback =
        std::bind(&StickySquare::applyPostMoveUpdates, StickySquare, _1, _2, _3);

    double maxTrialTranslation = 1.5;
    double maxTrialRotation    = (probTranslate < 1.0) ? M_PI : 0.0;
    double referenceRadius     = 0.5;
    bool   isRepulsive         = false;

    vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
        maxTrialTranslation, maxTrialRotation, probTranslate, referenceRadius,
        maxInteractions, &boxSize[0], isIsotropic, isRepulsive, callbacks, isLattice,
        nLatticeNeighbours);

    /* ----------  Create output files  ---------- */
    InputOutput io;
    io.appendXyzTrajectory(dimension, particles, box, true, n0, description, trajfile);

    stats = {0, StickySquare.getEnergy() * nParticles};
    nfrag = StickySquare.computeFragmentHistogram(n0, fragmenthist);
    stats.insert(stats.end(), fragmenthist.begin(), fragmenthist.end());
    io.appendStats(stats, true, description, statfile);

    /* ----------  Run the simulation  ---------- */
    clock_t start_time = clock();
    for(int i=0; i<nsteps; i++) {
        vmmc += nsweep * nParticles;
        io.appendXyzTrajectory(dimension, particles, false, trajfile);
        stats = {(double)i, StickySquare.getEnergy() * nParticles};
        nfrag = StickySquare.computeFragmentHistogram(n0, fragmenthist);
        stats.insert(stats.end(), fragmenthist.begin(), fragmenthist.end());
        io.appendStats(stats, false, "", statfile);
    }
    double time_s = (clock() - start_time) / (double)CLOCKS_PER_SEC;

    /* ----------  Report  ---------- */
    double efinal = StickySquare.getEnergy() * nParticles;
    cout << "Complete!" << endl;
    cout << "  Time = " << time_s << " s, " << time_s/60 << " min" << endl;
    cout << "  Acceptance ratio: " << (double)vmmc.getAccepts() / (double)vmmc.getAttempts() << endl;
    cout << "  Final energy = " << efinal << endl;
    cout << "  Number of fully completed fragments = " << stats.at(f0-1+2) << endl;
    cout << "  Fragments: ";
    for(int x : fragmenthist) cout << x << " ";
    cout << endl;

    return EXIT_SUCCESS;
}
