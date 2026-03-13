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


int main(int argc, char** argv)
{
    // Process input
    if(argc <= 2) {
        cout << "Not enough input arguments!\n"
             << "Usage: ./run_custom <inputfile> <bondfile>\n"
             << "  inputfile : same format as run_hier (filehead, n, ncopies, nsteps, nsweep, dens)\n"
             << "  bondfile  : one bond per line:  particle_i  particle_j  energy\n"
             << "              Lines starting with '#' are comments.\n";
        return 1;
    }

    string inputfile, bondfile;
    stringstream convert1 {argv[1]};
    stringstream convert2 {argv[2]};
    convert1 >> inputfile;
    convert2 >> bondfile;

    cout << "inputfile = " << inputfile << ", bondfile = " << bondfile << endl;


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

    int l0 = round(sqrt(n0));
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
    double interactionRange = 1.5;   // covers cardinal (1) and diagonal (sqrt(2)) neighbours
    unsigned int maxInteractions = 8; // up to 4 cardinal + 4 diagonal per particle
    double interactionEnergy = 0;

    MersenneTwister rng;


    /* ----------  Read custom bond file  ---------- */
    //
    // Format: one bond per non-comment line
    //   particle_i  particle_j  energy
    //
    // particle_i and particle_j must be adjacent (Manhattan distance = 1) in
    // the n0-particle target grid  (col = p % l0,  row = p / l0).
    // Bonds are stored omni-directionally: each bond fires regardless of which
    // direction the two particles happen to be adjacent at runtime, enabling
    // floppy/flexible polymer conformations on the lattice.
    // Bonds with energy <= 0 are skipped.

    vector<Triple> north0, east0;

    ifstream bfile;
    bfile.open(bondfile, ifstream::in);
    if(!bfile.is_open()) {
        cout << "Failed to open bond file: " << bondfile << endl;
        return 1;
    }

    int bonds_read = 0;
    while(getline(bfile, line)) {
        if(line.empty() || line[0] == '#') continue;
        istringstream ss(line);
        int pi, pj;
        double val;
        if(!(ss >> pi >> pj >> val)) continue;
        if(pi < 0 || pi >= n0 || pj < 0 || pj >= n0) {
            cout << "Warning: bond (" << pi << "," << pj << ") out of range for n0="
                 << n0 << "; skipping." << endl;
            continue;
        }
        if(val <= 0.0) continue;

        int col_i = pi % l0,  row_i = pi / l0;
        int col_j = pj % l0,  row_j = pj / l0;
        int manhattan = abs(col_j - col_i) + abs(row_j - row_i);

        if(manhattan != 1) {
            cout << "Warning: bond (" << pi << "," << pj << ") is not adjacent in "
                 << "the " << l0 << "x" << l0 << " target grid; skipping." << endl;
        } else {
            // Store in all four directional table entries so the bond fires
            // regardless of which direction the particles are adjacent at runtime.
            east0.push_back({pi, pj, val});
            east0.push_back({pj, pi, val});
            north0.push_back({pi, pj, val});
            north0.push_back({pj, pi, val});
            bonds_read++;
        }
    }
    bfile.close();
    cout << "Loaded " << bonds_read << " omni-directional bonds." << endl;


    /* ----------  Duplicate interactions for nCopies  ---------- */

    vector<Triple> north, east;
    for(int k=0; k<(int)north0.size(); k++) {
        for(int c1=0; c1<nCopies; c1++) {
            for(int c2=0; c2<nCopies; c2++) {
                north.push_back({north0[k].i + c1*n0, north0[k].j + c2*n0, north0[k].val});
            }
        }
    }
    for(int k=0; k<(int)east0.size(); k++) {
        for(int c1=0; c1<nCopies; c1++) {
            for(int c2=0; c2<nCopies; c2++) {
                east.push_back({east0[k].i + c1*n0, east0[k].j + c2*n0, east0[k].val});
            }
        }
    }

    Interactions interactions(nParticles, north, east);


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

    Initialise initialise;
    initialise.random(particles, cells, box, rng, false, isLattice);

    double coordinates[dimension*nParticles];
    double orientations[dimension*nParticles];

    for(int i=0; i<nParticles; i++) {
        for(int j=0; j<dimension; j++) {
            coordinates[dimension*i + j] = particles[i].position[j];
            orientations[dimension*i + j] = particles[i].orientation[j];
        }
        isIsotropic[i] = true;
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
    double maxTrialRotation    = 0.0;
    double probTranslate       = 1.0;
    double referenceRadius     = 0.5;
    bool   isRepulsive         = false;

    vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
        maxTrialTranslation, maxTrialRotation, probTranslate, referenceRadius,
        maxInteractions, &boxSize[0], isIsotropic, isRepulsive, callbacks, isLattice);


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
    double time = (clock() - start_time) / (double)CLOCKS_PER_SEC;


    /* ----------  Report  ---------- */

    double efinal = StickySquare.getEnergy() * nParticles;
    cout << "Complete!" << endl;
    cout << "  Time = " << time << " s, " << time/60 << " min" << endl;
    cout << "  Acceptance ratio: " << (double)vmmc.getAccepts() / (double)vmmc.getAttempts() << endl;
    cout << "  Final energy = " << efinal << endl;
    cout << "  Number of fully completed fragments = " << stats.at(f0-1+2) << endl;
    cout << "  Fragments: ";
    for(int x : fragmenthist) cout << x << " ";
    cout << endl;

    return EXIT_SUCCESS;
}
