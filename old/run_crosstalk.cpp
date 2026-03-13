#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>             // for timing operations
#include <sstream>
#include <fstream>      // std::ifstream

#include "Demo.h"
#include "VMMC.h"
#include "StickySquare.h"

using namespace std;

//typedef std::chrono::high_resolution_clock Clock;

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif


int num(int i, int j, int l0) {
    return l0*j +i;
}


int main(int argc, char** argv)
{

    // Process input
    if(argc <= 3) {    
        cout << "Not enough input arguments! \nUsage: <program name> <inputfile> <e1> <ctalk> " << endl;
        return 1;
    }

    // parameters to set in command line
    string inputfile;
    double e1;
    double ctalk;


    // extract arguments from command line
    stringstream convert1 {argv[1]};
    stringstream convert2 {argv[2]};
    stringstream convert3 {argv[3]};
    convert1 >> inputfile;
    convert2 >> e1;
    convert3 >> ctalk;


    cout << "Checking inputs: " << endl;
    cout << "inputfile = " << inputfile << ", e1 = " << e1 << ", ctalk = " << ctalk << endl;
    cout << "Crosstalk: " << ctalk << endl;


    /* ----------  Parse input file for parameters  ---------- */

    string filehead;  // header for output files
    int n0;           // number of particles in base cell
    int nCopies;      // how many copies of each square there are
    int nsteps;       // # of steps of simulation
    int nsweep;       // # of monte carlo sweeps per step
    double dens;      // density
    
    ifstream paramfile;
    string line;
    stringstream stream1;
    paramfile.open(inputfile,ifstream::in);   // open file containing parameters
    if(paramfile.is_open()) {                 // successfully opened; now read in parameters
        getline(paramfile, line, ' ');
        filehead = line;
        getline(paramfile, line);  // remainder of line
        getline(paramfile, line, '#');
        stream1 << line;
        stream1 >> n0;
        stream1.clear();
        getline(paramfile, line);  // remainder of line
        getline(paramfile, line, '#');
        stream1 << line;
        stream1 >> nCopies;
        stream1.clear();
        getline(paramfile, line);  // remainder of line
        getline(paramfile, line, '#');
        stream1 << line;
        stream1 >> nsteps;
        stream1.clear();
        getline(paramfile, line);  // remainder of line
        getline(paramfile, line, '#');
        stream1 << line;
        stream1 >> nsweep;
        stream1.clear();
        getline(paramfile, line);  // remainder of line
        getline(paramfile, line, '#');
        stream1 << line;
        stream1 >> dens;
        stream1.clear();
        getline(paramfile, line);  // remainder of line

        cout << "filehead=" << filehead << ", n0=" << n0 << ", nCopies=" << nCopies << 
             ", nsteps=" << nsteps << ", nsweep="<<nsweep<<", dens="<<dens<<endl;

        paramfile.close();

    } else {
        cout << "Failed to open input file; exiting program" << endl;
        return 1;
    }



    /* ----------  Set additional parameters  ---------- */

    int l0 = round(sqrt(n0));    
    int f0 = n0;   // for reporting; what is completed fragment size?
    double boxLength = round(sqrt(n0*nCopies / dens));  // length of simulation box (must be an integer)
    int nParticles = n0*nCopies;       // number of particles


    // Set energy values
    double e2 = e1/2.0;
    double e3 = e1/4.0;
    double e4 = e1/8.0;
    double e5 = e1/16.0;

    // create filenames
    string statfile, trajfile;
    statfile = filehead + "_stats.txt";
    trajfile = filehead + "_traj.txt";
    //statfile = filehead + "_e" + to_string((int)(e1*100+0.5)) + "_c" + to_string((int)(ctalk*100+0.5)) + "_stats.txt";
    //trajfile = filehead + "_e" + to_string((int)(e1*100+0.5)) + "_c" + to_string((int)(ctalk*100+0.5)) + "_traj.txt";

    cout << "l0 = " << l0 << endl;
    cout << "statfile = " << statfile << endl;
    cout << "trajfilefile = " << trajfile << endl;



    /* ----------  Set up code  ---------- */

    // create string describing simulation parameters
    string description;
    stringstream os;
    os << "n0=" << n0 << " nCopies=" << nCopies << " nParticles=" << nParticles << " dens=" << dens
       << " nsteps=" << nsteps << " nsweep=" << nsweep << " ctalk = " << ctalk 
       << " e1=" << e1  << " e2=" << e2 << " e3=" << e3 << " e4 = " << e4 << " e5 = " << e5 ;
    description = os.str();

    // Display experiment info
    cout << "-----------\n  " << description << endl;
    cout << statfile << ", " << trajfile << endl;


    vector<double> stats;       // statistics to compute
    vector<int> fragmenthist;   // histogram of fragment sizes
    int nfrag;                  // number of fragments


    // Parameters that generally shouldn't be changed
    bool isLattice = true;              // whether particles must stay on a lattice
    unsigned int dimension = 2;         // dimension of simulation box
    double interactionRange = 1.1;      // interaction range (used in CellList; not used in StickySquare)
    unsigned int maxInteractions = 6;   // maximum number of interactions per particle (orig 15)
    double interactionEnergy = 0;       // (not used; needed to set up StickySphere) pair interaction energy scale (in units of kBT)

    // Initialise random number generator
    MersenneTwister rng;


    /* ----------  Make interactions  ---------- */

    vector<Triple> north0;
    vector<Triple> east0;
    vector<Triple> north;
    vector<Triple> east;
    double val;  // value of energy of interaction
    int num1,num2;  // indices of interaction

    auto num{ [l0](int i, int j) { return (l0*j + i); } };  // calculates square's actual number, from col/row

    // East interactions
    for(int i=0; i<l0-1; i++) {     // horizontal coordinate
        for(int j=0; j<l0; j++) {   // vertical coordinate
            if((i+1)%16 == 0 && n0>=1024) val = e5;
            else if((i+1)%8 == 0 && n0>=256) val = e4;
            else if((i+1)%4 == 0 && n0>=64) val = e3;
            else if((i+1)%2 == 0 && n0>=16) val = e2;
            else val = e1;
            // add this interaction
            num1 = num(i,j);
            num2 = num(i+1,j);
            east0.push_back({num1,num2,val});
        }
    }

    // North interactions
    for(int i=0; i<l0; i++) {     // horizontal coordinate
        for(int j=0; j<l0-1; j++) {   // vertical coordinate
            if((j+1)%16 == 0 && n0>=1024) val = e5;
            else if((j+1)%8 == 0 && n0>=256) val = e4;
            else if((j+1)%4 == 0 && n0>=64) val = e3;
            else if((j+1)%2 == 0 && n0>=16) val = e2;
            else val = e1;
            // add this interaction
            num1 = num(i,j);
            num2 = num(i,j+1);
            north0.push_back({num1,num2,val});
        }
    }


    // Done 1st square. Make several copies of each square
    for(int k=0; k<north0.size(); k++) {
        for(int c1=0; c1<nCopies; c1++) {
            for(int c2=0; c2<nCopies; c2++) {
                north.push_back({north0[k].i+c1*n0, north0[k].j+c2*n0, north0[k].val});            
            }   
        }
    }
    for(int k=0; k<east0.size(); k++) {
        for(int c1=0; c1<nCopies; c1++) {
            for(int c2=0; c2<nCopies; c2++) {
                east.push_back({east0[k].i+c1*n0, east0[k].j+c2*n0, east0[k].val});
            }
        }   
    }
    
    Interactions interactions(nParticles,north,east,ctalk);

    // debug
    //cout << "Crosstalk: " << interactions.crosstalk << endl;
    /*cout << "North:" << endl;
    interactions.printInteractions(interactions.north);
    cout << "East:" << endl;
    interactions.printInteractions(interactions.east);*/
    


    /* ----------  Initialise data structures & classes  ---------- */

    // Data structures.
    std::vector<Particle> particles(nParticles);    // particle container
    bool isIsotropic[nParticles];                   // whether the potential of each particle is isotropic
    
    // Create simulation box object.
    std::vector<double> boxSize {boxLength,boxLength};          // simulation box sizes
    Box box(boxSize,isLattice);  // Initialise simulation box object.

    // Initialise cell list.
    CellList cells; 
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    // Initialise the sticky square potential model.
    StickySquare StickySquare(box, particles, cells,
        maxInteractions, interactionEnergy, interactionRange,
        interactions);



    /* ----------  Initialise Particles  ---------- */
    // Generate a random particle configuration using Initialise object & MersenneTwister object
    Initialise initialise;
    initialise.random(particles, cells, box, rng, false, isLattice);

    // Initialise data structures needed by the VMMC class.
    double coordinates[dimension*nParticles];
    double orientations[dimension*nParticles];

    // Copy particle coordinates and orientations into C-style arrays.
    for (int i=0;i<nParticles;i++) {
        for (int j=0;j<dimension;j++) {
            coordinates[dimension*i + j] = particles[i].position[j];
            orientations[dimension*i + j] = particles[i].orientation[j];
        }
        // Set all particles as isotropic.
        isIsotropic[i] = true;
    }

    /* ----------  Initialise VMMC functions & object  ---------- */

    // Initialise the VMMC callback functions.
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


    // Variables to intialise VMMC object; these shouldn't change
    double maxTrialTranslation = 1.5;
    double maxTrialRotation = 0.0;
    double probTranslate = 1.0;
    double referenceRadius = 0.5;
    bool isRepulsive = false;

    // Initialise VMMC object. 
    vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
        maxTrialTranslation, maxTrialRotation, probTranslate, referenceRadius, 
        maxInteractions, &boxSize[0], isIsotropic, isRepulsive, callbacks, isLattice);


    /* ----------  Create output file  ---------- */

    // Create output file & log initial condition
    InputOutput io;
    //if(rep <= 1) {
        io.appendXyzTrajectory(dimension, particles, box, true, n0, description,trajfile);  // "true" is for the first line
    //}

    // Initalise statistics and write to file
    stats = {0,StickySquare.getEnergy()*nParticles};
    nfrag = StickySquare.computeFragmentHistogram(n0,fragmenthist);
    stats.insert(stats.end(), fragmenthist.begin(), fragmenthist.end());
    // write to file (erase old contents)
    io.appendStats(stats,true,description,statfile);   


    /* ----------  Run the simulation!  ---------- */
    clock_t start_time = clock();  // time the loop
    for (int i=0;i<nsteps;i++)
    {
        // Increment simulation by Monte Carlo Sweeps.
        vmmc += nsweep*nParticles; 

        // Append particle coordinates to an xyz trajectory.
        //if(rep <= 1) {
            io.appendXyzTrajectory(dimension, particles, false,trajfile);
        //}

        // Compute statistics
        stats = {(double)i, StickySquare.getEnergy()*nParticles};
        nfrag = StickySquare.computeFragmentHistogram(n0,fragmenthist);
        stats.insert(stats.end(), fragmenthist.begin(), fragmenthist.end());
        io.appendStats(stats,false,"",statfile);
    }
    // save time
    double time = (clock() - start_time ) / (double) CLOCKS_PER_SEC;

    /* ----------  Report stuff  ---------- */
    double enative;
    if(n0==4) {
        enative = -(4.0*e1)*(double)nCopies;
    } else if(n0==16) {
        enative = -(16.*e1+8.*e2)*(double)nCopies;
    } else if(n0==64) {
        enative = ((16.*e1+8.*e2)*4.+16.*e3)*(double)nCopies;
    } else if(n0==256) {
        enative = -(((16.*e1+8.*e2)*4.+16.*e3)*4 + 32*e4)*(double)nCopies;
    } else if(n0==1024) {
        enative = -((((16.*e1+8.*e2)*4.+16.*e3)*4 + 32*e4)*4 + 64*e5)*(double)nCopies;
    }
    double efinal = StickySquare.getEnergy()*nParticles;

    cout << "Complete! " << endl;
    cout << "  Time = " << time << " seconds, = " << time/60 << " minutes, = " << time/60/60 << " hours." << endl;
    cout << "  Acceptance ratio: " << (double)vmmc.getAccepts() / (double)vmmc.getAttempts() << endl;
    cout << "  Final energy             = " << efinal << endl;
    cout << "  Native structure energy  = " << enative << endl;
    cout << "  Ratio = " << efinal / enative << endl;
    cout << "  Number of fully completed fragments = " << stats.at(f0-1+2) << endl;
    cout << "    Fragments: ";
    for (int x : fragmenthist)  cout << x << " "; cout << endl;


    // We're done!
    return (EXIT_SUCCESS);
}
