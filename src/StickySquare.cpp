/* 

Created March 4, 2021 by Miranda Holmes-Cerfon

Builds on "SquareWellium" class made by Lester Hedges. 

Square Well interaction potential; 
different interactions for different pairs & different sides.

Assumes particles are on a lattice. Hence,
MUST SET ifLattice = true, for Initialise class & VMMC class. 

*/

#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "StickySquare.h"
#include <iostream>

extern double INF;


void Neighbours::addVal(int j,double val) {
    inds.push_back(j);
    vals.push_back(val);
}

double Neighbours::getVal(int k) {
    if(inds.size() == 0) return cNone;
    for(int j=0; j<inds.size(); j++) {
        if(inds[j] == k) {
            return vals[j];
        }
    }
    return cNone;
}


Interactions::Interactions(int np,vector<Triple>& triplesNorth, vector<Triple>& triplesEast) {
    // reserve space to hold interactions for all particles
    north.resize(np);
    east.resize(np);

    // build North interactions
    for(Triple t : triplesNorth) {
        north[t.i].addVal(t.j,t.val);
    }
    // build East interactions
    for(Triple t : triplesEast) {
        east[t.i].addVal(t.j,t.val);
    }
}

Interactions::Interactions(int np,vector<Triple>& triplesNorth, vector<Triple>& triplesEast, double pcrosstalk)
: crosstalk {pcrosstalk} {
    // reserve space to hold interactions for all particles
    north.resize(np);
    east.resize(np);

    // build North interactions
    for(Triple t : triplesNorth) {
        north[t.i].addVal(t.j,t.val);
    }
    // build East interactions
    for(Triple t : triplesEast) {
        east[t.i].addVal(t.j,t.val);
    }
}

Interactions::Interactions(int np, int n0,
        vector<Triple>& triplesNorth, vector<Triple>& triplesEast,
        vector<vector<double>>& wD1, vector<vector<double>>& wDsq2,
        vector<vector<double>>& wD2,  vector<vector<double>>& wDsq5)
    : Interactions(np, triplesNorth, triplesEast)
{
    n0_size  = n0;
    weakD1   = wD1;
    weakDsq2 = wDsq2;
    weakD2   = wD2;
    weakDsq5 = wDsq5;
}


// Print out list of interactions; for debugging
void Interactions::printInteractions(vector<Neighbours>& interactions) {
    int i=0;
    for(Neighbours nbr : interactions) {
        if(nbr.inds.size() == 0) 
            cout << i << "  NULL" << endl;
        else {
            for(int j=0; j< nbr.inds.size(); j++) 
                cout << i << "  " << nbr.inds[j] << ", val=" << nbr.vals[j] << endl;
        }
        /*cout << "Particle " << i << ", " << "  size = " << nbr.inds.size() << endl;
        for(int j=0; j< nbr.inds.size(); j++) 
            cout << "    " << nbr.inds[j] << ", val = " << nbr.vals[j] << endl;*/
        i++;
    }
}


// Constructor
StickySquare::StickySquare(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_, 
    Interactions& interactions_):
    Model(box_, particles_, cells_, maxInteractions_, interactionEnergy_, interactionRange_), 
    interactions(interactions_)
{
    // Check dimensionality is valid.
    if (box.dimension != 2)
    {
        std::cerr << "[ERROR] StickySquare: dimension must be 2!\n";
        exit(EXIT_FAILURE);
    }
}

// Compute pair interactions
double StickySquare::computePairEnergy(unsigned int particle1, const double* position1,
    const double* orientation1, unsigned int particle2, const double* position2, const double* orientation2)
{
    // Separation vector.
    std::vector<double> sep(box.dimension);

    // Calculate separation.
    for (unsigned int i=0;i<box.dimension;i++)
        sep[i] = position1[i] - position2[i];

    // Enforce minimum image.
    box.minimumImage(sep);

    // compute norm of separation
    double normSqd = 0;
    for (unsigned int i=0;i<box.dimension;i++)  normSqd += sep[i]*sep[i];
    
    // reject if particles overlap, or are beyond sqrt(5) range
    // TOL, INF defined in Model.cpp
    if (normSqd < 1.-TOL) return INF;   // particles overlap
    if (normSqd > 5.+TOL) return 0;     // beyond sqrt(5) range

    // Only works for 2d squares
    double energy = 0.0;
    if(box.dimension == 2) {
        int id1 = (int)particle1 % (interactions.n0_size > 0 ? interactions.n0_size : 1);
        int id2 = (int)particle2 % (interactions.n0_size > 0 ? interactions.n0_size : 1);

        if (normSqd < 1.+TOL) {
            // Cardinal neighbour (distance 1): directional backbone lookup
            double backbone = Neighbours::cNone;
            if( mabs(sep[0]-1.) < TOL && mabs(sep[1]) < TOL )
                backbone = interactions.east[particle2].getVal(particle1);
            if( mabs(sep[0]+1.) < TOL && mabs(sep[1]) < TOL )
                backbone = interactions.east[particle1].getVal(particle2);
            if( mabs(sep[1]-1.) < TOL && mabs(sep[0]) < TOL )
                backbone = interactions.north[particle2].getVal(particle1);
            if( mabs(sep[1]+1.) < TOL && mabs(sep[0]) < TOL )
                backbone = interactions.north[particle1].getVal(particle2);
            if (backbone != Neighbours::cNone) energy += backbone;
            else energy += interactions.crosstalk;
            // Weak coupling (always added)
            if (!interactions.weakD1.empty()) energy += interactions.weakD1[id1][id2];

        } else if (normSqd < 2.+TOL) {
            // Diagonal neighbour (distance sqrt(2)): direction-agnostic backbone lookup
            double backbone = interactions.east[particle1].getVal(particle2);
            if (backbone == Neighbours::cNone) backbone = interactions.east[particle2].getVal(particle1);
            if (backbone != Neighbours::cNone) energy += backbone;
            else energy += interactions.crosstalk;
            if (!interactions.weakDsq2.empty()) energy += interactions.weakDsq2[id1][id2];

        } else if (normSqd < 4.+TOL) {
            // Distance 2: weak coupling only
            if (!interactions.weakD2.empty()) energy += interactions.weakD2[id1][id2];

        } else {
            // Distance sqrt(5): weak coupling only
            if (!interactions.weakDsq5.empty()) energy += interactions.weakDsq5[id1][id2];
        }
    }
    return -energy;
}

// Compute pair interactions -- but ignore crosstalk (used for computing histograms)
// identical computePairEnergy, except doesn't add in crosstalk
double StickySquare::computePairEnergyNative(unsigned int particle1, const double* position1,
    const double* orientation1, unsigned int particle2, const double* position2, const double* orientation2)
{
    // Separation vector.
    std::vector<double> sep(box.dimension);

    // Calculate separation.
    for (unsigned int i=0;i<box.dimension;i++)
        sep[i] = position1[i] - position2[i];

    // Enforce minimum image.
    box.minimumImage(sep);

    // compute norm of separation
    double normSqd = 0;
    for (unsigned int i=0;i<box.dimension;i++)  normSqd += sep[i]*sep[i];

    // reject if particles overlap, or are beyond diagonal range
    // TOL, INF defined in Model.cpp
    if (normSqd < 1.-TOL) return INF;   // particles overlap
    if (normSqd > 2.+TOL) return 0;     // beyond cardinal + diagonal range

    // Only works for 2d squares
    double energy = Neighbours::cNone;
    if(box.dimension == 2) {
        if( normSqd < 1.+TOL ) {
            // Cardinal neighbour (distance 1): directional bond lookup
            if( mabs(sep[0]-1.) < TOL && mabs(sep[1]) < TOL )   // (1,0) = (2-->1 East)
                energy = interactions.east[particle2].getVal(particle1);
            if( mabs(sep[0]+1.) < TOL && mabs(sep[1]) < TOL )   // (-1,0) = (1-->2 East)
                energy = interactions.east[particle1].getVal(particle2);
            if( mabs(sep[1]-1.) < TOL && mabs(sep[0]) < TOL )   // (0,1) = (2-->1 North)
                energy = interactions.north[particle2].getVal(particle1);
            if( mabs(sep[1]+1.) < TOL && mabs(sep[0]) < TOL )   // (0,-1) = (1-->2 North)
                energy = interactions.north[particle1].getVal(particle2);
        } else {
            // Diagonal neighbour (distance sqrt(2)): direction-agnostic lookup
            energy = interactions.east[particle1].getVal(particle2);
            if(energy == Neighbours::cNone)
                energy = interactions.east[particle2].getVal(particle1);
        }
        if(energy == Neighbours::cNone) energy = 0;
    }
    //cout << "energy = " << energy << endl;
    return -energy;
}



// Calculate the histogram of fragment sizes
//   Input: 
//         int maxFragmentSize       = maximum size of a fragment
//         vector<int> fragmentHist  = histogram of fragment sizes (1:maxFragmentSize)
//   Output: 
//         int = # of fragments
//
int StickySquare::computeFragmentHistogram(int maxFragmentSize, vector<int>& fragmentHist) {

    int np = particles.size();  // number of particles
    double energy;

    // Reset fragment histogram, to all 0s
    fragmentHist.resize(maxFragmentSize);
    fill(fragmentHist.begin(), fragmentHist.end(), 0);   // sets every element to 0

    // keeps track of which fragment each particle is in. Count goes from 0 to (np-1). -1 means not yet assigned
    vector<int> fragmentID(np,-1);  // vector of size np, with elements initialized to -1

    // total number of fragments so far
    int nfrag = 0;

    // Loop through particles
    for(unsigned int i=0; i<np; i++) {  

        /* // debug
        cout << "i = " << i << ", ID = " << fragmentID.at(i) << endl;
        cout << "  fragmentHist: ";
        for (int x : fragmentHist)  cout << x << " ";
        cout << endl;
        cout << "  fragmentID: ";
        for (int x : fragmentID)  cout << x << " ";
        cout << endl;
        */

        // check if it's already in a fragment
        if(fragmentID.at(i) == -1) {  // it's not yet in a fragment

            vector<unsigned int> cluster {i};  // lists indices of particles in this fragment; initalise with i
            fragmentID.at(i) = nfrag;   // all particles in this fragment will have ID = nfrag

            // Loop through all particles in this cluster
            int icluster = 0;
            while(icluster < cluster.size()) {

                unsigned int j = cluster.at(icluster);  // the next particle in the cluster to consider

                // Get neighbours of j (see VMMC.cpp lines 514-556, and Model.cpp lines 111-167)
                unsigned int neighbours[maxInteractions];
                int nnbrs = computeInteractions(j, &particles[j].position[0], &particles[j].orientation[0], 
                                                neighbours);  // sets neighbours

                // Loop through neighbours
                for (int inbr=0; inbr < nnbrs; inbr++) {

                    unsigned int k = neighbours[inbr];   // index of neighbour particle

                    // check if it's in the cluster already
                    if(fragmentID.at(k) == -1)  {  // it's not yet in the cluster
                        // calculate energy between j & k
                        energy = computePairEnergyNative(j, &particles[j].position[0], &particles[j].orientation[0], 
                                                   k, &particles[k].position[0], &particles[k].orientation[0]);
                        if(energy < 0) {   // interaction is favourable
                            // add k to cluster, and record its ID
                            cluster.push_back(k);
                            fragmentID.at(k) = nfrag;
                        }
                    }
                }  // end loop through neighbours

                icluster++;   // move to next particle in cluster

            }   // end loop through cluster  

            nfrag++;  // update number of fragments
            fragmentHist.at(cluster.size()-1)++;   

        }  // end if particle not yet in a fragment

    }  // end loop through particles

    return nfrag;
}




// Calculate the histogram of fragment sizes, where fragments are bound with energy in some range
//   Input: 
//         int maxFragmentSize       = maximum size of a fragment
//         vector<int> fragmentHist  = histogram of fragment sizes (1:maxFragmentSize)
//   Output: 
//         int = # of fragments
//
int StickySquare::computeFragmentHistogramEnergy(int maxFragmentSize, double emin, double emax, vector<int>& fragmentHist) {

    double etol =  1e-5;  // tolerance for computing energies

    int np = particles.size();  // number of particles
    double energy;

    // Reset fragment histogram, to all 0s
    fragmentHist.resize(maxFragmentSize);
    fill(fragmentHist.begin(), fragmentHist.end(), 0);   // sets every element to 0

    // keeps track of which fragment each particle is in. Count goes from 0 to (np-1). -1 means not yet assigned
    vector<int> fragmentID(np,-1);  // vector of size np, with elements initialized to -1

    // total number of fragments so far
    int nfrag = 0;

    // Loop through particles
    for(unsigned int i=0; i<np; i++) {  

        /* // debug
        cout << "i = " << i << ", ID = " << fragmentID.at(i) << endl;
        cout << "  fragmentHist: ";
        for (int x : fragmentHist)  cout << x << " ";
        cout << endl;
        cout << "  fragmentID: ";
        for (int x : fragmentID)  cout << x << " ";
        cout << endl;
        */

        // check if it's already in a fragment
        if(fragmentID.at(i) == -1) {  // it's not yet in a fragment

            vector<unsigned int> cluster {i};  // lists indices of particles in this fragment; initalise with i
            fragmentID.at(i) = nfrag;   // all particles in this fragment will have ID = nfrag

            // Loop through all particles in this cluster
            int icluster = 0;
            while(icluster < cluster.size()) {

                unsigned int j = cluster.at(icluster);  // the next particle in the cluster to consider

                // Get neighbours of j (see VMMC.cpp lines 514-556, and Model.cpp lines 111-167)
                unsigned int neighbours[maxInteractions];
                int nnbrs = computeInteractions(j, &particles[j].position[0], &particles[j].orientation[0], 
                                                neighbours);  // sets neighbours

                // Loop through neighbours
                for (int inbr=0; inbr < nnbrs; inbr++) {

                    unsigned int k = neighbours[inbr];   // index of neighbour particle

                    // check if it's in the cluster already
                    if(fragmentID.at(k) == -1)  {  // it's not yet in the cluster
                        // calculate energy between j & k
                        energy = computePairEnergyNative(j, &particles[j].position[0], &particles[j].orientation[0], 
                                                   k, &particles[k].position[0], &particles[k].orientation[0]);
                        if(energy <= -(emin-etol) && energy >= -(emax+etol)) {   // interaction is favourable
                            // add k to cluster, and record its ID
                            cluster.push_back(k);
                            fragmentID.at(k) = nfrag;
                        }
                    }
                }  // end loop through neighbours

                icluster++;   // move to next particle in cluster

            }   // end loop through cluster  

            nfrag++;  // update number of fragments
            fragmentHist.at(cluster.size()-1)++;   

        }  // end if particle not yet in a fragment

    }  // end loop through particles

    return nfrag;
}
