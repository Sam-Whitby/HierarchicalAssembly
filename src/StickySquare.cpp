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

Interactions::Interactions(int np, int n0,
        vector<Triple>& triplesNorth, vector<Triple>& triplesEast,
        vector<vector<double>>& wD1, vector<vector<double>>& wDsq2,
        vector<vector<double>>& wD2,  vector<vector<double>>& wDsq5,
        double springK, vector<vector<int>>& bbPartners)
    : Interactions(np, n0, triplesNorth, triplesEast, wD1, wDsq2, wD2, wDsq5)
{
    backboneSpringK  = springK;
    backbonePartners = bbPartners;
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


// Given a particle's orientation vector (ox, oy) and an integer cardinal
// direction (dx, dy) in the world frame, return which local patch slot
// of the particle faces that direction.
//   Slot 0: local-east  (same as orientation vector)
//   Slot 1: local-north (90° CCW from orientation)
//   Slot 2: local-west  (opposite to orientation)
//   Slot 3: local-south (90° CW from orientation)
// Returns -1 if (dx,dy) is not a valid cardinal direction.
static int getPatchSlot(const double* ori, int dx, int dy) {
    // Rotate world (dx,dy) into local frame: R^{-T} = R(-theta)
    // local_x = dx*ox + dy*oy,  local_y = -dx*oy + dy*ox
    int lx = (int)round((double)dx * ori[0] + (double)dy * ori[1]);
    int ly = (int)round(-(double)dx * ori[1] + (double)dy * ori[0]);
    if (lx ==  1 && ly ==  0) return 0;
    if (lx ==  0 && ly ==  1) return 1;
    if (lx == -1 && ly ==  0) return 2;
    if (lx ==  0 && ly == -1) return 3;
    return -1;  // not a cardinal direction
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

    // Backbone spring: applies at any distance (replaces wDsq2/wD2/wDsq5 for backbone pairs)
    if (interactions.backboneSpringK > 0.0 && !interactions.backbonePartners.empty()) {
        int n0 = interactions.n0_size > 0 ? interactions.n0_size : 1;
        int id1_ = (int)particle1 % n0;
        int id2_ = (int)particle2 % n0;
        bool isBackbone = false;
        for (int bp : interactions.backbonePartners[id1_]) {
            if (bp == id2_) { isBackbone = true; break; }
        }
        if (isBackbone) {
            double d = sqrt(normSqd);
            double physE = interactions.backboneSpringK * (d - 1.0) * (d - 1.0);
            // At d=1, also apply level coupling from wD1 (attractive term)
            if (normSqd < 1.0 + TOL && !interactions.weakD1.empty()) {
                physE -= interactions.weakD1[id1_][id2_];  // coupling → physical: -coupling
            }
            return physE;
        }
    }

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
            // Weak d=1 coupling, gated by directional patches when enabled.
            if (!interactions.weakD1.empty()) {
                bool patchOk = true;
                if (interactions.patchesEnabled && !interactions.patchSlots.empty()) {
                    // sep = pos1 - pos2, so direction from 1→2 is -sep, from 2→1 is +sep.
                    int dx12 = (int)round(-sep[0]);
                    int dy12 = (int)round(-sep[1]);
                    int s1 = getPatchSlot(orientation1,  dx12,  dy12);
                    int s2 = getPatchSlot(orientation2, -dx12, -dy12);
                    patchOk = (s1 >= 0 && s2 >= 0 &&
                               interactions.patchSlots[id1][s1] &&
                               interactions.patchSlots[id2][s2]);
                }
                if (patchOk) energy += interactions.weakD1[id1][id2];
            }

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
            // Patch gate on weakD1 (only active in patch mode; backbone excluded)
            if (energy == Neighbours::cNone && !interactions.weakD1.empty()) {
                double w = interactions.weakD1[
                    (int)particle1 % (interactions.n0_size > 0 ? interactions.n0_size : 1)][
                    (int)particle2 % (interactions.n0_size > 0 ? interactions.n0_size : 1)];
                if (interactions.patchesEnabled && !interactions.patchSlots.empty()) {
                    int id1n = (int)particle1 % (interactions.n0_size > 0 ? interactions.n0_size : 1);
                    int id2n = (int)particle2 % (interactions.n0_size > 0 ? interactions.n0_size : 1);
                    int dx12 = (int)round(-sep[0]);
                    int dy12 = (int)round(-sep[1]);
                    int s1 = getPatchSlot(orientation1,  dx12,  dy12);
                    int s2 = getPatchSlot(orientation2, -dx12, -dy12);
                    if (s1 >= 0 && s2 >= 0 &&
                        interactions.patchSlots[id1n][s1] &&
                        interactions.patchSlots[id2n][s2])
                        energy = w;
                } else {
                    energy = w;
                }
            }
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



// Override computeInteractions to always include backbone spring partners
unsigned int StickySquare::computeInteractions(unsigned int particle,
    const double* position, const double* orientation, unsigned int* nbrs)
{
    // Standard cell-list neighbors
    unsigned int n = Model::computeInteractions(particle, position, orientation, nbrs);

    // Add backbone partners not already found (regardless of current distance)
    if (interactions.backboneSpringK > 0.0 && !interactions.backbonePartners.empty()) {
        int n0 = interactions.n0_size > 0 ? interactions.n0_size : 1;
        int id1 = (int)particle % n0;
        int copyBase = (int)(particle / n0) * n0;
        for (int partner_id : interactions.backbonePartners[id1]) {
            unsigned int partner = (unsigned int)(copyBase + partner_id);
            bool found = false;
            for (unsigned int k = 0; k < n; k++) {
                if (nbrs[k] == partner) { found = true; break; }
            }
            if (!found) {
                if (n == maxInteractions) {
                    std::cerr << "[ERROR] StickySquare::computeInteractions: maxInteractions exceeded!\n";
                    exit(EXIT_FAILURE);
                }
                nbrs[n++] = partner;
            }
        }
    }
    return n;
}

// Override computeEnergy to use our computeInteractions (includes backbone partners at any range)
double StickySquare::computeEnergy(unsigned int particle, const double* position, const double* orientation)
{
    double energy = 0;
    unsigned int nbrs[maxInteractions];
    unsigned int n = computeInteractions(particle, position, orientation, nbrs);
    for (unsigned int k = 0; k < n; k++) {
        unsigned int nbr = nbrs[k];
        energy += computePairEnergy(particle, position, orientation,
                    nbr, &particles[nbr].position[0], &particles[nbr].orientation[0]);
    }
    return energy;
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
