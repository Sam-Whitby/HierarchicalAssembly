/* 

Created March 4, 2021 by Miranda Holmes-Cerfon

Builds on "SquareWellium" class made by Lester Hedges. 

Square Well interaction potential; 
different interactions for different pairs & different sides.

Assumes particles are on a lattice. Hence,
MUST SET ifLattice = true, for Initialise class & VMMC class. 

*/


#ifndef _STICKYSQUARE_H
#define _STICKYSQUARE_H

#include "Model.h"
#include "Box.h"

using namespace std;

// Triples, used to construct specific pair energies
struct Triple {
  int i;
  int j;
  double val;
};

// Holds indices & energies of interacting neighbours for a given particle
class Neighbours {
public: 
    vector<int> inds;      // indices of particles interact with
    vector<double> vals;   // energies of interactions
    void addVal(int,double);  // add ind/val pair
    double getVal(int);       // search for index; return corresponding value if found (INF if not)
    static constexpr double cNone = -99;  // value when interaction not found
};

// Holds interactions for a given particle
class Interactions {
public:
    vector<Neighbours> north;   // north neighbours
    vector<Neighbours> east;    // east neighbours
    double crosstalk = 0.0;     // crosstalk parameter; default = 0 (no crosstalk)
    int n0_size = 0;  // target structure size (for particle identity = index % n0_size)
    // Weak coupling matrices [id1][id2], symmetric, additive on top of backbone
    vector<vector<double>> weakD1;     // distance 1 (cardinal)
    vector<vector<double>> weakDsq2;   // distance sqrt(2) (diagonal)
    vector<vector<double>> weakD2;     // distance 2
    vector<vector<double>> weakDsq5;   // distance sqrt(5) (knight's move)
    // Backbone spring (unlimited range)
    double backboneSpringK = 0.0;            // spring stiffness k; 0 = disabled
    vector<vector<int>> backbonePartners;    // backbonePartners[id] = list of backbone partner IDs
    Interactions(int,vector<Triple>&,vector<Triple>&);   // constructor: nParticples, TriplesNorth, TriplesEast
    Interactions(int,vector<Triple>&,vector<Triple>&,double);   // constructor: nParticples, TriplesNorth, TriplesEast, crosstalk
    Interactions(int nParticles, int n0,
                 vector<Triple>& north, vector<Triple>& east,
                 vector<vector<double>>& wD1, vector<vector<double>>& wDsq2,
                 vector<vector<double>>& wD2,  vector<vector<double>>& wDsq5);
    Interactions(int nParticles, int n0,
                 vector<Triple>& north, vector<Triple>& east,
                 vector<vector<double>>& wD1, vector<vector<double>>& wDsq2,
                 vector<vector<double>>& wD2,  vector<vector<double>>& wDsq5,
                 double springK, vector<vector<int>>& bbPartners);

    void printInteractions(vector<Neighbours>&);
};


//! Class defining the square-well potential.
class StickySquare : public Model
{
public:

    // Data to keep track of specific pair interactions; size is nParticles
    Interactions interactions;

    //! Constructor.
    /*! \param box_ : A reference to the simulation box object.

        \param particles_ : A reference to the particle list.

        \param cells_ : A reference to the cell list object.

        \param maxInteractions_ : The maximum number of interactions per particle.

        \param interactionEnergy_ : The square well interaction energy (in units of kBT).

        \param interactionRange_ : The square well interaction range (in units of the particle diameter).

        \param Interactions : class containing north&east interactions 
     */
    StickySquare(Box&, std::vector<Particle>&, CellList&, unsigned int, double, double,
                 Interactions&);

    //! Calculate the pair energy between two particles.
    /*! \param particle1 : The index of the first particle.

        \param position1 : The position vector of the first particle.

        \param orientation1 : The orientation vector of the first particle.

        \param particle2 : The index of the second particle.

        \param position2 : The position vector of the second particle.

        \param orientation2 : The orientation vector of the second particle.

        \return : The pair energy between particles 1 and 2.
     */
    double computePairEnergy(unsigned int, const double*, const double*, unsigned int, const double*, const double*);

    // Same as above, but doesn't include crosstalk in pair interaction
    double computePairEnergyNative(unsigned int, const double*, const double*, unsigned int, const double*, const double*);

    unsigned int computeInteractions(unsigned int, const double*, const double*, unsigned int*);
    double computeEnergy(unsigned int, const double*, const double*);

    // Calculate the histogram of fragment sizes
    //   Input: 
    //         int maxFragmentSize       = maximum size of a fragment
    //         vector<int> fragmentHist  = histogram of fragment sizes (1:maxFragmentSize)
    //   Output: 
    //         int = # of fragments
    //
    int computeFragmentHistogram(int, vector<int>&);


    // Calculate the histogram of fragment sizes
    //   Input: 
    //         int maxFragmentSize       = maximum size of a fragment
    //         double emin               = minimum (absolute) energy to consider particles bound
    //         double emax               = maximum (absolute) energy to consider particles bound
    // 
    //         vector<int> fragmentHist  = histogram of fragment sizes (1:maxFragmentSize)
    //   Output: 
    //         int = # of fragments
    //
    int computeFragmentHistogramEnergy(int, double, double, vector<int>&);
    

    // Compute absolute value
    double mabs(double x) {
        if(x >= 0) return x;
        else return -x;
    }

};

#endif  /* _StickySquare_H */
