/*
  NucleolusModel.h

  A spatial-gradient variant of the StickySquare model for simulating
  ribosome assembly within a condensate column (nucleolus model).

  Coupling between particles i and j is scaled by γ(x_i) * γ(x_j) where
      γ(x) = min(x / L, 1.0)
  and L is the condensate length.  The backbone bonds (strong intra-polymer
  cohesion, value ≈ 1000) are NOT scaled – they are always at full strength.
  Hard-core repulsion (d < 1) is also never scaled.

  With hasGradient = false every γ = 1 (uniform coupling, no gradient).
*/

#ifndef _NUCLEOLUS_MODEL_H
#define _NUCLEOLUS_MODEL_H

#include "StickySquare.h"

class NucleolusModel : public StickySquare
{
public:
    double columnLength;   //!< L: condensate length in lattice units.
    bool   hasGradient;    //!< Whether a spatial gradient is active.
    bool   denatured;      //!< When true γ=0 everywhere (β→0 denaturation phase).

    //! Constructor – identical to StickySquare plus L and gradient flag.
    NucleolusModel(Box&, std::vector<Particle>&, CellList&,
                   unsigned int maxInteractions, double interactionEnergy,
                   double interactionRange, Interactions&,
                   double L_, bool hasGradient_);

    //! γ(x) – chemical coupling parameter at position x.
    double gamma(double x) const;

    //! Override: pair energy scaled by γ(x_i) * γ(x_j) for weak terms.
    double computePairEnergy(unsigned int, const double*, const double*,
                             unsigned int, const double*, const double*);

    //! Override: total particle energy (calls our computePairEnergy).
    double computeEnergy(unsigned int, const double*, const double*);
};

#endif /* _NUCLEOLUS_MODEL_H */
