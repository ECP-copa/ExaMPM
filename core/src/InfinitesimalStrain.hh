//---------------------------------------------------------------------------//
/*!
 * \file InfinitesimalStrain.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_INFINITESIMALSTRAIN_HH
#define EXAMPM_INFINITESIMALSTRAIN_HH

#include "StrainModel.hh"
#include "Particle.hh"

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*
 * \class InfinitesimalStrain
 * \brief Implements a small deformation theory strain.
 *
 * E = 0.5 * (F^T + F) - I
 */
class InfinitesimalStrain : public StrainModel
{
  public:

    // Given a particle state calculate the particle strain.
    void calculateStrain( ExaMPM::Particle& p ) const override;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_INFINITESIMALSTRAIN_HH

//---------------------------------------------------------------------------//
// end InfinitesimalStrain.hh
//---------------------------------------------------------------------------//
