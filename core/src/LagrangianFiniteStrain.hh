//---------------------------------------------------------------------------//
/*!
 * \file LagrangianFiniteStrain.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_LAGRANGIANFINITESTRAIN_HH
#define EXAMPM_LAGRANGIANFINITESTRAIN_HH

#include "StrainModel.hh"
#include "Particle.hh"

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*
 * \class LagrangianFiniteStrain
 * \brief Implements a Lagrangian finite strain model:
 *
 * E = 0.5 * ( F*F^T - I )
 */
class LagrangianFiniteStrain : public StrainModel
{
  public:

    // Given a particle state calculate the particle strain.
    void calculateStrain( ExaMPM::Particle& p ) const override;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_LAGRANGIANFINITESTRAIN_HH

//---------------------------------------------------------------------------//
// end LagrangianFiniteStrain.hh
//---------------------------------------------------------------------------//
