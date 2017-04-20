//---------------------------------------------------------------------------//
/*!
 * \file EulerianAlmansiFiniteStrain.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_EULERIANALMANSIFINITESTRAIN_HH
#define EXAMPM_EULERIANALMANSIFINITESTRAIN_HH

#include "StrainModel.hh"
#include "Particle.hh"

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*
 * \class EulerianAlmansiFiniteStrain
 * \brief Implements a Eulerian-Almansi finite strain model:
 *
 * E = 0.5 * ( F*F^T - I )
 */
class EulerianAlmansiFiniteStrain : public StrainModel
{
  public:

    // Given a particle state calculate the particle strain.
    void calculateStrain( ExaMPM::Particle& p ) const override;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_EULERIANALMANSIFINITESTRAIN_HH

//---------------------------------------------------------------------------//
// end EulerianAlmansiFiniteStrain.hh
//---------------------------------------------------------------------------//
