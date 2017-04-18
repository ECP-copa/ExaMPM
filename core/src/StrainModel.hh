//---------------------------------------------------------------------------//
/*!
 * \file StrainModel.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_STRAINMODEL_HPP
#define EXAMPM_STRAINMODEL_HPP

#include "Particle.hh"

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*
 * \class Strain model.
 * \brief Strain model for a given material.
 */
class StrainModel
{
  public:

    // Destructor
    virtual ~StrainModel() = default;

    // Given a particle state calculate the particle strain.
    virtual void calculateStrain( ExaMPM::Particle& p ) const = 0;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_STRAINMODEL_HPP

//---------------------------------------------------------------------------//
// end StrainModel.hh
//---------------------------------------------------------------------------//
