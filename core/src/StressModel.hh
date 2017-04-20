//---------------------------------------------------------------------------//
/*!
 * \file StressModel.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_STRESSMODEL_HPP
#define EXAMPM_STRESSMODEL_HPP

#include "Particle.hh"

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*
 * \class Stress model.
 * \brief Specific stress model for a given material. The specific stress is
 * the stress divided by the density.
 */
class StressModel
{
  public:

    // Destructor
    virtual ~StressModel() = default;

    // Given a particle state calculate the specific particle stress.
    virtual void calculateStress( ExaMPM::Particle& p ) const = 0;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_STRESSMODEL_HPP

//---------------------------------------------------------------------------//
// end StressModel.hh
//---------------------------------------------------------------------------//
